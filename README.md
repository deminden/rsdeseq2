# rsdeseq2

[![Rust CI](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml)
[![R CI](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml)
[![crates.io](https://img.shields.io/crates/v/rsdeseq2.svg)](https://crates.io/crates/rsdeseq2)
[![API documentation](https://img.shields.io/docsrs/rsdeseq2)](https://docs.rs/rsdeseq2)
[![License: MIT + BSD notices](https://img.shields.io/badge/license-MIT%20%2B%20BSD-blue.svg)](THIRD_PARTY_NOTICES.md)

`rsdeseq2` independently implements DESeq2-compatible count-data
workflows in Rust: normalization, dispersion estimation, negative-binomial
Wald and likelihood-ratio tests, transformations, and result assembly. The
Rust crate and CLI require neither R nor DESeq2 at runtime.

The Rust core uses DESeq2 only for validation references. `rsdeseq2` is not
a drop-in replacement for `DESeqDataSet` or the complete Bioconductor API.

## Choose an Interface

| interface | intended use | scope |
| --- | --- | --- |
| [CLI](docs/cli.md) | File-based normalization, Wald/LRT, VST, and rlog workflows | Explicit numeric TSV inputs; no R formula parser |
| [Rust crate](crates/rsdeseq2/README.md) | Embedding analysis stages in Rust applications | Broadest implemented API |
| [R package](r-pkg/rsdeseq2/README.md) | Rust-backed matrix, list, formula, and limited `DESeqDataSet` workflows; VST/rlog and selected shrinkage | Does not implement the complete Bioconductor object API |

Version 0.5.0 runs supported formula, transform and shrinkage
calculations in the shared Rust core. The [supported contract](docs/native-supported-contract.md)
records methods and limits; [CI validation](docs/testing-parity.md#ci-status)
records passing Rust, R, platform, full-GTEx and sanitizer checks.

## Validation and Performance

The recorded process comparison uses native R 4.6.2 / DESeq2 1.52.0 on eleven
GTEx v8 Kidney Cortex contrasts: biological comparisons and reproducible
blocked/unblocked random splits. Every contrast uses all 56,200 genes.

### Precision

Maximum absolute differences from matched DESeq2 results, including estimated
dispersions. Group values are maxima across distinct contrasts; gene coverage
and finite/NA patterns match throughout.

| Contrast group | Contrasts | Max LFC error | Max SE error | Max p-value error | Max adjusted p-value error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original blocked | 1 | `7.84e-07` | `1.38e-06` | `3.73e-07` | `1.14e-06` |
| Male vs female | 1 | `1.41e-09` | `1.35e-05` | `2.98e-05` | `3.42e-11` |
| Male age 4 vs 3 | 1 | `1.47e-07` | `2.34e-05` | `2.16e-05` | `8.25e-07` |
| Blocked random splits | 5 | `7.15e-06` | `2.75e-05` | `3.74e-05` | `1.26e-05` |
| Unblocked random splits | 3 | `5.28e-09` | `1.2e-05` | `1.72e-05` | `2.87e-05` |

### Speed

| Contrast group | Rust (s) | DESeq2 (s) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 5.96 | 31.03 | 5.20× |
| Male vs female | 2.75 | 17.69 | 6.42× |
| Male age 4 vs 3 | 2.38 | 14.54 | 6.10× |
| Blocked random splits | 4.39 | 23.71 | 5.40× |
| Unblocked random splits | 3.43 | 19.63 | 5.72× |

One process run per tool and contrast on Apple M5/macOS, with one fitting/BLAS
thread. Multiple-contrast groups show medians across distinct contrasts. Timings
include startup, loading, fitting and completed output. Summed time is
**43.33 seconds for Rust versus 242.12 seconds for DESeq2**, an observed **5.59×** ratio.
These are single-pass observations, without discarded outliers. Other user and
system workloads were active; the measurements do not establish small speed
changes between development revisions.

### Memory

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 304.62 | 1571.23 | 5.16× |
| Male vs female | 243.05 | 1689.53 | 6.95× |
| Male age 4 vs 3 | 213.61 | 1600.56 | 7.49× |
| Blocked random splits | 304.41 | 1681.58 | 5.52× |
| Unblocked random splits | 244.09 | 1643.09 | 6.73× |

Peak whole-process RSS from the same runs; MiB = 2^20 bytes, including the R
interpreter and libraries. Groups show medians across distinct contrasts.
Rust peaks at **214–307 MiB**, using **80.6–86.7% less peak RSS than DESeq2**.

Per-contrast measurements, provenance, precision checks and reproducible scripts
are in the [current comparison](docs/reports/gtex-native-wald-benchmark.md).
Rust and R tests also cover documented normalization, Wald/LRT, transforms,
shrinkage, fitted-object reuse and batched results. This establishes behavior for
the tested inputs and settings; the complete Bioconductor API is outside scope.

## Installation

The library and CLI are distributed through [crates.io](https://crates.io/crates/rsdeseq2/0.5.0).
Install the published CLI with `cargo install rsdeseq2 --version 0.5.0 --locked`,
or add the library with `cargo add rsdeseq2@0.5.0`. Prebuilt CLI archives are
attached to the [versioned GitHub release](https://github.com/deminden/rsdeseq2/releases/tag/v0.5.0).

For source installation, check out `main`:

```bash
git clone --branch main https://github.com/deminden/rsdeseq2.git
cd rsdeseq2
```

Source installation requires **Rust 1.97.1 or newer, with Cargo**. Installing
the R package also compiles Rust and requires R and its native build tools.
R 4.6.2 is the pinned comparison environment; see the
[R guide](docs/r-wrapper.md) for supported interfaces and optional dependencies.

| Interface | Install from the checkout root | Usage |
| --- | --- | --- |
| CLI | `cargo install --path crates/rsdeseq2 --locked` | [Commands and input schemas](docs/cli.md) |
| Rust library | `cargo add rsdeseq2 --path /absolute/path/to/rsdeseq2/crates/rsdeseq2` from your application | [Library usage](crates/rsdeseq2/README.md#rust-usage) |
| R package | `R CMD INSTALL r-pkg/rsdeseq2` | [Fitting and results](docs/r-wrapper.md) |

For a Wald analysis with your own count and numeric design tables:

```bash
rsdeseq2 wald \
  --counts counts.tsv \
  --design design.tsv \
  --fit-type mean \
  --coefficient 1 \
  --output results.tsv
```

Here `--fit-type mean` explicitly selects a constant mean dispersion trend;
it is a modeling choice. Select the trend appropriate to your analysis.
Counts must be unnormalized non-negative integers. The leading columns are
`gene` in the count table and `sample` in the design table; sample labels align
the two tables. Coefficient indices are zero-based. The [CLI guide](docs/cli.md)
describes complete schemas, contrasts and outputs.

## Supported Workflows

| Area | Implemented scope |
| --- | --- |
| Normalization and fitting | Ratio/`poscounts`, normalization factors, dispersion estimation, Wald and LRT |
| Results and diagnostics | Contrasts, Cook's replacement/refitting, independent filtering and adjusted p-values |
| Transformations and shrinkage | VST, rlog, and supported normal, `ashr` and `apeglm` methods |
| R integration | Matrix/list and supported formula/`DESeqDataSet` workflows; selected Bioconductor result objects |

Full Bioconductor object semantics, `glmGamPoi`, plotting and custom formula
functions are outside scope. Method-specific restrictions are maintained in
[compatibility and limits](docs/compatibility.md) and the
[supported contract](docs/native-supported-contract.md).
The public Rust API is pre-1.0 and can change between minor releases.

## Documentation

- **Run analyses:** [CLI inputs and commands](docs/cli.md)
- **Check feature support:** [compatibility and limits](docs/compatibility.md)
- **Understand numerical behavior:** [algorithms](docs/algorithms.md)
- **Inspect evidence:** [benchmarks](docs/benchmarks.md) and
  [reproducibility](docs/reproducibility.md)
- **Use the R layer:** [R wrapper status](docs/r-wrapper.md)
- **Contribute or release:** [development](docs/development.md) and
  [release checklist](docs/release.md)

## Support and Citation

Report bugs and parity differences through
[GitHub Issues](https://github.com/deminden/rsdeseq2/issues). Include the
rsdeseq2 version, command or API call, input dimensions, and the R/DESeq2
versions used for comparison.

For the statistical method, cite Love, Huber, and Anders (2014),
[“Moderated estimation of fold change and dispersion for RNA-seq data with
DESeq2”](https://doi.org/10.1186/s13059-014-0550-8). Cite the exact rsdeseq2
release used for software provenance.

## License

Project code is MIT-licensed; adapted numerical components also carry
BSD-3-Clause and BSD-2-Clause-Views terms. See [LICENSE](LICENSE) and
[third-party notices](THIRD_PARTY_NOTICES.md), consistent with the Cargo
license expression `MIT AND BSD-3-Clause AND BSD-2-Clause-Views`.

The R-compatible arithmetic is independently implemented from published
probability identities and validated against saved black-box R outputs. No R
or DESeq2 implementation is copied, linked, or distributed.
