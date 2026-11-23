# rsdeseq2

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

The shared Rust core runs supported formula, transform and shrinkage
calculations. The [supported contract](docs/native-supported-contract.md)
records methods and limits. The [0.6.0 release notes](docs/releases/0.6.0.md)
describe the changes and current validation status.

## Validation and Performance

The 0.6.0 release panels and golden outputs passed on macOS arm64 and
Linux arm64. The full 56,200-gene GTEx replay passed in both numeric
profiles. Both profiles also passed all 130 precision gates across the
eleven full GTEx v8 Kidney Cortex Wald contrasts against R 4.6.2 / DESeq2
1.52.0. The contrast runs estimate normalization and dispersions with Cook's
masking and independent filtering disabled. Both profiles had zero
significance flips at adjusted p-value cutoffs of 0.05 and 0.1.

| Numeric profile | Max absolute LFC error | Max absolute SE error | Max absolute p-value error | Max absolute adjusted p-value error |
| --- | ---: | ---: | ---: | ---: |
| Portable | `3.65e-6` | `6.18e-5` | `7.25e-5` | `3.57e-5` |
| Host R | `7.15e-6` | `2.75e-5` | `3.74e-5` | `2.87e-5` |

A single four-worker macOS arm64 run of the portable CLI took 15.85 s
across the eleven contrasts and peaked at 259.41 MiB RSS.
[Validation and measurement details](docs/reports/hardening-progress.md#batched-acceptance-after-the-portable-arithmetic-changes)
record the checked revision, comparison settings and output checks.

## Installation

This checkout contains the **0.6.0** source. The published version is **0.5.0**
on [crates.io](https://crates.io/crates/rsdeseq2/0.5.0) and
[GitHub Releases](https://github.com/deminden/rsdeseq2/releases/tag/v0.5.0).
Use the checkout for the 0.6.0 interfaces described here.

Source installation requires **Rust 1.96.0 or newer, with Cargo**. Installing
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
| Transformations and shrinkage | VST, rlog and normal shrinkage; `ashr` and `apeglm` methods are **experimental** (see below) |
| R integration | Matrix/list and `DESeqDataSet` workflows; formula evaluation is **experimental**; selected Bioconductor result objects |

**Experimental:** `ashr`- and `apeglm`-style shrinkage and R-style formula
evaluation have independent reference tests. The formula parser and its
transport also have fuzz targets. Stability will require broader property
testing and validation on independent datasets.

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
