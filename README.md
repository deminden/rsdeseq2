# rsdeseq2

[![Rust CI](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml)
[![R CI](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml)
[![crates.io](https://img.shields.io/crates/v/rsdeseq2.svg)](https://crates.io/crates/rsdeseq2)
[![API documentation](https://img.shields.io/docsrs/rsdeseq2)](https://docs.rs/rsdeseq2)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

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

The R wrapper’s [Rust-only migration](docs/rust-only-wrapper-migration.md) is
in progress. Supported formula calculations use the shared Rust core; frozen
local VST reuse and optional ashr/apeglm posterior calculations still involve R.

## Install and Run

Install the 0.4.0 CLI from crates.io:

```bash
cargo install rsdeseq2 --version 0.4.0 --locked
```

From a repository checkout, run the complete small Wald example:

```bash
cargo run --release -p rsdeseq2 -- wald \
  --counts examples/quickstart/counts.tsv \
  --design examples/quickstart/design.tsv \
  --fit-type mean \
  --coefficient 1 \
  --output results.tsv

sed -n '1,5p' results.tsv
```

The result table contains `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
`pvalue`, `padj`, fitted dispersion, convergence, and filtering fields.

The CLI expects unnormalized non-negative integer counts. `counts.tsv` has a
leading `gene` column and one column per sample; `design.tsv` has a leading
`sample` column followed by numeric design-matrix columns. Design rows are
aligned to count columns by sample label. See the [CLI reference](docs/cli.md)
for complete schemas, contrasts, normalization inputs, and output sidecars.

### Rust library

Add the library to a Rust project:

```bash
cargo add rsdeseq2
```

A normalization primitive:

```rust
use rsdeseq2::prelude::*;

fn main() -> Result<(), DeseqError> {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 12, 20, 24,
            0, 0, 5, 7,
            100, 80, 90, 120,
        ],
    )?;

    let fit = DeseqBuilder::new()
        .size_factor_method(SizeFactorMethod::Ratio)
        .execution_mode(ExecutionMode::Strict)
        .fit_size_factors_and_base_means(&counts)?;

    println!("{:?}", fit.size_factors);
    Ok(())
}
```

### R access layer

Install the source package from a repository checkout:

```bash
R CMD INSTALL r-pkg/rsdeseq2
```

Then call a supported primitive:

```r
library(rsdeseq2)

counts <- matrix(
  c(10L, 12L, 20L, 24L,
    5L,  7L,  6L,  8L,
    100L, 80L, 90L, 120L),
  nrow = 3L,
  byrow = TRUE
)
estimateSizeFactorsRust(counts, native = TRUE)
```

## Supported Workflows

| need | status | notes |
| --- | --- | --- |
| Ratio and `poscounts` normalization | Supported | Includes supplied geometric means, control genes, normalized counts, base metadata, and normalization-factor offsets |
| Wald and LRT result workflows | Supported for documented paths | Supplied-dispersion GLMs and implemented native-dispersion branches; coefficient, list, and numeric contrasts; Normal and t tails |
| Diagnostics and result assembly | Supported for documented paths | Cook's diagnostics and replacement/refitting, independent filtering, adjusted p-values, and result metadata |
| Transformations and priors | Supported for documented paths | Mean, parametric, and local VST; rlog; normal, `ashr`, and `apeglm` shrinkage; selected formula/model-frame workflows. Normal shrinkage with formula interactions and some transformed-factor contrast forms are unsupported |
| R package | Partial | Matrix, list, and limited `DESeqDataSet` Wald/LRT fits with fitted assays, row/result metadata, and eligible Cook's replacement/refits; dataset transforms return Bioconductor transform objects and selected `lfcShrink` methods return `DESeqResults`. Additional Bioconductor semantics remain unsupported |
| Complete DESeq2/Bioconductor surface | Not supported | No `glmGamPoi`, plotting, arbitrary R expressions, or complete object-mutation and metadata semantics |

The public Rust API is pre-1.0 and can change between minor releases. See the
[compatibility matrix](docs/compatibility.md) and
[DESeq2 gap analysis](docs/deseq2-gap-analysis.md) for feature-level detail.

## Validation and Performance

Version 0.4.0 is compared with native R 4.6.2 / DESeq2 1.52.0 on eleven
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
| Original blocked | 3.75 | 19.08 | 5.09× |
| Male vs female | 2.58 | 14.49 | 5.61× |
| Male age 4 vs 3 | 2.23 | 13.46 | 6.03× |
| Blocked random splits | 3.90 | 19.60 | 5.02× |
| Unblocked random splits | 2.54 | 14.71 | 5.78× |

One process run per tool and contrast on Apple M5/macOS, with one fitting/BLAS
thread. Multiple-contrast groups show medians across distinct contrasts. Timings
include startup, loading, fitting and completed output. Summed time is
**36.27 seconds for Rust versus 189.25 seconds for DESeq2**, an observed **5.22×** ratio.
These are single-pass observations, without discarded outliers.

### Memory

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 302.06 | 1763.86 | 5.84× |
| Male vs female | 240.94 | 1689.47 | 7.01× |
| Male age 4 vs 3 | 213.02 | 1600.69 | 7.51× |
| Blocked random splits | 301.84 | 1757.19 | 5.82× |
| Unblocked random splits | 240.89 | 1632.77 | 6.78× |

Peak whole-process RSS from the same runs; MiB = 2^20 bytes, including the R
interpreter and libraries. Groups show medians across distinct contrasts.
Rust peaks at **213–302 MiB**, using **82.8–86.7% less peak RSS than DESeq2**.

Per-contrast measurements, provenance, precision checks and reproducible scripts
are in the [current comparison](docs/gtex-native-wald-benchmark.md).
Rust and R tests also cover documented normalization, Wald/LRT, transforms,
shrinkage, fitted-object reuse and batched results. This establishes behavior for
the tested inputs and settings; the complete Bioconductor API is outside scope.

## Documentation

- **Run analyses:** [CLI inputs and commands](docs/cli.md)
- **Check feature support:** [compatibility matrix](docs/compatibility.md) and
  [DESeq2 gap analysis](docs/deseq2-gap-analysis.md)
- **Understand numerical behavior:** [algorithms](docs/algorithms.md)
- **Inspect evidence:** [benchmarks](docs/benchmarks.md) and
  [reproducibility](docs/reproducibility.md)
- **Use the R layer:** [R wrapper status](docs/r-wrapper.md)
- **Contribute or release:** [development](docs/development.md) and
  [release checklist](docs/release.md)

## Development

Python tooling/tests require Python 3.10 or newer; the validation run used Python 3.14.7.

The workspace requires Rust 1.97.1 or newer and uses the Rust 2024 edition, as
declared by `rust-version` and `edition` in `Cargo.toml`. The R package is tested
with native R 4.6.2. Each saved fixture records its generating environment.

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
```

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

MIT. The R-compatible arithmetic is independently implemented from published
probability identities and validated against saved black-box R outputs. No R
or DESeq2 implementation is copied, linked, or distributed.
