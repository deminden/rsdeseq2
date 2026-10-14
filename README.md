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

DESeq2 is used only to generate saved validation references. `rsdeseq2` is not
a drop-in replacement for `DESeqDataSet` or the complete Bioconductor API.

## Choose an Interface

| interface | intended use | scope |
| --- | --- | --- |
| [CLI](docs/cli.md) | File-based normalization, Wald/LRT, VST, and rlog workflows | Explicit numeric TSV inputs; no R formula parser |
| [Rust crate](crates/rsdeseq2/README.md) | Embedding analysis stages in Rust applications | Broadest implemented API |
| [R package](r-pkg/rsdeseq2/README.md) | Rust-backed matrix, list, formula, and limited `DESeqDataSet` workflows; VST/rlog and selected shrinkage | Does not implement the complete Bioconductor object API |

## Install and Run

Install the released CLI from crates.io:

```bash
cargo install rsdeseq2 --locked
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
| Transformations and priors | Supported for documented paths | Mean, parametric, and local VST; rlog; normal, `ashr`, and `apeglm` shrinkage; selected formula/model-frame workflows. Formula interactions and some transformed-factor contrast forms are unsupported |
| R package | Partial | Matrix, list, and limited `DESeqDataSet` Wald/LRT fits; dataset transforms return Bioconductor transform objects and selected `lfcShrink` methods return `DESeqResults`. Full diagnostic and refit semantics are not implemented |
| Complete DESeq2/Bioconductor surface | Not supported | No `glmGamPoi`, plotting, arbitrary R expressions, or complete object-mutation and metadata semantics |

The public Rust API is pre-1.0 and can change between minor releases. See the
[compatibility matrix](docs/compatibility.md) and
[DESeq2 gap analysis](docs/deseq2-gap-analysis.md) for feature-level detail.

## Validation

Automated tests compare Rust stages with saved DESeq2 references for
normalization, dispersion trends, GLM fitting, Wald/LRT results, and
transformations. The real-data Wald test runs all 56,200 GTEx v8 genes across
61 Kidney Cortex samples and compares Rust outputs with DESeq2 1.52.0. For the
48,683 rows with finite Wald results, maximum absolute differences are
`2.386e-7` for log2 fold change, `4.140e-10` for standard error, `3.734e-7`
for the statistic, and `4.955e-8` for the p-value; size factors differ by at
most `5.551e-15` relative. R integration tests also check `DESeqDataSet` fits,
seeded full/fast VST and rlog workflows, and supported `lfcShrink()` methods.
These tests establish behavior for their documented inputs, not complete
DESeq2 workflow parity. See the [real-data fixture notes](crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md)
for provenance and reproduction.

See [Benchmarks](docs/benchmarks.md), [Compatibility](docs/compatibility.md),
and [Reproducibility](docs/reproducibility.md) for fixture definitions,
commands, and detailed results.

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

The workspace requires Rust 1.97.1 or newer and uses the Rust 2024 edition, as
declared by `rust-version` and `edition` in `Cargo.toml`. The R package is tested
with R 4.6.1. Each saved fixture records its generating environment.

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
