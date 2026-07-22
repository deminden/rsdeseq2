# rsdeseq2

[![Rust CI](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/rust.yml)
[![R CI](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml/badge.svg)](https://github.com/deminden/rsdeseq2/actions/workflows/r.yml)
[![crates.io](https://img.shields.io/crates/v/rsdeseq2.svg)](https://crates.io/crates/rsdeseq2)
[![API documentation](https://img.shields.io/docsrs/rsdeseq2)](https://docs.rs/rsdeseq2)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`rsdeseq2` independently implements selected DESeq2-compatible count-data
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
| [R package](r-pkg/rsdeseq2/README.md) | Selected primitive helpers and native bridges | Does not fit complete `DESeqDataSet` workflows |

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
| Transformations and priors | Partly supported | Beta-prior refits, supported formula/model-frame workflows, `normTransform`, VST, and rlog building blocks |
| R package | Limited | Primitive helpers and selected native bridges; no complete `DESeqDataSet` fitting |
| Complete DESeq2/Bioconductor surface | Not supported | No `glmGamPoi`, `lfcShrink`, plotting, arbitrary R expressions, or full object-mutation and metadata semantics |

The public Rust API is pre-1.0 and can change between minor releases. See the
[compatibility matrix](docs/compatibility.md) and
[DESeq2 gap analysis](docs/deseq2-gap-analysis.md) for feature-level detail.

## Measured Validation

Measurements recorded on 2026-07-22 compare a release-mode rsdeseq2 0.2.5
binary with saved R 4.6.1 / DESeq2 1.52.0 outputs unless stated otherwise.

- **Broad real-data Wald checks:** four full-blocked contrasts across three
  tissues covered 278,257 result rows with zero missing-row or finite/NA-pattern
  mismatches. Maximum absolute Wald-statistic errors by contrast were
  `1.526e-3`, `1.768e-3`, `5.615e-5`, and `2.803e-4`.
- **Frozen high-error regression set:** for 100 rows selected by the largest
  rsdeseq2 0.2.4 errors, the measured 0.2.5 median, mean, and maximum selected
  result-column absolute errors were `6.063612945084174e-10`,
  `1.260818774570247e-4`, and `1.5259081158007781e-3`. The 0.2.4 measurements
  were `1.4637657972313423e-4`, `3.793917566690452e-4`, and
  `3.0938714191082184e-3`: the 0.2.5 values were 241401.589x, 3.00909032x, and
  50.6796531% lower. Of the 100 rows, 89 improved and 78 improved by at least
  10x. This selected set is a regression stress test, not a distribution-wide
  accuracy estimate.
- **Normalization sweep:** 17 tissues and 8,731 samples covered 612,699,575
  normalized-count values with zero finite/NA mismatches. The measured maximum
  absolute difference was `1.937e-7`, and the maximum relative difference was
  `9.887e-15`.
- **Optimizer isolation:** replaying 512 bounded negative-binomial objectives
  with `rcompat-lbfgsb` 0.2.1 produced 512/512 exact endpoints, objective
  values, and evaluation counts. The 0.1.6 measurements were 0/512 exact
  endpoint-plus-objective matches and 311/512 exact evaluation counts.
- **Primitive process benchmark:** for matched size-factor and base-mean CLI
  stages over 10,000/50,000 genes and 16 samples, rsdeseq2 measured
  `0.010–0.140 s` and `6.01–16.54 MiB`; DESeq2 measured `3.865–4.625 s` and
  `601.64–638.58 MiB`. The measured ratios were 33.04x–386.5x lower elapsed
  time and 38.62x–100.05x lower peak RSS.

The 0.2.5 precision work did not establish a release-to-release performance
gain. In the controlled two-run heart replay, the measured elapsed medians for
size factors, normalized counts, and base mean were `2.108 s`, `14.317 s`, and
`2.107 s`; the 0.2.4 measurements were `2.090 s`, `14.477 s`, and `2.099 s`.
The relative changes were +0.87%, -1.11%, and +0.39%, within the observed
run-to-run variation. The larger speed and memory ratios above apply only to
the explicitly matched primitive stages, not to a complete `DESeq()` workflow.

Ratio size-factor estimation uses compensated accumulation of log counts
before downstream fitting. Rare optimizer fallback rows use independently
implemented R-compatible arithmetic and a reference-independent stability
check. Full measurements, interpretation, fixtures, and commands are in
[Benchmarks](docs/benchmarks.md), [Algorithms](docs/algorithms.md), and
[Reproducibility](docs/reproducibility.md).

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
