# rsdeseq2 R Package

This package provides R access to the `rsdeseq2` Rust core. It builds and
statically links the core from source, so Rust and Cargo are required to install
the package. It supports matrix and formula-driven Wald/LRT workflows and
selected DESeq2-shaped helpers, but not every Bioconductor object operation.

Version 0.2.6 adds dataset-based VST/rlog and normal, `ashr`, and `apeglm`
shrinkage workflows on top of the Rust-backed fit. See the sections below for
the supported inputs and known boundaries.

The R package is tested with R 4.6.1. Numerical reference fixtures record the
R/DESeq2 versions that generated them.

## Install and Run an Analysis

From a repository checkout:

```bash
R CMD INSTALL r-pkg/rsdeseq2
```

```r
library(rsdeseq2)

counts <- matrix(
  c(10L, 12L, 20L, 24L,
    5L,  7L,  6L,  8L,
    100L, 80L, 90L, 120L),
  nrow = 3L,
  byrow = TRUE
)
design <- cbind(Intercept = 1, condition = c(0, 0, 1, 1))
fit <- DESeq(list(counts = counts, design = design), fitType = "mean")
results(fit, independentFiltering = FALSE)
```

`DESeq()` also accepts a numeric model matrix through `design`, and a
Bioconductor `DESeqDataSet` when `SummarizedExperiment` is installed. Formula
designs are expanded by R and fitted by Rust. LRT calls require a reduced design
matrix or formula. The returned `rsdeseq2RustFit` stores coefficient estimates,
covariance, size factors, the Rust result table, count matrix, and design
metadata. `results()` supports named coefficients and DESeq2-style contrasts.
Existing dataset size factors, normalization factors, observation-weight
assays, and formula offsets are forwarded into the Rust fit.
When input is a `DESeqDataSet`, `DESeq()` returns that dataset with fit state
attached in metadata; matrix and list inputs return `rsdeseq2RustFit`.
`estimateDispersionsRust()` returns the final per-gene dispersion estimates
from the same Rust workflow.

For `DESeqDataSet` inputs, `vst()` and
`varianceStabilizingTransformation()` return `DESeqTransform` objects, while
`rlog()` also supports reuse of supplied frozen gene intercepts and a shared
sample-effect prior variance. These transformations run in Rust and honor
design/blind mode, normalization factors, observation weights, and formula
offsets. `lfcShrink()` supports `type="normal"` for a named or indexed Wald
coefficient or an expanded-design factor-level, numeric, or list contrast from
a dataset fitted with `rsdeseq2::DESeq()`; the contrast route uses the Rust
expanded-model refit. Numeric contrasts require a supplied result table; list
contrasts average selected expanded coefficients and need supplied results
when their names are absent from the unshrunk design. The returned
`DESeqResults` records the estimated beta-prior variances. With optional
package `ashr` installed, `type="ashr"` returns posterior means/SDs and
s-values. With optional package `apeglm` installed, `type="apeglm"` supports
named/indexed coefficients and returns posterior SDs and s-values. Formula
interactions are rejected for normal shrinkage. Cleaned coefficient aliases
work for transformed-factor formulas; factor-level contrast syntax remains
unsupported when DESeq2 cannot resolve its coefficient names.

## Supported Helpers

The package also exposes `estimateSizeFactorsRust()`, `normalizedCountsRust()`,
`baseMeanRust()`, `baseMetadataRust()`, `applyCooksCutoffRust()`,
`resultsTableRust()`, `waldFitRust()`, `lrtFitRust()`, `nbinomWaldTestRust()`,
`nbinomLRTRust()`, `resultsNames()`, `resolveResultsContrastRust()`, and
`applyIndependentFilteringRust()`.

These helpers call Rust through registered `.Call` bridges. Legacy `native`
arguments remain accepted for source compatibility; they no longer select an R
implementation. Full DESeq2 diagnostic assays and every DESeq2 fitting option
remain outside this package's implemented scope.

The Rust core validation summary and benchmark scope are documented in the
[repository README](../../README.md#validation) and
[benchmark report](../../docs/benchmarks.md).
