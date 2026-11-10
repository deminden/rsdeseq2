# rsdeseq2 R Package

This package provides R access to the `rsdeseq2` Rust core. It builds and
statically links the core from source, so Rust and Cargo are required to install
the package. It supports matrix and formula-driven Wald/LRT workflows and
selected DESeq2-shaped helpers, but not every Bioconductor object operation.

For `DESeqDataSet` inputs, the fit attaches Rust `mu` and original `cooks`
assays, fitted row metadata, and per-coefficient MLE/SE columns. It performs
Cook's outlier replacement and refits eligible rows by default when the design
has samples with at least seven replicates. Set `minReplicatesForReplace = Inf`
to disable replacement, or lower the threshold (minimum 3) for smaller designs.
Original counts and Cook's distances are retained. When genes are replaced,
the dataset also gets `replaceCounts` and `replaceCooks` assays.
`results()` metadata includes the replacement threshold, eligible samples, and
refit rows.

The R package is tested with R 4.6.2 / DESeq2 1.52.0. Numerical reference fixtures record the
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
design metadata is parsed in R; supported formula expressions, frames, offsets,
and model matrices are calculated in the shared Rust core. Custom formula
functions are rejected before execution. LRT calls require a reduced design
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
sample-effect prior variance. The main transformation fits run in Rust and honor
design/blind mode, normalization factors, observation weights, and formula
offsets. `lfcShrink()` supports `type="normal"` for a named or indexed Wald
coefficient or an expanded-design factor-level, numeric, or list contrast from
a dataset fitted with `rsdeseq2::DESeq()`; the contrast route uses the Rust
expanded-model refit. Numeric contrasts require a supplied result table; list
contrasts average selected expanded coefficients and need supplied results
when their names are absent from the unshrunk design. The returned
`DESeqResults` records the estimated beta-prior variances. Local development
`type="ashr"` returns native posterior means/SDs and s-values, while
`type="apeglm"` uses native negative-binomial fitting for
named/indexed coefficients and returns posterior SDs and s-values. Formula
interactions are rejected for normal shrinkage. Cleaned coefficient aliases
work for transformed-factor formulas; factor-level contrast syntax remains
unsupported when DESeq2 cannot resolve its coefficient names.

The default optional shrinkage workflows, general apeglm likelihoods,
`nbinomR` and beta-binomial `betabinR`/`betabinCR`/`betabinC` run in Rust.
General likelihood functions are translated from
restricted expressions without R execution. MAP-only `nbinomC*` and `betabinC*`
use native random starts with portable R state import/export. Advanced method options
remain unfinished and fail explicitly; see the
[migration plan](../../docs/rust-only-wrapper-migration.md). Ordered factors use
polynomial coordinates; select their named coefficients or supply numeric
contrasts rather than treatment-style factor-level contrasts.

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
[repository README](../../README.md#validation-and-performance) and
[benchmark report](../../docs/benchmarks.md).

Local development also supports native apeglm grid sign summaries, scalar
parameter uncertainty and covariance-based posterior contrasts. The primitive
engine exposes interval/contrast matrices; `lfcShrink()` keeps its selected-
coefficient result-table shape. See the
[posterior analysis scope](../../docs/native-posterior-analysis-development.md)
for controls, validation and documented original-package grid corrections.
