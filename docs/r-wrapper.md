# R Wrapper

The R package currently exposes primitive matrix helpers for the implemented
early normalization stages:

- `estimateSizeFactorsRust(counts, method = c("ratio", "poscounts"), geoMeans = NULL, controlGenes = NULL, native = FALSE)`
- `normalizedCountsRust(counts, sizeFactors = NULL, normalizationFactors = NULL, native = FALSE)`
- `baseMeanRust(counts, sizeFactors = NULL, normalizationFactors = NULL, native = FALSE)`
- `baseMetadataRust(counts, sizeFactors = NULL, normalizationFactors = NULL, weights = NULL, native = FALSE)`
- `applyCooksCutoffRust(pvalue, maxCooks, cooksCutoff, counts = NULL, cooks = NULL, lowCountHeuristic = FALSE, native = FALSE)`
- `resultsTableRust(baseMean, log2FoldChange = NULL, lfcSE = NULL, stat = NULL, pvalue = NULL, padj = NULL, dispersion = NULL, converged = NULL, rowNames = NULL)`
- `applyIndependentFilteringRust(results, alpha = 0.1, theta = NULL, enabled = TRUE)`
- `rsdeseq2DiagnosticSchemaRust(native = FALSE)`
- `deseq2McolsDiagnosticsRust(nGenes, test = c("none", "Wald", "LRT"), ...)`

These helpers validate ordinary R matrices with genes in rows and samples in
columns. They use an R fallback implementation that mirrors the Rust-supported
size-factor, normalized-count, normalization-factor, baseMean, and early
base-metadata algorithms while the native Rust bridge is still being wired.
When `normalizationFactors` are supplied to `normalizedCountsRust()`,
`baseMeanRust()`, or `baseMetadataRust()`, they preempt `sizeFactors`, matching
DESeq2's normalized-count behavior.

`estimateSizeFactorsRust(native = TRUE)`, `normalizedCountsRust(native = TRUE)`,
and `baseMeanRust(native = TRUE)` first try registered `.Call` bridges for
primitive normalization summaries and fall back to the R implementation if the
package shared library is not loaded.

`baseMetadataRust()` returns primitive row metadata columns `baseMean`,
`baseVar`, and `allZero`. When `weights` are supplied, they must be a finite
non-negative matrix with the same dimensions as `counts`; the helper multiplies
normalized counts by those raw observation weights before calculating row means
and sample variances. This mirrors the early DESeq2 weighted metadata shape but
does not modify DESeqDataSet `mcols()`. With `native = TRUE`, the helper first
tries the registered `.Call` bridge for this primitive and falls back to the R
implementation if the package shared library is not loaded.

`applyCooksCutoffRust()` is a primitive result-table helper. It masks p-values
where `maxCooks > cooksCutoff`, recomputes BH-adjusted p-values, and can apply
the DESeq2 two-group low-count Cook's heuristic when the caller supplies raw
`counts`, per-sample `cooks`, and has already established the one-factor
two-level formula condition. It does not inspect DESeqDataSet objects or infer
formula semantics. With `native = TRUE`, the masking and low-count heuristic
step can use the registered `.Call` bridge; R still assembles the output table
and computes BH-adjusted p-values.

`resultsTableRust()` assembles already-computed primitive vectors into a
DESeq2-shaped result data frame with `baseMean`, `log2FoldChange`, `lfcSE`,
`stat`, `pvalue`, and `padj`, plus optional `dispersion` and `converged`
diagnostic columns. It computes BH-adjusted p-values from `pvalue` when `padj`
is omitted. It does not run Wald/LRT fitting, contrast construction, Cook's
masking, outlier replacement, or independent filtering.

`applyIndependentFilteringRust()` applies baseMean-driven filtered BH
adjustment to a primitive result data frame. It uses the same default theta-grid
shape as the Rust core, candidate cutoffs from R type-7 quantiles, BH adjustment
after filtering, and `stats::lowess(..., f = 1/5)` metadata for threshold
selection. The returned data frame gets a `filtered` column and an
`independentFiltering` attribute with theta, rejection counts, selected
threshold, lowess fit, and alpha. This helper does not inspect DESeqDataSet
objects or integrate with high-level `results()`.

`deseq2McolsDiagnosticsRust()` is a pure-R shape contract for the Rust fit-state
diagnostic alias view. It validates primitive vectors and returns a data frame
with DESeq2-style row metadata names such as `betaConv`, `fullBetaConv`,
`reducedBetaConv`, `betaIter`, `reducedBetaIter`, `deviance`, and `maxCooks`.
`rsdeseq2DiagnosticSchemaRust()` returns the same schema and can call the
registered native C stub when the shared library is loaded. These helpers do
not inspect or modify DESeqDataSet objects.

Unsupported high-level entry points remain explicit:

- `DESeq()`
- `results()`
- `estimateDispersionsRust()`
- `nbinomWaldTestRust()`
- `nbinomLRTRust()`

Full DESeqDataSet integration is not implemented yet.
