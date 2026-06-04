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
columns. They include an R-level helper implementation that mirrors the
Rust-supported size-factor, normalized-count, normalization-factor, baseMean,
and early base-metadata algorithms while the native Rust bridge is still being
wired.
When `normalizationFactors` are supplied to `normalizedCountsRust()`,
`baseMeanRust()`, or `baseMetadataRust()`, they preempt `sizeFactors`, matching
DESeq2's normalized-count behavior.

`estimateSizeFactorsRust(native = TRUE)`, `normalizedCountsRust(native = TRUE)`,
and `baseMeanRust(native = TRUE)` first try registered `.Call` bridges for
primitive normalization summaries. If the package shared library is not loaded,
the helper uses the R-level implementation.

`baseMetadataRust()` returns primitive row metadata columns `baseMean`,
`baseVar`, and `allZero`. When `weights` are supplied, they must be a finite
non-negative matrix with the same dimensions as `counts`; the helper multiplies
normalized counts by those raw observation weights before calculating row means
and sample variances. This mirrors the early DESeq2 weighted metadata shape but
does not modify DESeqDataSet `mcols()`. With `native = TRUE`, the helper first
tries the registered `.Call` bridge for this primitive. If the package shared
library is not loaded, the helper uses the R-level implementation.

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
is omitted. It does not run Wald/LRT fitting, Cook's masking, outlier
replacement, or independent filtering.

`waldFitRust()` packages already-computed primitive Wald beta matrices,
coefficient covariance arrays, optional coefficient SEs, and optional contrast
diagnostics such as counts, sample levels, factor levels, and model matrices.
`results()` can consume this primitive fit object for default/name-selected
coefficients and DESeq2-style `contrast = ...` forms, including
`lfcThreshold`/`altHypothesis` Wald p-values with Normal or supplied t tails.
For convenience, `results()` can also consume the same list-like primitive Wald
fields directly and internally package them through `nbinomWaldTestRust()`.
Already-fitted DESeq2-shaped list or S4-like objects are also recognized when
they expose row metadata in `mcols` with coefficient columns and matching `SE_`
columns, plus optional nested `assays$counts`, `sampleLevels`, and
`modelMatrix` diagnostics. When `factorLevels` or factor-valued `colData` is
available, character triplet contrasts use the first factor level as the
reference unless `reference` is supplied explicitly.
`resolveResultsContrastRust()` resolves character triplet, list/listValues, and
numeric contrast forms against `resultsNames`, preserving whether character or
numeric all-zero handling applies. This remains an already-fitted route; it does
not run GLM fitting from DESeqDataSet objects. Default and `name = ...`
coefficient results do not apply contrast all-zero handling.

`resultsNames()` returns coefficient/result names from primitive Wald/LRT fit
objects, list-like primitive objects, DESeq2-shaped fitted list/S4 metadata,
beta matrices, or already-supplied character vectors. It is the R-wrapper helper
for feeding `resolveResultsContrastRust()` without manually duplicating
coefficient names;
`resultsNamesRust()` provides the same primitive helper under the Rust-suffixed
name.

`nbinomWaldTestRust()` currently accepts the same already-computed primitive
Wald inputs, a list-like primitive object containing those fields, or an
already-fitted DESeq2-shaped list/S4-like object with compatible `mcols` fields,
and returns an `rsdeseq2PrimitiveWaldFit` object. For unfitted
`DESeqDataSet`-style inputs, use the Rust crate or CLI workflow entry points and
then pass the fitted diagnostics through this wrapper surface.

`lrtFitRust()` packages already-computed primitive full-model beta matrices,
coefficient covariance arrays, optional coefficient SEs, and LRT statistic /
p-value diagnostics. If `lrtPvalue` is omitted, callers can supply `lrtDf` and
the wrapper computes chi-square tail p-values. `results()` consumes
`rsdeseq2PrimitiveLrtFit` objects with the same `name` and `contrast = ...`
selection rules as primitive Wald fits: the selected coefficient or contrast
defines the reported LFC and SE, while `stat` and `pvalue` remain the supplied
LRT model-comparison statistic and p-value. `nbinomLRTRust()` accepts those
primitive fields, a list-like primitive object containing them, or an
already-fitted DESeq2-shaped list/S4-like object with compatible `mcols` fields,
and returns an `rsdeseq2PrimitiveLrtFit` object. For unfitted
`DESeqDataSet`-style inputs, use the Rust crate or CLI workflow entry points and
then pass the fitted diagnostics through this wrapper surface. `results()` can
also consume list-like primitive or fitted-object
LRT fields directly and internally package them through `nbinomLRTRust()`.

`applyIndependentFilteringRust()` applies baseMean-driven filtered BH
adjustment to a primitive result data frame. It uses the same default theta-grid
shape as the Rust core, candidate cutoffs from R type-7 quantiles, BH adjustment
after filtering, and `stats::lowess(..., f = 1/5)` metadata for threshold
selection. The returned data frame gets a `filtered` column and an
`independentFiltering` attribute with theta, rejection counts, selected
threshold, lowess fit, and alpha. This helper does not inspect DESeqDataSet
objects.

`deseq2McolsDiagnosticsRust()` is a pure-R shape contract for the Rust fit-state
diagnostic alias view. It validates primitive vectors and returns a data frame
with DESeq2-style row metadata names such as `betaConv`, `fullBetaConv`,
`reducedBetaConv`, `betaIter`, `reducedBetaIter`, `deviance`, and `maxCooks`.
`rsdeseq2DiagnosticSchemaRust()` returns the same schema and can call the
registered native C stub when the shared library is loaded. These helpers do
not inspect or modify DESeqDataSet objects.

Unsupported high-level entry points remain explicit:

- `DESeq()`
- `estimateDispersionsRust()`

High-level `DESeqDataSet` fitting/mutation integration is the remaining wrapper
boundary; the implemented wrapper surface covers primitive and already-fitted
result objects.
