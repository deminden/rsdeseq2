# R Wrapper

The R package exposes primitive matrix helpers for the supported early
normalization stages:

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
and early base-metadata algorithms. The R API does not expose full Rust-core
workflow execution. Registered primitive `.Call` bridges are
available for the selected operations described below.
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
step can use the registered `.Call` bridge; R assembles the output table
and computes BH-adjusted p-values. When p-values carry gene names, named
`maxCooks`, count rows, and Cook's rows are aligned to that order before
masking; named Cook's columns are aligned to named count columns before the
low-count heuristic is evaluated.

`resultsTableRust()` assembles already-computed primitive vectors into a
DESeq2-shaped result data frame with `baseMean`, `log2FoldChange`, `lfcSE`,
`stat`, `pvalue`, and `padj`, plus optional `dispersion` and `converged`
diagnostic columns. It computes BH-adjusted p-values from `pvalue` when `padj`
is omitted. When output row names are known from `rowNames` or named
`baseMean`, named result vectors are aligned to that row order before the table
is assembled. It does not run Wald/LRT fitting, Cook's masking, outlier
replacement, or independent filtering.

Rust result-table metadata retains the resolved numeric contrast vector for
implemented contrast and replacement/refit routes. Wrapper code can use this
structured provenance instead of recovering the contrast from display labels;
primitive `resultsTableRust()` callers can also pass finite numeric
coefficient-shaped `contrast` metadata directly.

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
columns, plus optional nested `assays$counts` or `assays[["counts"]]`,
`sampleLevels`, and `modelMatrix` diagnostics. Extracted row metadata can be a
base `data.frame` or a DataFrame-like object that can be safely coerced with
`as.data.frame()`, which covers common `mcols`/`rowData`-style containers.
When explicit `factorReferences` metadata is available, character triplet
contrasts use it as the treatment/base reference unless `reference` is supplied
explicitly. Otherwise, when `factorLevels` or factor-valued `colData` is
available, character triplet contrasts use the first factor level. When
counts are available and `sampleLevels` is not supplied, character contrast
all-zero handling can use the matching factor column from `colData` as the
per-sample level vector; if count columns and `colData` rows are both named,
`colData` is aligned to the count column order before that per-sample metadata
is used. Extracted `colData` can be a base `data.frame` or a DataFrame-like
object that can be safely coerced with `as.data.frame()`. Named `sampleLevels`
vectors are aligned to named count columns with the same unique-name and
missing-name validation. `factorLevels` and `colData` factor lookup is exact
first, then an
unambiguous R-cleaned alias, so non-syntactic metadata names can be addressed
like cleaned coefficient names. `factorLevels` names must themselves have
unique R-cleaned aliases; explicit metadata and levels inferred from `colData`
are rejected at fit packaging time when two factor fields would clean to the
same name. Character contrast levels used for all-zero
handling follow the same exact-first, unambiguous R-cleaned alias rule against
observed sample levels, so cleaned level requests can select raw
non-syntactic factor values without guessing through collisions. Explicit
`reference` values for character contrasts are canonicalized through the same
level-alias rule before coefficient contrast resolution. `mcols` beta
and `SE_` coefficient columns use the same exact-first, unambiguous R-cleaned
alias lookup when explicit
`resultsNames` are available, which lets ordinary cleaned data-frame column
names back non-syntactic result names without silently picking ambiguous
columns. Object-style design matrices can be exposed as fields, slots, or
attributes named `modelMatrix`, `model.matrix`, `designMatrix`, or
`fullModelMatrix`; after extraction they are aligned to named count samples
when possible. Named model-matrix columns must be unique and are also aligned
to result-name order with the same exact-first, unambiguous R-cleaned alias
lookup before numeric contrast all-zero handling uses them. Coefficient
matrices and names can also use common object aliases such as `coef`,
`coefficients`, `coefNames`,
`coefficientNames`, or `resultNames`, while covariance arrays can use
`betaCov`, `coefCovariance`, or `coefficientCovariance`. Factor-level metadata can be exposed as
`factorLevels`, `factorLevelNames`, `factorLevelsList`, or `levelNames`.
Explicit factor-reference metadata can be exposed as `factorReferences`,
`factorReference`, `factorReferenceLevels`, or `referenceLevels`; references
are exact-first, then unambiguous R-cleaned aliases against `factorLevels` when
those levels are available, and against observed sample levels when declared
factor levels are absent. Primitive fit constructors infer factor levels from
factor-valued or character `colData` before canonicalizing supplied factor
references, so cleaned reference metadata can be stored under the canonical raw
factor and level names. The canonical raw reference level is stored in the fit
object after resolution. The
same alias names are checked on object attributes after accessors and direct
fields, which matches DESeq2-shaped fitted objects that retain supplied model
matrices outside ordinary row metadata. When explicit `resultsNames` and named
coefficient-matrix columns are supplied together, the coefficient matrix is
aligned to result-name order through exact-first, unambiguous R-cleaned alias
lookup before canonical names are restored. For primitive fit constructors, named
row-shaped inputs are aligned to named beta rows with duplicate and missing
gene-name validation. This covers `baseMean`, coefficient SE matrices,
coefficient covariance arrays, count matrices, and LRT statistic / p-value /
degrees-of-freedom vectors, so row metadata follows gene names when both sides
expose them. Named coefficient SE columns and named coefficient axes of
covariance arrays are also aligned to result-name order with the exact-first,
unambiguous R-cleaned alias resolver before canonical `resultsNames` are
restored. When scalar row diagnostics are extracted from `mcols`, explicit
`mcols` row names are preserved on extracted vectors and matrix-like payloads
so the same alignment rules apply to object-shaped metadata.
`resolveResultsContrastRust()` resolves character triplet, list/listValues, and
numeric contrast forms against `resultsNames`, preserving whether character or
numeric all-zero handling applies. This remains an already-fitted route; it does
not run GLM fitting from DESeqDataSet objects. Default and `name = ...`
coefficient results do not apply contrast all-zero handling. Coefficient names
created from supported formula numeric transforms, such as
`cell type_as.double`, use the same exact-first, unambiguous R-cleaned alias
lookup for `name = ...` and coefficient-list contrasts. Exact duplicate
`resultsNames` are rejected at packaging time, while raw names that merely
share a cleaned alias remain valid so exact lookup can win before alias lookup.

`resultsNames()` returns coefficient/result names from primitive Wald/LRT fit
objects, list-like primitive objects, DESeq2-shaped fitted list/S4 metadata,
beta matrices, or already-supplied character vectors. It is the R-wrapper helper
for feeding `resolveResultsContrastRust()` without manually duplicating
coefficient names. List-like objects use the same coefficient/name aliases as
`results()`, including `coef`, `coefficients`, `coefNames`,
`coefficientNames`, and `resultNames`, so formula-derived non-syntactic names
are preserved through helper calls. When alias-provided coefficient names and a
coefficient matrix are both present, their dimensions are validated together;
`resultsNamesRust()` provides the same primitive helper under the Rust-suffixed
name.

`nbinomWaldTestRust()` accepts the same already-computed primitive
Wald inputs, a list-like primitive object containing those fields, or an
already-fitted DESeq2-shaped list/S4-like object with compatible `mcols` fields,
and returns an `rsdeseq2PrimitiveWaldFit` object. For unfitted
`DESeqDataSet`-style inputs, use the Rust crate or CLI workflow entry points and
then pass the fitted diagnostics through the R API.

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
then pass the fitted diagnostics through the R API. `results()` can
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

`deseq2McolsDiagnosticsRust()` is a pure-R schema for the Rust fit-state
diagnostic aliases. It validates primitive vectors and returns a data frame
with DESeq2-style row metadata names such as `betaConv`, `fullBetaConv`,
`reducedBetaConv`, `betaIter`, `reducedBetaIter`, `deviance`, and `maxCooks`.
`rsdeseq2DiagnosticSchemaRust()` returns the same schema and can call the
registered native C stub when the shared library is loaded. These helpers do
not inspect or modify DESeqDataSet objects.

Unsupported high-level entry points remain explicit:

- `DESeq()`
- `estimateDispersionsRust()`

High-level `DESeqDataSet` fitting and mutation are unsupported. The R API
accepts primitive and already-fitted result objects.
