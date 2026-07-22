# rsdeseq2 R Package

This package provides R access to selected `rsdeseq2` primitives and native
bridges. It does not fit complete `DESeqDataSet` workflows or replace the full
DESeq2/Bioconductor interface.

The package is tested with R 4.6.1. Each numerical reference fixture records
the R version that generated it.

## Install and Try a Primitive

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
estimateSizeFactorsRust(counts, native = TRUE)
```

The Rust core is validated against saved R 4.6.1 / DESeq2 1.52.0 outputs. The
repository [validation summary](../../README.md#measured-validation) records
the absolute measurements and their scope; the
[benchmark documentation](../../docs/benchmarks.md) contains unrounded data,
baselines, and reproduction commands. Those workflow-level measurements
describe the native Rust core, not additional surface implemented by this R
package.

## Supported Primitive Helpers

- `estimateSizeFactorsRust()`
- `normalizedCountsRust()`
- `baseMeanRust()`
- `baseMetadataRust()`
- `applyCooksCutoffRust()`
- `resultsTableRust()`
- `waldFitRust()`
- `lrtFitRust()`
- `nbinomWaldTestRust()`
- `nbinomLRTRust()`
- `resultsNames()`
- `resultsNamesRust()`
- `resolveResultsContrastRust()`
- `applyIndependentFilteringRust()`
- `rsdeseq2DiagnosticSchemaRust()`
- `deseq2McolsDiagnosticsRust()`

These helpers include an R-level implementation that mirrors the
Rust-supported normalization algorithms, including gene/sample normalization
factors, weighted base metadata, Cook's cutoff result masking, and diagnostic
metadata shape. The R API does not expose the full Rust workflow. Registered
primitive `.Call` bridges cover selected normalization, metadata, Cook's, and
diagnostic operations.
`baseMetadataRust(native = TRUE)` can use the registered package `.Call` bridge
for that primitive and uses the R-level implementation when the shared library
is unavailable. `estimateSizeFactorsRust(native = TRUE)`,
`normalizedCountsRust(native = TRUE)`, and `baseMeanRust(native = TRUE)` follow
the same pattern for primitive normalization summaries.
`applyCooksCutoffRust(native = TRUE)` can use the registered native bridge for
Cook's masking while R retains BH adjustment and output assembly.
`resultsTableRust()` assembles already-computed primitive vectors into a
DESeq2-shaped result data frame and computes BH-adjusted p-values when needed.
`waldFitRust()` packages primitive Wald beta/covariance outputs.
`nbinomWaldTestRust()` accepts those same already-computed primitive
Wald inputs, or a list-like primitive object, and returns the same fit object.
`lrtFitRust()` and `nbinomLRTRust()` do the corresponding primitive packaging
for already-computed LRT statistics: `results()` reports the selected
coefficient or contrast LFC/SE while preserving the LRT statistic and p-value.
`results()` can report selected coefficients or DESeq2-style
`contrast = ...` requests from primitive Wald or LRT fit objects, or from
list-like objects containing the same primitive fields. It also recognizes
already-fitted DESeq2-shaped list/S4 metadata with coefficient columns and
matching `SE_` columns in `mcols`. These routes include
`lfcThreshold`/`altHypothesis` Wald p-values with Normal or supplied t tails.
When fitted `factorLevels` or factor-valued `colData` metadata is available,
character triplet contrasts use the first factor level as the reference unless
`reference` is supplied explicitly.
Character, list, and numeric contrast forms are resolved by
`resolveResultsContrastRust()` with matching all-zero handling metadata.
Named coefficient results do not apply contrast all-zero handling.
`resultsNames()` exposes coefficient names from primitive Wald/LRT fit objects,
list-like primitive objects, fitted-object metadata, beta matrices, or validated
character vectors for the same contrast-resolution workflow.
`resultsNamesRust()` is the same primitive helper under the Rust-suffixed name.
`applyIndependentFilteringRust()` applies baseMean-driven filtered BH
adjustment to primitive result tables and records filtering metadata as an
attribute.
DESeqDataSet fitting integration is unsupported and returns clear errors.
