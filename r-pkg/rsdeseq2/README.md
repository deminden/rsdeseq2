# rsdeseq2 R Package Scaffold

This package is the R-first access layer for `rsdeseq2`.

Currently implemented primitive matrix helpers:

- `estimateSizeFactorsRust()`
- `normalizedCountsRust()`
- `baseMeanRust()`
- `baseMetadataRust()`
- `applyCooksCutoffRust()`
- `resultsTableRust()`
- `applyIndependentFilteringRust()`
- `rsdeseq2DiagnosticSchemaRust()`
- `deseq2McolsDiagnosticsRust()`

These helpers use an R fallback implementation that mirrors the Rust-supported
normalization algorithms, including gene/sample normalization factors,
weighted base metadata, Cook's cutoff result masking, and diagnostic metadata
shape while the native Rust bridge is still being wired.
`baseMetadataRust(native = TRUE)` can use the registered package `.Call` bridge
for that primitive and falls back to the R implementation when the shared
library is unavailable. `estimateSizeFactorsRust(native = TRUE)`,
`normalizedCountsRust(native = TRUE)`, and `baseMeanRust(native = TRUE)` follow
the same pattern for primitive normalization summaries.
`applyCooksCutoffRust(native = TRUE)` can use the registered native bridge for
Cook's masking while R retains BH adjustment and output assembly.
`resultsTableRust()` assembles already-computed primitive vectors into a
DESeq2-shaped result data frame and computes BH-adjusted p-values when needed.
`applyIndependentFilteringRust()` applies baseMean-driven filtered BH
adjustment to primitive result tables and records filtering metadata as an
attribute.
DESeqDataSet integration and full DESeq2-style result generation remain
unsupported and return clear errors.
