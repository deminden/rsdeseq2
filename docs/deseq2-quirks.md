# DESeq2 compatibility quirks

Deliberate behaviours that reproduce DESeq2 1.52.0 even where they look
surprising. Each row names the DESeq2 source, the Rust location and the test
that pins it. Do not "fix" a row without a NUMERIC task and new references.

| Behaviour | DESeq2 source | Rust location | Pinning test |
| --- | --- | --- | --- |
| Gene-wise dispersion with observation weights: `fitDisp` receives the full non-all-zero weight matrix while counts and means are subset by `fitidx`, and the C++ indexes weights by compact row position. Rust passes `weights.row(fitting_gene_order[compact_row])` to reproduce the same pairing. | `R/core.R` `estimateDispersionsGeneEst` → `fitDispWrapper`; `src/DESeq2.cpp` `fitDisp` | `crates/rsdeseq2/src/dispersion/gene_est/implementation/estimate.rs` (refit loop) | `tests/observation_weights.rs`, golden `wald_weights` / `lrt_weights` |
| Geometric means for ratio and poscounts size factors use a compensated sum of logs to reproduce `rowMeans(log(x))`'s extended-precision accumulation, so near-boundary genes keep their optimizer routes. | `R/core.R` `estimateSizeFactorsForMatrix` (`rowMeans(log(counts))`) | `crates/rsdeseq2/src/normalization.rs` `checked_compensated_sum` | `tests/size_factors.rs`, GTEx compact panel |
| Local VST with normalization factors uses DESeq2's mean-of-inverse-size-factor approximation rather than per-cell factors. | `R/vst.R` `getVarianceStabilizedData` (`fitType == "local"`) | `crates/rsdeseq2/src/transform/vst.rs` | `tests/transform.rs` |
| Cook's replacement refits use the replaced counts with the original size factors and do not rebuild normalization, as `refitWithoutOutliers` does. | `R/core.R` `refitWithoutOutliers` | `crates/rsdeseq2/src/results/implementation/cooks_replacement.rs` | `tests/replacement_execution.rs` |
