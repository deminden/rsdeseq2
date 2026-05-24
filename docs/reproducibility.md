# Reproducibility

DESeq2 compatibility requires more than matching final p-values. Intermediate
outputs should be generated and compared stage by stage.

## Reference Generation

Run:

```bash
Rscript scripts/generate_deseq2_references.R
```

The script writes references under:

```text
crates/rsdeseq2/tests/data/deseq2_reference/
```

Generated files should include:

- `metadata.tsv`
- `sessionInfo.txt`
- `counts.tsv`
- `col_data.tsv`
- `design_full.tsv`
- `design_reduced.tsv`
- `fixed_dispersions.tsv`
- `size_factors_ratio.tsv`
- `size_factors_poscounts.tsv`
- `normalized_counts_ratio.tsv`
- `base_mean_ratio.tsv`
- `base_metadata_ratio.tsv`
- `normalization_factors.tsv`
- `normalized_counts_nf.tsv`
- `base_metadata_nf.tsv`
- `native_nf_dispersion_reference.tsv`, when DESeq2 internals are available
- `native_nf_mu.tsv`, when DESeq2 internals are available
- `parametric_trend_reference.tsv`
- `dispersion_prior_variance_reference.tsv`
- `map_dispersion_reference.tsv`
- `results_wald_ratio.tsv`, when full DESeq2 succeeds on the fixture
- `fixed_wald_reference.tsv`, when DESeq2 internals are available
- `fixed_wald_t_reference.tsv`, when DESeq2 internals are available
- `fixed_lrt_reference.tsv`, when DESeq2 internals are available
- `fixed_mu_full.tsv`
- `fixed_hat_full.tsv`
- `fixed_cooks_full.tsv`
- `native_glm_mu_cr_reference.tsv`, when DESeq2 internals are available
- `native_glm_mu_cr_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_reference.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_lrt_reference.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_wald_mu.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_wald_hat.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_cr_map_reference.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_cr_lrt_reference.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_cr_map_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_cr_wald_mu.tsv`, when DESeq2 internals are available
- `native_glm_mu_mean_cr_wald_hat.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_reference.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_lrt_reference.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_cr_reference.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_cr_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_mean_cr_map_reference.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_mean_cr_lrt_reference.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_mean_cr_map_dispersion_mu.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_mean_cr_wald_mu.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_mean_cr_wald_hat.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_wald_mu.tsv`, when DESeq2 internals are available
- `native_weighted_glm_mu_wald_hat.tsv`, when DESeq2 internals are available
- `cooks_replacement_counts.tsv`
- `cooks_replacement_design.tsv`
- `cooks_replacement_cooks.tsv`
- `cooks_replacement_size_factors.tsv`
- `cooks_replacement_options.tsv`
- `cooks_replacement_candidate_counts.tsv`
- `cooks_replacement_replaced_counts.tsv`
- `cooks_replacement_rows.tsv`

The fixed-dispersion files use `DESeq2:::fitNbinomGLMs` with supplied
dispersions, default `1e-6` beta ridge, `useQR=FALSE`, and `useOptim=FALSE`.
They exist to validate the current Rust fixed-dispersion GLM path. They are
not a substitute for full DESeq2 dispersion-estimation parity. The default
reference set includes the unweighted fixed Wald/LRT, fitted `mu`, hat
diagonal, Cook's distance files, and weighted fixed Wald/LRT files because
these are numerically reproduced. The GLM-mu native files include the default
unweighted and weighted Cox-Reid gene-wise branches,
`estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE)`, the default
unweighted Cox-Reid mean-trend MAP/Wald/LRT branch, plus the current narrow
mean-trend MAP/Wald/LRT path:
`estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE)`, mean trend
or local trend fitting, `estimateDispersionsMAP(useCR=FALSE)`, and final
full/reduced `fitNbinomGLMs(useQR=FALSE,useOptim=FALSE)`. The matched GLM-mu Wald/LRT
reference rows include compact DESeq2-shaped result tables with BH-adjusted
p-values for result-table parity. The unweighted GLM-mu local-trend fixture
also records MAP, Wald, LRT, and result rows for the tiny-data case where one
row is usable for the local fit; the weighted GLM-mu local fixture covers the
same MAP/Wald/LRT/result-row surface with `weightsFail` expansion. A separate
unweighted GLM-mu Cox-Reid local-trend fixture records the MAP dispersion
intermediate and stored dispersion means for that combination.

Rust golden tests skip automatically when these files are absent. After running
the R script, those same tests compare size factors, normalized counts,
normalization-factor normalized counts, baseMean/baseVar/allZero,
fixed-dispersion Wald/LRT fields, fitted means, hat diagonals, Cook's
distances, matched GLM-mu result-row p-values and adjusted p-values, and
compact matched GLM-mu Wald/LRT result rows, and Cook's replacement/refit
bookkeeping against the generated references.

The current native Wald pipeline is covered by Rust self-consistency tests:
it preserves dispersion intermediates, expands all-zero rows, and produces the
same GLM/result fields as the supplied-dispersion Wald path when fed its own
final MAP dispersions. The reference generator now emits normalization-factor
native dispersion anchors for `roughDispEstimate`, `momentsDispEstimate`,
bounded starts, and post-`minmu` fitted means from
`estimateDispersionsGeneEst`, plus weighted GLM-mu anchors for the current
mean-trend MAP/Wald/LRT branch. Future R references should extend this path to
broader DESeq2 stage-by-stage comparisons once the remaining dispersion
branches are implemented.
Thresholded selected-coefficient Wald alternatives are currently covered by
hand/R-formula tests; future references should add `results(lfcThreshold=...)`
tables for the supported alternatives.

The current Rust dispersion tests are hand-computable. Future reference files
should add separate columns for `roughDispEstimate`, `momentsDispEstimate`,
linear fitted means, Cox-Reid objective values, log-dispersion prior objective
values, GLM-refit mean matrices from the non-`linearMu` branch, `fitidx`
decisions across `niter`, parametric trend coefficients, mean trend constants,
`dispFit`, fallback grid estimates, line-search estimates, weighted objective
values, and MAP dispersions. Line-search references should include final
score, first derivative, and second derivative so drift can be localized before
full result-table comparison.

## Version Recording

Every reference output set should record:

- DESeq2 version
- Bioconductor version
- R version
- platform
- package versions that affect fitting

## Why Intermediates Matter

Final p-values can agree accidentally or disagree for reasons that are hard to
diagnose. Comparing `sizeFactors`, normalized counts, dispersion estimates,
coefficients, standard errors, Cook's distances, and filtering decisions gives a
clear path to localizing drift.
