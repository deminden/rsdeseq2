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
- `fixed_mu_reduced.tsv`
- `fixed_hat_reduced.tsv`
- `fixed_cooks_full.tsv`
- `fixed_weighted_mu_reduced.tsv`
- `fixed_weighted_hat_reduced.tsv`
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
They validate the implemented Rust fixed-dispersion GLM path. They are
not a substitute for full DESeq2 dispersion-estimation parity. The default
reference set includes the unweighted fixed Wald/LRT, fitted `mu`, hat
diagonal, Cook's distance files, and weighted fixed Wald/LRT files including
reduced-model fitted means and hat diagonals because these are numerically
reproduced. The same supplied-dispersion MLE beta matrix
also drives `beta_prior_variance_reference.tsv`, which records DESeq2
`estimateBetaPriorVar` weighted and quantile beta-prior variance outputs for
the primitive Rust beta-prior estimator, plus `beta_prior_refit_reference.tsv`,
`beta_prior_refit_mu.tsv`, and `beta_prior_refit_hat.tsv` for the corresponding
supplied-dispersion ridge refit. The forced-optim fallback reference includes
both fitted means and hat diagonals. The GLM-mu native files include the default
unweighted and weighted Cox-Reid gene-wise branches,
`estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE)`, the default
unweighted Cox-Reid mean-trend MAP/Wald/LRT branch, plus the supported
mean-trend MAP/Wald/LRT path:
`estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE)`, mean trend
or local trend fitting, `estimateDispersionsMAP(useCR=FALSE)`, and final
full/reduced `fitNbinomGLMs(useQR=FALSE,useOptim=FALSE)`. The matched GLM-mu Wald/LRT
reference rows include compact DESeq2-shaped result tables with BH-adjusted
p-values for result-table parity. The unweighted GLM-mu local-trend fixture
also records MAP, Wald, LRT, and result rows for the tiny-data case where one
row is usable for the local fit; the weighted GLM-mu local fixture covers the
same MAP/Wald/LRT results with `weightsFail` expansion. A separate
GLM-mu Cox-Reid local-trend fixture family records unweighted and weighted MAP
dispersion intermediates, stored dispersion means, Wald fitted means and hat
diagonals, LRT full/reduced likelihoods, and compact Wald/LRT result rows for
that combination.

Rust golden tests skip automatically when these files are absent. After running
the R script, those same tests compare size factors, normalized counts,
normalization-factor normalized counts, baseMean/baseVar/allZero,
fixed-dispersion Wald/LRT fields, fitted means, hat diagonals, Cook's
distances, matched GLM-mu result-row p-values and adjusted p-values, and
compact matched GLM-mu Wald/LRT result rows, and Cook's replacement/refit
bookkeeping against the generated references.

The separate L-BFGS-B stress reference is generated with the recorded R 4.6.1
environment:

```bash
OPENBLAS_NUM_THREADS=1 \
  Rscript scripts/generate_lbfgsb_synthetic_stress_fixtures.R
```

The generator refuses other R versions and serializes doubles with 17
significant digits. The 2026-07-22 fixture contains 512 cases and records x86_64
Linux, OpenBLAS 0.3.32, one thread. Replaying it with `rcompat-lbfgsb` 0.2.1
produces 512/512 exact endpoints, objective values, and evaluation counts. The
0.1.6 baseline produced 493/512 practical endpoint-plus-objective matches,
507/512 practical objective matches, 0/512 exact endpoint-plus-objective
matches, and 311/512 exact count matches.

## Versioned High-Error Benchmark

Optimizer isolation does not prove end-to-end parity, so the repository also
stores the 100 rows with the largest v0.2.4 errors from a versioned 69,045-gene
validation contrast. Each row records the result column with the largest
baseline error. The fixture is
[`data/wald_frozen_worst100_r461.tsv`](data/wald_frozen_worst100_r461.tsv).

Build in release mode, generate all 69,045 diagnostics from the saved real-data
reference, and run the accuracy check:

```bash
cargo build --release -p rsdeseq2
python3 scripts/real_data_parity.py \
  --study-root /path/to/study-inputs-and-reference-outputs \
  --binary target/release/rsdeseq2 \
  --contrast VALIDATION_CONTRAST_ID \
  --contrast-size-factors estimate \
  --output results/benchmarks/frozen_worst100.tsv \
  --diagnostics-output results/benchmarks/frozen_worst100_diagnostics.tsv \
  --diagnostics-limit 69045
python3 scripts/score_frozen_worst_genes.py \
  --fixture docs/data/wald_frozen_worst100_r461.tsv \
  --diagnostics results/benchmarks/frozen_worst100_diagnostics.tsv \
  --report-only
```

`VALIDATION_CONTRAST_ID` is the contrast identifier used by the selected study
input and saved-reference bundle.

The full diagnostics table lets the scorer evaluate every fixed gene from its
measured row. The v0.2.5 release-mode run measured median, mean, and maximum
absolute errors of `6.063612945084174e-10`, `1.260818774570247e-4`, and
`1.5259081158007781e-3`. The frozen v0.2.4 baseline measured
`1.4637657972313423e-4`, `3.793917566690452e-4`, and
`3.0938714191082184e-3`, respectively. The v0.2.5 median and mean are
241401.589x and 3.00909032x lower; the maximum is 50.6796531% lower, with a
v0.2.5-to-v0.2.4 ratio of `0.493203469`. Of the 100 rows, 89 improved and 78
improved by at least 10x.

The ratio size-factor path uses compensated accumulation of log counts when
computing per-gene geometric means. This reduces floating-point loss in
normalization before the dispersion and fitting stages.

The native Wald pipeline is covered by Rust self-consistency tests:
it preserves dispersion intermediates, expands all-zero rows, and produces the
same GLM/result fields as the supplied-dispersion Wald path when fed its own
final MAP dispersions. The reference generator emits normalization-factor
native dispersion references for `roughDispEstimate`, `momentsDispEstimate`,
bounded starts, and post-`minmu` fitted means from
`estimateDispersionsGeneEst`, plus weighted GLM-mu references for the supported
mean-trend MAP/Wald/LRT branch. Generated reference coverage excludes broader
DESeq2 stage-by-stage comparisons for unimplemented dispersion branches.
Thresholded selected-coefficient Wald alternatives are covered by
hand/R-formula tests; generated `results(lfcThreshold=...)` tables are absent.

Rust dispersion unit tests use hand-computable cases. A complete generated
stage-reference set would include separate columns for `roughDispEstimate`, `momentsDispEstimate`,
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
