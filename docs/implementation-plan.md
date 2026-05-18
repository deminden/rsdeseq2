# Implementation Plan

This plan keeps `rsdeseq2` close to DESeq2 while preserving a clean Rust
implementation. The local ignored clone at `external/DESeq2` is the inspection
reference; no DESeq2 source is vendored or translated line by line.

## Active TODOs

- [x] Create Rust workspace, core crate, R scaffold, Python placeholder, docs,
  scripts, and CI skeleton.
- [x] Implement row-major count and numeric matrix types.
- [x] Implement DESeq2-shaped `ratio` and `poscounts` size-factor estimation.
- [x] Implement normalized counts, `baseMean`, and BH adjusted p-values.
- [x] Clone DESeq2 into ignored `external/DESeq2` for local inspection.
- [x] Expose early DESeq2 row metadata: `baseMean`, `baseVar`, and `allZero`.
- [x] Add DESeq2-style weighted `baseMean`/`baseVar` helpers that multiply
  normalized counts by observation weights before ordinary row summaries.
- [x] Add builder-owned observation weights, design-aware `weights_fail`
  fit-state metadata, and supplied-dispersion Wald/LRT weighted GLM wiring.
- [x] Add builder-owned `geoMeans` and `controlGenes` options.
- [x] Add builder-owned caller-supplied size factors for fixed-size-factor
  parity tests and external R wrapper integration.
- [x] Add hand tests for DESeq2 size-factor error cases.
- [x] Add skip-safe DESeq2 golden-reference tests for generated normalization,
  supplied-dispersion Wald/LRT, and Cook's diagnostic files.
- [x] Add DESeq2-style gene/sample normalization factors for normalized counts
  and supplied-dispersion fixed Wald/LRT offsets.
- [x] Extend the current native linear-mu dispersion/Wald subset to use
  DESeq2-style normalization-factor moments starts and fitted raw means.
- [x] Add skip-safe DESeq2 golden-reference tests for normalization-factor
  native rough/moments starts and `linearModelMuNormalized` fitted means.
- [x] Extend `scripts/generate_deseq2_references.R` to write a deterministic
  tiny fixture, model matrices, version metadata, and fixed-dispersion GLM
  references from DESeq2 internals.
- [x] Extend generated references and skip-safe Rust tests for weighted
  base metadata, normalized observation weights, and supplied-dispersion
  weighted Wald/LRT outputs.
- [x] Extend generated references and skip-safe Rust tests for the current
  native weighted GLM-mu mean-trend LRT branch.
- [x] Add R package primitive matrix helpers for size factors, normalized
  counts, gene/sample normalization factors, baseMean, and weighted base
  metadata using an explicit R fallback while the native bridge is pending.
- [x] Add an opt-in registered `.Call` bridge for primitive size-factor
  estimation with fallback to the R implementation when the shared library is
  unavailable.
- [x] Add an opt-in registered `.Call` bridge for primitive base metadata with
  fallback to the R implementation when the shared library is unavailable.
- [x] Add an opt-in registered `.Call` bridge for primitive normalized counts
  with the same fallback behavior.
- [x] Add an opt-in registered `.Call` bridge for primitive baseMean with the
  same fallback behavior.
- [x] Add R package primitive result helper for Cook's cutoff masking and the
  explicit two-group low-count heuristic.
- [x] Add an opt-in registered `.Call` bridge for primitive Cook's cutoff
  masking with fallback to the R implementation when the shared library is
  unavailable.
- [x] Add R package primitive result-table assembly helper with DESeq2-shaped
  columns and BH-adjusted p-value calculation.
- [x] Add R package primitive independent-filtering helper for baseMean-driven
  filtered BH adjustment and threshold metadata.
- [x] Port selected DESeq2 size-factor validation cases into the R wrapper
  tests, citing the original test file as the behavioral source.
- [x] Add standard R package testthat entrypoint and CI coverage for source-tree
  wrapper tests plus `R CMD check`.
- [ ] Generate and commit small DESeq2 reference outputs when the R environment
  has Bioconductor DESeq2 installed.

## DESeq2 Reference Anchors

- Size factors: `external/DESeq2/R/core.R`, `estimateSizeFactorsForMatrix`.
- Size-factor tests: `external/DESeq2/tests/testthat/test_size_factor.R`.
- Early row metadata: `external/DESeq2/R/core.R`,
  `getBaseMeansAndVariances`.
- Normalization factors: `external/DESeq2/R/methods.R`,
  `counts.DESeqDataSet` and `normalizationFactors`; `external/DESeq2/R/core.R`,
  `getSizeOrNormFactors`.
- Result table shape: `external/DESeq2/R/results.R`.
- All-zero row handling in downstream stages: `external/DESeq2/R/core.R`,
  `buildMatrixWithNARows` and calls using `mcols(object)$allZero`.
- Parametric dispersion trend: `external/DESeq2/R/core.R`,
  `parametricDispersionFit` and `estimateDispersionsFit`.
- Dispersion objective derivatives and optimizer diagnostics:
  `external/DESeq2/src/DESeq2.cpp`, `log_posterior`,
  `dlog_posterior`, `d2log_posterior`, and `fitDisp`.

## Phase 1: Implemented Foundation

- [x] Count matrix validation.
- [x] Row-major numeric matrix storage.
- [x] Design matrix wrapper for R-generated model matrices.
- [x] Deterministic design matrix rank helpers and DESeq2-style full-rank
  guards for GLM-facing builder paths.
- [x] Execution modes and statistical option enums.
- [x] Error type with explicit unsupported-feature handling.
- [x] Size factors: `ratio`, `poscounts`, supplied geometric means, control
  gene index subset, and caller-supplied size factors.
- [x] Normalized counts.
- [x] Gene/sample normalization factors for normalized counts and base row
  metadata.
- [x] `baseMean`.
- [x] BH adjustment with missing-value support.
- [x] Inspectable `DeseqFit` skeleton.

## Phase 1.1: Closer Early Metadata Parity

- [x] Store `allZero` per gene in `DeseqFit`.
- [x] Store `baseVar` per gene in `DeseqFit`.
- [x] Make `baseVar` use sample variance, matching `matrixStats::rowVars`.
- [x] Return `NaN` for `baseVar` when only one sample is present, matching R
  variance behavior.
- [x] Add weighted base metadata helpers matching `getBaseMeansAndVariances`
  weighted-count preprocessing.
- [x] Add builder APIs for DESeq2-like `geoMeans` and `controlGenes`.
- [x] Support both control-gene indices and logical masks in Rust API.

## Phase 2: Fixed-Dispersion GLM

- [x] Negative-binomial log PMF matching
  `dnbinom(x, mu=mu, size=1/disp, log=TRUE)`.
- [x] Row and matrix negative-binomial log-likelihood helpers matching
  DESeq2 `nbinomLogLike`.
- [x] DESeq2-style `-2 * logLik` helper for fitted rows.
- [x] Intercept-only fixed-dispersion shortcut matching `fitNbinomGLMs`.
- [x] Initial fixed-dispersion IRLS beta fitting using the standard
  design-matrix branch.
- [x] DESeq2-style augmented QR weighted least-squares solver option for
  `fitBeta(useQR=TRUE)` foundations.
- [x] IRLS convergence criterion:
  `abs(dev - dev_old)/(abs(dev) + 0.1) < betaTol`.
- [x] Natural-log beta fitting with explicit log2 conversion for result fields.
- [x] Coefficient standard errors for the initial fixed-dispersion IRLS branch.
- [x] Hat diagonals for the initial fixed-dispersion IRLS branch.
- [x] QR branch foundation matching DESeq2's augmented least-squares update
  shape for fixed-dispersion IRLS.
- [x] Observation weights for general IRLS, matching low-level `fitBeta`
  working-weight and deviance weighting semantics.
- [x] DESeq2-style observation-weight preprocessing helper with row-max
  normalization, weighted design-rank checks, thresholded Cox-Reid sub-design
  checks, and `weights_fail` flags.
- [x] Wire builder observation weights into weighted base metadata and
  supplied-dispersion fixed Wald/LRT compact GLM fitting.
- [x] Per-coefficient natural-log-scale ridge support for IRLS, matching
  DESeq2's `diag(lambda)` shape after log2-to-natural-scale conversion.
- [x] Log2-scale beta covariance storage and primitive numeric Wald contrast
  helper using `c' beta` and `sqrt(c' Sigma c)`.
- [ ] Full DESeq2 beta-prior variance estimation, expanded model-matrix
  handling, and R-style contrast numerator/denominator construction.
- [x] Default Wald statistic/p-value for a selected coefficient.
- [x] t-distribution Wald p-values for `useT=TRUE`, including residual,
  scalar, and per-gene degrees of freedom.
- [x] Selected-coefficient LFC-threshold Wald alternatives for the current
  primitive matrix result path.
- [x] Supplied-dispersion Wald pipeline for primitive numeric contrasts with
  result rows, Cook's cutoff masking, and independent filtering.
- [x] Selected-coefficient result-row assembly with `baseMean`, LFC, SE, stat,
  p-value, and BH-adjusted p-value.
- [x] Supplied-dispersion fixed-dispersion Wald pipeline for one coefficient.
- [x] Supplied-dispersion fixed-dispersion LRT pipeline for full vs reduced
  design matrices.
- [x] DESeq2-style full-model deviance diagnostic (`-2 * logLike`) in
  `DeseqFit` for GLM-backed pipelines.
- [x] LRT fit-state diagnostics for reduced-model log likelihood, beta
  convergence, and beta iteration counts.
- [x] Gene/sample normalization factors as GLM offsets for supplied-dispersion
  fixed Wald/LRT pipelines.
- [x] All-zero row expansion for the supplied-dispersion Wald/LRT pipelines, using
  `NaN` in internal matrices and `None` in result rows.
- [x] DESeq2-style Cook's distance matrix for the supplied-dispersion Wald
  pipeline using robust method-of-moments dispersion.
- [x] `samplesForCooks` and `maxCooks` behavior for model-matrix cells with at
  least three replicates.
- [x] Cook's cutoff p-value masking and BH recomputation for result rows.
- [x] Base-mean independent filtering metadata for current Wald result rows.
- [x] Optional DESeq2-internal fixed-dispersion reference checks for Wald, LRT,
  and Cook's distances.
- [x] Lightweight DESeq2 `mcols(dds)`-style diagnostic alias view for
  implemented Wald/LRT fit-state fields.
- [x] Pure-R diagnostic metadata shape helper for future R wrapper bindings.
- [x] Minimal native diagnostic schema contract for the R bridge scaffold.
- [x] Add initial linear-mu gene-wise dispersion foundation:
  `linearModelMu`, `roughDispEstimate`, `momentsDispEstimate`, bounded starts,
  Cox-Reid objective scoring, Armijo line search, and grid fallback.
- [ ] Full result-table assembly with Bioconductor-style metadata.
- [ ] General all-zero row expansion helpers for future dispersion and full
  result-table outputs.
- [ ] Optim fallback for non-converged or unstable rows.
- [ ] Stage-by-stage comparison against DESeq2 internals on tiny datasets.

## Phase 3: Dispersion Estimation

- [x] Linear-mu projection helper matching DESeq2 `linearModelMu` shape.
- [x] Rough dispersion starts matching `roughDispEstimate`.
- [x] Moments dispersion starts matching `momentsDispEstimate`.
- [x] Normalization-factor moments starts using
  `mean(1 / colMeans(normalizationFactors))` on non-all-zero rows.
- [x] Linear-mu fitted raw means using normalization factors through
  `linearModelMuNormalized`/`getSizeOrNormFactors` semantics.
- [x] Bounded initial alpha values using `min(rough, moments)`.
- [x] Fixed-mean two-pass log-alpha grid search, matching the shape of
  DESeq2's `fitDispGridWrapper` fallback.
- [x] Cox-Reid log determinant adjustment for unweighted fixed-mean dispersion
  objectives.
- [x] DESeq2 alpha-dependent negative-binomial log-likelihood kernel for
  dispersion scoring.
- [x] First derivative of the unweighted Cox-Reid-adjusted profile likelihood.
- [x] DESeq2-shaped Armijo line-search optimizer for fixed-mean dispersions.
- [x] Grid fallback for non-converged line-search estimates above
  `minDisp * 10`.
- [x] Builder method for the current linear-mu gene-wise dispersion stage.
- [x] Log-dispersion prior objective/first-derivative/second-derivative
  support for future MAP dispersion fitting, including prior-aware line-search
  and grid entry points.
- [x] Parametric dispersion trend foundation:
  `asymptDisp + extraPois / mean`, DESeq2 row-selection rule, robust residual
  filtering, and Gamma identity-link IRLS.
- [x] Mean dispersion trend:
  shared `100 * minDisp` viability gate, `10 * minDisp` filtered trimmed mean,
  constant fitted trend expansion, and `FitType::Mean` builder dispatch.
- [x] Dispersion prior variance branches:
  MAD-squared log residual variance, `trigamma((m - p) / 2)` subtraction, and
  `0.25` floor for residual df greater than 3, no subtraction for saturated
  designs, and deterministic low-df histogram/KL matching for residual df 1
  through 3.
- [x] Initial MAP dispersion stage:
  DESeq2 `dispInit`, `log(dispFit)` prior means, line search, grid fallback,
  `dispMAP`, `dispOutlier`, and final dispersion outputs.
- [x] Initial GLM-mu mean/dispersion alternation for non-linear-mu
  designs using fixed-dispersion IRLS, fixed-mean dispersion optimization,
  DESeq2's `niter` count, and the `.05` log-dispersion `fitidx` rule.
- [x] Parametric/mean trend fitting, prior variance, and MAP shrinkage on top
  of the GLM-mu gene-wise dispersion branch.
- [x] Native Wald wiring on top of the GLM-mu MAP dispersion branch.
- [x] Second derivative of the Cox-Reid-adjusted profile likelihood.
- [x] Observation weights in the Cox-Reid objective, first derivative, and
  second derivative.
- [ ] Local and glmGamPoi dispersion trend types.
- [x] DESeq2-shaped prior variance branch for residual df 1 through 3.
- [x] Low-level observation-weighted MAP dispersion fitting.
- [x] Builder-level observation weights through the GLM-mu gene-wise
  dispersion, MAP, and native Wald path.
- [x] Optional R reference generation and skip-safe Rust checks for weighted
  GLM-mu dispersion/MAP/Wald intermediates.
- [ ] Complete weighted dispersion parity for DESeq2's non-linear-mu GLM mean
  fitting path, including generated-reference validation and remaining edge
  cases.
- [ ] glmGamPoi MAP dispersion path.

## Phase 4: Full Wald Pipeline

- [x] Run size factors, base statistics, supplied dispersion, beta fitting,
  Wald tests, p-values, and adjusted p-values in a single inspectable state for
  one selected coefficient.
- [x] Add native dispersion estimation to that pipeline for the limited
  linear-mu, parametric-trend, deterministic-prior, no-weight MAP branch.
- [x] Generalize the limited native MAP/Wald branch to builder-selected
  parametric or mean dispersion trends.
- [x] Generalize the limited native MAP/Wald branch to the GLM-mu mean-refit
  dispersion branch, including builder-level observation weights.
- [ ] Generalize the native Wald pipeline to DESeq2's full dispersion and GLM
  fitting behavior.
- [x] Add R wrapper around primitive matrices for implemented normalization
  and early base-metadata helpers.
- [ ] Add R wrapper around DESeqDataSet inputs without faking unsupported S4
  behavior.

## Phase 5: Results Compatibility

- [x] Primitive numeric Wald linear contrasts.
- [x] Result-row assembly for precomputed primitive numeric Wald contrasts.
- [x] Builder-level supplied-dispersion Wald pipeline for primitive numeric
  contrasts.
- [x] Primitive coefficient-name, positive/negative coefficient-list, and
  common factor-level contrast resolution against design coefficient names.
- [x] Numeric/expanded `contrastAllZero` behavior for primitive Wald contrasts,
  matching DESeq2's `contrastAllZeroNumeric` model-matrix selection rule.
- [x] Character/factor-level `contrastAllZero` behavior for primitive
  factor-level Wald contrasts with caller-supplied sample levels, matching
  DESeq2's `contrastAllZeroCharacter` selected-sample rule.
- [ ] Full DESeq2 contrast handling for `results(contrast=...)`, including
  colData/formula-aware factor-level semantics, complete coefficient-name
  cleanup, and contrast-aware Cook's/refit edge cases.
- [x] Initial Cook's distance and `maxCooks` diagnostics.
- [x] Cook's outlier p-value filtering.
- [x] Explicit primitive helper for DESeq2's two-group low-count Cook's
  heuristic, to be called only when R-side formula/colData semantics establish
  the one-factor two-level condition.
- [x] Primitive Cook's outlier count replacement transform: trimmed normalized
  means, size-factor/normalization-factor rescaling, integer truncation,
  replaceable-sample mask, and `replace` flags.
- [x] Replacement-refit planning metadata: replacement-count base metadata,
  `nrefit`, `refitReplace`, `newAllZero`, and post-refit `maxCooks` masking.
- [x] Limited Cook's replacement-refit path for GLM-mu native Wald and LRT:
  original fit, replacement counts, replacement-count refit with original size
  factors, `refitReplace` merge, `newAllZero` result clearing, and final
  filtering.
- [ ] Full Cook's outlier replacement behavior with DESeq2-style replacement
  refit for contrasts, beta priors, and Bioconductor object metadata.
- [ ] Full formula-aware outlier handling and R wrapper integration.
- [x] Initial independent filtering.
- [x] R `stats::lowess`-shaped independent-filter threshold selection for the
  DESeq2 default theta grid, including fitted-curve metadata.
- [x] R `lowess` parity for dense custom theta grids where R's
  `delta` interpolation shortcut skips fitted points.
- [x] Initial fixed-dispersion LRT.
- [x] Limited native-dispersion LRT for current linear-mu and GLM-mu MAP
  dispersion branches.
- [x] Optional DESeq2-internal native weighted GLM-mu LRT reference hook for
  the current mean-trend branch.
- [ ] Full LRT parity with native dispersion reference outputs, local/glmGamPoi
  trends, optim fallback, and remaining edge cases.

## Phase 6: Secondary Features

- [ ] Local dispersion trend.
- [x] Mean fit type.
- [ ] glmGamPoi-like mode if feasible.
- [ ] VST.
- [ ] rlog.
- [ ] Python package.
- [ ] Mature CLI.

## Engineering Rules

- Keep unsupported stages explicit with `DeseqError::UnsupportedFeature`.
- Prefer strict, deterministic behavior before adding fast paths.
- Preserve row-major gene-contiguous storage in Rust.
- Add hand-computable tests before golden-reference tests.
- Document every intentional deviation from DESeq2 semantics.
