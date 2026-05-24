# Implementation Plan

This plan keeps `rsdeseq2` close to DESeq2 while preserving a clean Rust
implementation. The local ignored clone at `external/DESeq2` is the inspection
reference; no DESeq2 source is vendored or translated line by line.

## Active TODOs

- [x] Create Rust workspace, core crate, docs, scripts, and CI skeleton.
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
  parity tests and external caller integration.
- [x] Add hand tests for DESeq2 size-factor error cases.
- [x] Add fitted dispersion trend type labels to DESeq2-shaped fit diagnostics.
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
- [x] Promote generated DESeq2 1.46.0 references for normalization-factor
  native dispersion, weighted fixed-dispersion Wald/LRT, and weighted GLM-mu
  mean-trend MAP/Wald/LRT into the default passing fixture set.
- [x] Add DESeq2-generated BH-adjusted p-value columns and result-row padj
  parity checks for the matched GLM-mu Wald/LRT fixture branches.
- [x] Add compact DESeq2-shaped result-table fixtures for the matched GLM-mu
  Wald/LRT branches and compare public Rust result rows against them.
- [x] Add unweighted GLM-mu local-trend MAP/Wald/LRT/result-table fixtures and
  handle the single-usable-row local fit edge case as a constant local trend.
- [x] Add weighted GLM-mu local-trend MAP/Wald/LRT/result-table fixtures with
  `weightsFail` row expansion.
- [x] Add primitive result-table column schema helpers for Rust APIs.
- [x] Remove current Python wrapper scaffold from the active workspace.
- [x] Restore the experimental R package scaffold and R CI surface.
- [ ] Mature the R package wrapper after core parity improves. Mature wrapper
  paths must call the Rust implementation and must not fall back to
  R/Bioconductor DESeq2 for runtime computation.
- [x] Generate and commit small DESeq2 reference outputs when the R environment
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
- [x] Design matrix wrapper for caller-supplied model matrices.
- [x] Deterministic design matrix rank helpers and DESeq2-style full-rank
  guards for GLM-facing builder paths.
- [x] Execution modes and statistical option enums.
- [x] Error type with explicit unsupported-feature handling.
- [x] Size factors: `ratio`, `poscounts`, supplied geometric means, control
  gene index subset, and caller-supplied size factors.
- [x] Normalized counts.
- [x] TSV export for raw and normalized count matrices.
- [x] Gene/sample normalization factors for normalized counts and base row
  metadata.
- [x] TSV export for gene/sample normalization-factor matrices.
- [x] `baseMean`.
- [x] BH adjustment with missing-value support.
- [x] Inspectable `DeseqFit` skeleton.

## Phase 1.1: Closer Early Metadata Parity

- [x] Store `allZero` per gene in `DeseqFit`.
- [x] Store `baseVar` per gene in `DeseqFit`.
- [x] Make `baseVar` use sample variance, matching `matrixStats::rowVars`.
- [x] Return `NaN` for `baseVar` when only one sample is present, matching R
  variance behavior.
- [x] TSV export for implemented early row metadata: `baseMean`, `baseVar`,
  and `allZero`.
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
- [x] Public math distribution aliases and helper namespace for the
  implemented DESeq2-parameterized negative-binomial primitives.
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
- [x] Public `fit_irls` dispatcher that uses the DESeq2-style intercept-only
  shortcut when eligible and general fixed-dispersion IRLS otherwise.
- [x] Public fixed-dispersion `estimate_beta` wrapper over the implemented beta
  fitting dispatcher.
- [x] Public supplied-dispersion GLM wrapper over the implemented
  fixed-dispersion IRLS path.
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
- [x] DESeq2-style beta-prior variance estimation for primitive MLE beta
  matrices, including quantile and weighted quantile methods, finite-beta
  filtering, and wide intercept priors.
- [x] DESeq2 `estimateBetaPriorVar` fixture checks for primitive beta-prior
  variance estimation on supplied-dispersion MLE beta matrices.
- [x] DESeq2 beta-prior ridge refit fixture checks for supplied-dispersion GLM
  betas, SEs, log-likelihoods, fitted means, and hat diagonals.
- [x] Combined estimated-prior beta refit fixture check that runs MLE fitting,
  estimates beta-prior variance, and compares the refit against DESeq2 anchors.
- [x] Primitive beta-prior GLM refit using supplied or estimated log2-scale
  beta-prior variances, with DESeq2's `1 / betaPriorVar / log(2)^2`
  natural-log ridge conversion, including size-factor, normalization-factor,
  and observation-weight fixed-dispersion GLM inputs.
- [ ] Expanded model-matrix beta-prior averaging, high-level beta-prior
  workflow plumbing, and DESeq2-style contrast numerator/denominator
  construction.
- [x] Default Wald statistic/p-value for a selected coefficient.
- [x] t-distribution Wald p-values for `useT=TRUE`, including residual,
  scalar, and per-gene degrees of freedom.
- [x] Selected-coefficient LFC-threshold Wald alternatives for the current
  primitive matrix result path.
- [x] Public selected-coefficient Wald wrapper over the implemented
  coefficient-with-options path.
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
  convergence, beta iteration counts, fitted means, and hat diagonals.
- [x] Gene/sample normalization factors as GLM offsets for supplied-dispersion
  fixed Wald/LRT pipelines.
- [x] All-zero row expansion for the supplied-dispersion Wald/LRT pipelines, using
  `NaN` in internal matrices and `None` in result rows.
- [x] Add DESeq2-style optim fallback row routing for unstable IRLS rows,
  non-positive coefficient variances, and optional non-converged rows.
- [x] Add bounded pure-Rust optim fallback refits for routed fixed-dispersion
  IRLS rows, including optimized betas, SE/covariance, fitted means, and row
  log likelihoods.
- [x] Add optional DESeq2 reference-generation and skip-safe Rust test hooks
  for fixed-dispersion force-optim fallback rows.
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
- [x] Expose gene-wise dispersion iteration diagnostics (`dispGeneIter`) in
  high-level fit state and the `mcols(dds)`-style alias view.
- [x] Expose dispersion-stage `mcols(dds)`-style diagnostic aliases for
  `dispGeneEst`, `dispFit`, `dispersion`, `dispIter`, and `dispOutlier`.
- [x] Add stable present-column listing for the `mcols(dds)`-style diagnostic
  alias view.
- [x] Add typed data-frame assembly for present `mcols(dds)`-style diagnostic
  aliases.
- [x] Add TSV export for present `mcols(dds)`-style diagnostic aliases.
- [x] Keep matrix-valued GLM diagnostics (`mu`, hats, and reduced-model
  matrices) as explicit `DeseqFit` fields rather than `mcols(dds)` row
  metadata columns.
- [x] Add primitive result-table metadata carrier for test type, reported
  coefficient/contrast, column descriptions, p-value adjustment method, and
  independent-filtering metadata.
- [x] Carry Wald `lfcThreshold` and `altHypothesis` settings into high-level
  result-table metadata.
- [x] Add DESeq2-style comparison-aware descriptions for implemented result
  columns.
- [x] Add a typed data-frame view for implemented result rows, with row names,
  numeric/logical columns, and per-column metadata.
- [x] Factor shared all-zero row expansion helpers for compact GLM matrices,
  compact per-gene vectors, and full-length masked vectors.
- [x] Add initial linear-mu gene-wise dispersion foundation:
  `linearModelMu`, `roughDispEstimate`, `momentsDispEstimate`, bounded starts,
  Cox-Reid objective scoring, Armijo line search, and grid fallback.
- [x] Full result-table assembly with DESeq2-style metadata for implemented
  primitive result rows.
- [x] Optim fallback for non-converged or unstable fixed-dispersion IRLS rows.
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
  filtering, Gamma identity-link IRLS, and offline prediction reference
  coverage.
- [x] Mean dispersion trend:
  shared `100 * minDisp` viability gate, `10 * minDisp` filtered trimmed mean,
  constant fitted trend expansion, `FitType::Mean` builder dispatch, and
  offline DESeq2 reference coverage.
- [x] Initial local dispersion trend:
  `10 * minDisp` fit rule, base-mean weights, all-near-minimum floor behavior,
  pure-Rust adaptive local polynomial smoothing, `FitType::Local` builder
  dispatch, optional DESeq2 local-trend reference check, and explicit
  all-near-minimum local floor, mixed-threshold, and out-of-sample prediction
  fixtures.
- [x] Default `fit_dispersion_trend` dispatcher for implemented `FitType`
  values: parametric, local, and mean.
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
- [x] DESeq2-style weighted Cox-Reid threshold-subset objective, first
  derivative, and second derivative.
- [x] Initial local dispersion trend type.
- [ ] glmGamPoi dispersion trend type and exact local `locfit` edge-case parity.
- [x] DESeq2-shaped prior variance branch for residual df 1 through 3.
- [x] Public `estimate_dispersion_prior` stage wrapper over the implemented
  prior-variance estimator.
- [ ] Exact DESeq2 seeded Monte Carlo/R `loess` numerical identity for
  residual-df 1 through 3 prior variance.
- [x] Low-level observation-weighted MAP dispersion fitting.
- [x] Builder-level observation weights through the GLM-mu gene-wise
  dispersion, MAP, and native Wald path.
- [x] Default R reference generation and skip-safe Rust checks for unweighted
  GLM-mu mean-trend MAP/Wald/LRT intermediates.
- [x] Default R reference generation and skip-safe Rust checks for weighted
  GLM-mu dispersion/MAP/Wald intermediates.
- [x] Match DESeq2's `minmu`-floored stored means in the GLM-mu dispersion
  estimation intermediate.
- [x] Match DESeq2's weighted `fitDisp` row-indexing behavior during GLM-mu
  `fitidx` mean/dispersion alternation.
- [x] Match DESeq2's weighted Cox-Reid behavior for the non-linear-mu
  gene-wise path: observation weights multiply likelihood terms, while the
  Cox-Reid determinant uses the `weightThreshold` sample subset without
  multiplying by the weights.
- [x] Add default DESeq2-internal unweighted GLM-mu `useCR=TRUE` gene-wise
  reference checks for Cox-Reid dispersion and stored-mean parity.
- [x] Add default DESeq2-internal unweighted GLM-mu mean-trend MAP/Wald/LRT
  reference checks with Cox-Reid enabled through gene-wise and MAP dispersion
  fitting.
- [x] Add default DESeq2-internal weighted GLM-mu `useCR=TRUE` gene-wise
  reference checks for weighted Cox-Reid dispersion parity.
- [x] Add default DESeq2-internal weighted GLM-mu mean-trend MAP/Wald/LRT
  reference checks with Cox-Reid enabled through gene-wise and MAP dispersion
  fitting.
- [x] Assert result-table BH-adjusted p-values against DESeq2 for the matched
  unweighted and weighted GLM-mu mean-trend Wald/LRT branches.
- [x] Assert compact public result-row parity for matched GLM-mu Wald/LRT
  branches, covering LFC, SE, statistic, p-value, adjusted p-value,
  dispersion, and convergence.
- [x] Assert MAP, Wald, LRT, and compact public result-row parity for the
  current unweighted GLM-mu local-trend fixture.
- [x] Assert MAP, Wald, LRT, and compact public result-row parity for the
  current weighted GLM-mu local-trend fixture.
- [x] Assert MAP dispersion intermediate parity for the current unweighted
  GLM-mu Cox-Reid local-trend fixture.
- [x] Assert MAP dispersion intermediate parity for the current weighted
  GLM-mu Cox-Reid local-trend fixture.
- [x] Assert compact Wald/LRT result-row parity for the current unweighted and
  weighted GLM-mu Cox-Reid local-trend fixtures.
- [x] Assert detailed Wald/LRT intermediate parity for GLM-mu Cox-Reid
  local-trend fixtures beyond the compact public result rows.
- [x] Add native GLM-mu Wald contrast entry points that reuse MAP dispersions,
  including numeric, coefficient-name/list, and primitive factor-level
  contrast metadata paths.
- [ ] Complete broader weighted dispersion parity for DESeq2's non-linear-mu
  GLM mean fitting path beyond the current deterministic mean/local-trend
  fixtures.
- [ ] glmGamPoi MAP dispersion path.

## Phase 4: Full Wald Pipeline

- [x] Run size factors, base statistics, supplied dispersion, beta fitting,
  Wald tests, p-values, and adjusted p-values in a single inspectable state for
  one selected coefficient.
- [x] Add native dispersion estimation to that pipeline for the limited
  linear-mu, parametric-trend, deterministic-prior, no-weight MAP branch.
- [x] Generalize the limited native MAP/Wald branch to builder-selected
  parametric, local, or mean dispersion trends.
- [x] Generalize the limited native MAP/Wald branch to the GLM-mu mean-refit
  dispersion branch, including builder-level observation weights.
- [x] Wire `DeseqBuilder::fit()` to the implemented top-level GLM-mu Wald
  workflow, using the last design coefficient by default.
- [x] Add top-level Wald result-table workflow via
  `DeseqBuilder::fit_with_results()`.
- [x] Add top-level GLM-mu Wald result-table workflow for primitive numeric,
  named/list, and caller-supplied factor-level contrasts.
- [x] Add fit-only top-level GLM-mu Wald contrast helpers mirroring the
  result-table contrast routes.
- [ ] Generalize the native Wald pipeline to DESeq2's full dispersion and GLM
  fitting behavior.
- [ ] Keep mature R wrapper workflow exposure deferred until the Rust pipeline
  is complete enough to expose without DESeq2 runtime fallback.

## Phase 5: Results Compatibility

- [x] Primitive numeric Wald linear contrasts.
- [x] Result-row assembly for precomputed primitive numeric Wald contrasts.
- [x] TSV writer for assembled result tables, preserving DESeq2-shaped column
  order and R-style `NA` output for missing numeric or logical values.
- [x] Builder-level supplied-dispersion Wald pipeline for primitive numeric
  contrasts.
- [x] Primitive coefficient-name, positive/negative coefficient-list, and
  common factor-level contrast resolution against design coefficient names.
- [x] Stable result-table names and comparison labels for named primitive
  contrast specifications.
- [x] Numeric/expanded `contrastAllZero` behavior for primitive Wald contrasts,
  matching DESeq2's `contrastAllZeroNumeric` model-matrix selection rule.
- [x] Character/factor-level `contrastAllZero` behavior for primitive
  factor-level Wald contrasts with caller-supplied sample levels, matching
  DESeq2's `contrastAllZeroCharacter` selected-sample rule.
- [x] Native linear-mu and GLM-mu Wald contrast wrappers for primitive numeric,
  coefficient-name/list, and caller-supplied factor-level requests.
- [ ] Full DESeq2 contrast handling for `results(contrast=...)`, including
  colData/formula-aware factor-level semantics, complete coefficient-name
  cleanup, and contrast-aware Cook's/refit edge cases.
- [x] Initial Cook's distance and `maxCooks` diagnostics.
- [x] Cook's outlier p-value filtering.
- [x] Explicit primitive helper for DESeq2's two-group low-count Cook's
  heuristic, to be called only when caller or future-wrapper formula metadata
  establishes the one-factor two-level condition.
- [x] Primitive Cook's outlier count replacement transform: trimmed normalized
  means, size-factor/normalization-factor rescaling, integer truncation,
  replaceable-sample mask, and `replace` flags.
- [x] Replacement-refit planning metadata: replacement-count base metadata,
  `nrefit`, `refitReplace`, `newAllZero`, and post-refit `maxCooks` masking.
- [x] Limited Cook's replacement-refit path for GLM-mu native Wald and LRT:
  original fit, replacement counts, replacement-count refit with original size
  factors, `refitReplace` merge, `newAllZero` result clearing, and final
  filtering.
- [x] Limited Cook's replacement-refit path for GLM-mu native Wald contrasts,
  including primitive numeric, named/list, and caller-supplied factor-level
  contrast routes.
- [x] Top-level GLM-mu Wald/LRT result helpers for limited Cook's replacement
  refit, including default coefficient and primitive contrast result routes.
- [ ] Full Cook's outlier replacement behavior with DESeq2-style replacement
  refit for beta priors, wrapper object metadata, and remaining contrast edge
  cases.
- [ ] Full formula-aware outlier handling and future wrapper integration
  without DESeq2 runtime fallback.
- [x] Initial independent filtering.
- [x] R `stats::lowess`-shaped independent-filter threshold selection for the
  DESeq2 default theta grid, including fitted-curve metadata.
- [x] R `lowess` parity for dense custom theta grids where R's
  `delta` interpolation shortcut skips fitted points.
- [x] Initial fixed-dispersion LRT.
- [x] Limited native-dispersion LRT for current linear-mu and GLM-mu MAP
  dispersion branches.
- [x] Limited native-dispersion LRT dispatch through the initial local
  dispersion trend for current linear-mu and GLM-mu branches.
- [x] Wire `DeseqBuilder::fit_lrt()` to the implemented top-level GLM-mu LRT
  workflow, using the last full-design coefficient by default.
- [x] Add top-level LRT result-table workflow via
  `DeseqBuilder::fit_lrt_with_results()`.
- [x] DESeq2-internal native weighted GLM-mu LRT reference hook for
  the current mean-trend branch.
- [x] Default DESeq2-internal native weighted GLM-mu Wald/LRT reference checks
  for the current mean-trend branch.
- [ ] Full LRT parity with native dispersion reference outputs, exact local
  `locfit` behavior, glmGamPoi trends, optim fallback, and remaining edge
  cases.

## Phase 6: Secondary Features

- [x] Initial local dispersion trend.
- [x] Mean fit type.
- [ ] glmGamPoi-like mode if feasible.
- [x] `normTransform` log2 normalized-count-plus-one transform.
- [x] Mean-fit VST closed-form transform for normalized counts.
- [x] Parametric-trend VST closed-form transform for normalized counts.
- [x] Local-trend VST numerical-integration transform for normalized counts.
- [x] VST dispatch from an already-fitted `DispersionTrendFit`.
- [x] Local VST `mean(1 / sizeFactors)` and normalization-factor `xim`
  helpers.
- [x] Factor-aware VST dispatch helpers for size factors and normalization
  factors.
- [x] Store implemented fitted dispersion trends in `DeseqFit` and expose
  fit-level normalized-count, `normTransform`, and VST helpers.
- [x] Fit-level VST helper that applies an external fitted trend to the full
  count matrix, enabling the fast-subset trend/full-data transform split.
- [x] Add a fit-level `vst` alias and LRT workflow coverage for fitted-trend
  VST reuse.
- [x] Fast-VST deterministic subset index helper matching DESeq2's
  `baseMean > 5`, ordered row selection, and R-style rounding rule.
- [x] Public fast-VST default `nsub=1000` constant and default-size builder
  convenience methods.
- [x] Public fast-VST eligible-row count helper, plus fit-level eligibility
  query from stored `baseMean`.
- [x] Fast-VST normalized-count row-subset helper for the selected trend-fit
  subset.
- [x] Fast-VST aligned gene/sample matrix row-subset helper for normalization
  factors and observation weights.
- [x] Public `CountMatrix::select_rows` helper for fast-VST raw count subsets,
  preserving gene and sample names.
- [x] Public fast-VST row-aligned subset bundle for raw counts, normalized
  counts, optional normalization factors, optional observation weights, and
  original row indices.
- [x] Fast-VST subset metadata view with subset shape, original row indices,
  and factor/weight presence flags.
- [x] Fit-level fast-VST subset helper using stored `baseMean`, normalization
  factors, and preprocessed observation weights.
- [x] Builder-level GLM-mu fast-VST subset trend fitting that preserves
  full-data size factors and subset normalization factors.
- [x] Builder-level GLM-mu fast-VST transform applying the subset-fitted trend
  to the full normalized count matrix with subset diagnostics.
- [x] Named fast-VST GLM-mu output object for transformed counts, subset fit,
  and row-aligned subset diagnostics.
- [x] Explicit fast-VST output metadata view with transform shape, subset
  shape, subset row indices, and trend-fit shape.
- [x] Public fast-VST builder-level `nsub > 0` validation before branch-specific
  unsupported-feature checks.
- [x] Automatic GLM-mu VST helper that uses the deterministic fast subset when
  enough rows are eligible and otherwise falls back to a full-data Rust trend
  fit.
- [x] Automatic VST trend-source metadata recording fast-subset vs full-data
  trend fitting, requested `nsub`, and eligible-row count.
- [x] Accessor helpers for automatic VST trend-source metadata.
- [x] Automatic VST fast-subset trend fitting with observation weights carried
  through the deterministic row subset.
- [x] Full-data automatic VST reason metadata for insufficient eligible rows.
- [x] Stable string labels for automatic VST trend-source and full-data reason
  metadata.
- [x] Automatic VST output metadata view with source labels, subset sizing,
  transform shape, trend-fit row count, and optional fast-subset row count.
- [x] Automatic VST output metadata includes trend-fit sample count and
  original fast-subset row indices for parity diagnostics.
- [x] Explicit and automatic VST output metadata includes stable trend-fit type
  labels for parametric, local, and mean dispersion trends.
- [x] Blind automatic GLM-mu VST helper using a named intercept-only design,
  matching the implemented part of DESeq2's `blind=TRUE` transform shape.
- [ ] Full VST with automatic trend estimation, frozen dispersion-function
  reuse, fast-subset trend fitting, exact local `splinefun` parity, and object
  metadata.
- [ ] rlog implementation. The transform namespace currently exports an
  explicit `UnsupportedFeature` marker for callers that probe support.
- [ ] Mature R package wrapper backed only by Rust, with no fallback to
  R/Bioconductor DESeq2.
- [ ] Mature CLI.

## Engineering Rules

- Keep unsupported stages explicit with `DeseqError::UnsupportedFeature`.
- Prefer strict, deterministic behavior before adding fast paths.
- Preserve row-major gene-contiguous storage in Rust.
- Add hand-computable tests before golden-reference tests.
- Keep runtime statistical computation in Rust. Future wrappers may prepare
  inputs and present outputs, but must not call DESeq2 as a fallback.
- Document every intentional deviation from DESeq2 semantics.
