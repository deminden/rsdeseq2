# Compatibility

For a higher-level comparison against original DESeq2, see
[`deseq2-gap-analysis.md`](deseq2-gap-analysis.md).

## Numerical Parity Snapshot

The strongest current matches are stage-level primitives with generated
DESeq2 fixtures or direct formula checks. They are suitable for
apples-to-apples validation and benchmarking, but they are not a claim of full
`DESeq()` workflow parity.

| area | current numerical status | evidence |
| --- | --- | --- |
| Size factors and normalized counts | Matches DESeq2 `ratio` and `poscounts` behavior for covered fixtures, including supplied geometric means and control genes. | Hand tests, generated DESeq2 normalization fixtures, and real/synthetic CLI benchmark diffs. |
| Base row metadata | Matches `baseMean`, `baseVar`, `allZero`, normalization-factor metadata, and weighted base metadata for implemented inputs. | Generated DESeq2 metadata fixtures and unit tests. |
| Negative-binomial likelihood/deviance | Matches DESeq2's `mu`/dispersion parameterization and `-2 * logLike` convention. | Hand-formula and fixed-dispersion GLM fixture checks. |
| Fixed-dispersion GLM | Matches implemented `fitNbinomGLMs` fields for supplied dispersions: betas, SEs/covariance, fitted means, hats, log likelihood, Wald/LRT, weighted paths, and forced optim fallback fixtures. | Optional DESeq2-internal Wald/LRT/Cook's/optim reference tests plus deterministic Rust tests. |
| Beta prior primitives | Implements DESeq2-shaped beta-prior variance and fixed-dispersion refit math, including Hmisc-style weighted quantiles and log2-to-natural ridge conversion. | Source-matched formula tests; expanded-model averaging and high-level plumbing still need fixture coverage. |
| Dispersion trend and MAP pieces | Matches or closely tracks parametric/mean trend fixtures, initial local-trend fixtures including a single-usable-row edge case, prior variance, MAP shrinkage, unweighted GLM-mu Cox-Reid mean MAP/Wald/LRT and local MAP, unweighted GLM-mu mean and local MAP/Wald/LRT, weighted GLM-mu Cox-Reid mean MAP/Wald/LRT, weighted GLM-mu local MAP/Wald/LRT, and weighted GLM-mu deterministic anchors. | Generated DESeq2 trend/prior/MAP/GLM-mu fixtures and finite-difference objective tests. |
| Results, Cook's, filtering | Matches implemented result-table assembly, including DESeq2-shaped result rows and BH-adjusted p-values for the matched GLM-mu Wald/LRT fixture branches, Cook's distance/masking/replacement planning, selected replacement-refit paths, and independent-filtering lowess fixtures. | Unit tests plus generated Cook's, GLM-mu result-row, and independent-filtering fixtures. |
| Transform primitives | Matches closed-form `normTransform`, mean VST, parametric VST, deterministic fast-subset selection, and implemented local numerical-integration helpers. | Formula tests and stage-level dispatch tests; full Bioconductor object workflow remains future work. |

## Implemented

- Non-negative integer count matrix represented as genes x samples.
- Row-major storage for per-gene operations.
- DESeq2-like `ratio` size-factor estimation.
- DESeq2-style `poscounts` size-factor estimation.
- Optional caller-supplied size factors.
- Optional supplied geometric means.
- Optional control-gene subset by row index.
- Size-factor normalized counts, plus raw and normalized count-matrix TSV
  export with gene and sample labels.
- Gene/sample normalization factors that preempt size factors for normalized
  counts, base row metadata, supplied-dispersion fixed Wald/LRT GLM offsets,
  and the current native linear-mu dispersion/Wald subset, plus genes x samples
  TSV export for normalization-factor matrices.
- `baseMean`.
- `baseVar` sample variance of normalized counts.
- DESeq2-shaped TSV export for implemented early row metadata:
  `baseMean`, `baseVar`, and `allZero`.
- DESeq2-style weighted base metadata helpers that multiply normalized counts
  by observation weights before ordinary row means and variances.
- Builder-owned observation weights with design-aware row-normalization,
  `weights_fail` fit-state output, and all-zero-style skipping of failed rows
  in supplied-dispersion fixed Wald/LRT pipelines.
- `allZero` raw-count row flags.
- Benjamini-Hochberg adjusted p-values with missing-value preservation.
- Negative-binomial log PMF and row log-likelihood using DESeq2's
  `mu`/dispersion parameterization.
- Intercept-only fixed-dispersion NB GLM shortcut with beta estimates, fitted
  means, log likelihood, DESeq2-style full deviance, standard errors, iteration
  flags, and hat diagonals.
- Initial fixed-dispersion IRLS for supplied design matrices, with optional
  observation weights, scalar or per-coefficient natural-log-scale ridge
  values, and selectable normal-equation or DESeq2-style augmented QR solvers.
  Inspectable GLM fit state includes full/reduced log likelihoods,
  DESeq2-style full deviance, beta convergence, and beta iteration counts.
- DESeq2-style optim fallback row routing for unstable IRLS rows, non-positive
  coefficient variances, and optional non-converged rows, followed by bounded
  pure-Rust refits that refresh beta estimates, SE/covariance, fitted means,
  row log likelihoods, and convergence flags.
- DESeq2-style beta-prior variance estimation for primitive MLE beta matrices:
  unweighted type-7 quantile matching, DESeq2-vendored Hmisc weighted
  upper-quantile matching using `1 / (1 / baseMean + dispFit)` row weights,
  finite-beta filtering, and wide intercept priors.
- Primitive beta-prior fixed-dispersion GLM refits from supplied or estimated
  log2-scale beta-prior variances, including DESeq2's conversion to
  natural-log-scale ridge values before IRLS.
- DESeq2-style observation-weight preprocessing helper with row-max
  normalization, weighted design-rank checks, thresholded Cox-Reid sub-design
  checks, and `weights_fail` flags.
- DESeq2-style full-rank design matrix checks for supplied-dispersion Wald/LRT
  pipelines and native dispersion stages.
- DESeq2-shaped fit diagnostics expose stored fitted dispersion trend type
  labels and implemented dispersion-stage aliases, including `dispGeneEst`,
  `dispFit`, `dispersion`, `dispIter`, and `dispOutlier`, plus stable present
  column names, typed data-frame assembly, and TSV export for table exporters.
- Initial linear-mu gene-wise dispersion estimator with rough and moments
  dispersion starts.
- Initial GLM-mu gene-wise dispersion estimator with rough/moments
  starts, fixed-dispersion IRLS mean refits, fixed-mean dispersion
  optimization, `niter`, DESeq2's `.05` log-dispersion `fitidx` update rule,
  and optional builder observation weights.
- Fixed-mean Cox-Reid dispersion objective plus first and second derivatives,
  including DESeq2-style observation-weighted variants where weights multiply
  likelihood terms and threshold the Cox-Reid design subset.
- DESeq2-style log-dispersion prior objective plus first and second derivative
  terms for MAP dispersion fitting.
- DESeq2-shaped Armijo line-search optimizer for fixed-mean gene-wise
  dispersion estimates, with grid fallback for non-converged rows and
  prior-aware and observation-weighted optimizer entry points.
- Parametric dispersion trend fitting with DESeq2's
  `asymptDisp + extraPois / mean` form, `dispGeneEst > 100 * minDisp`
  row-selection rule, robust residual screen, and Gamma identity-link IRLS.
- Mean dispersion trend fitting with DESeq2's preliminary
  `dispGeneEst > 100 * minDisp` viability gate and the constant trimmed-mean
  fit over estimates above `10 * minDisp`.
- `normTransform` and mean-fit, parametric, and local numerical-integration VST
  primitives for normalized-count matrices, including dispatch from an
  already-fitted `DispersionTrendFit` and DESeq2-shaped local VST size-factor
  summary helpers. Factor-aware VST dispatch helpers now accept ordinary size
  factors or normalization factors directly for the local branch. `DeseqFit`
  retains the implemented fitted trend object and can produce normalized
  counts, `normTransform`, and VST output for its source count matrix, with a
  short `vst` method alias for the transform. It can also apply a caller-supplied
  fitted trend to its full count matrix, which supports the fast-VST split
  between full-data normalization and subset-fitted trend estimation. The
  deterministic row-subset helper used by DESeq2's fast `vst()` wrapper
  includes the DESeq2 default `nsub=1000` constant and explicit-size lower-level
  building blocks: eligible-row counting under the same `baseMean > 5` rule,
  normalized-count subsetting, aligned gene/sample matrix subsetting, raw count
  `CountMatrix` row selection with name preservation, and a row-aligned subset
  bundle for counts, normalized counts, optional normalization factors,
  optional observation weights, and original row indices. `DeseqFit` exposes
  the eligible-row count and that bundle using its stored `baseMean`,
  normalization factors, and preprocessed observation weights. The builder can
  fit the selected GLM-mu dispersion trend on the deterministic fast-VST subset
  while preserving full-data size factors and subset normalization factors;
  the paired fast-VST helper applies that subset-fitted trend to the full
  normalized count matrix through size-factor or normalization-factor dispatch
  and returns a named output with the transformed matrix, subset fit, and
  subset diagnostics, plus a metadata view with full transform shape, subset
  shape, subset row indices, trend-fit shape, and stable trend-fit type label.
  The automatic GLM-mu VST helper uses the fast subset when
  enough rows are eligible and otherwise fits the selected Rust trend on all
  rows, with trend-source metadata for the chosen path, requested `nsub`, and
  eligible-row count. Full-data trend metadata distinguishes too-few-eligible
  rows from observation-weighted input, and exposes stable labels for source
  and reason fields. The output metadata view packages those labels with
  transform shape, trend-fit shape, trend-fit type label, optional fast-subset
  row count, and optional original fast-subset row indices.
  When observation weights are present,
  automatic VST uses the full-data trend path because observation-weighted
  fast-VST trend fitting is still held back until the weight preprocessing
  semantics are wired without double-normalizing rows. The blind automatic
  helper uses a named intercept-only design for the same decision.
- Initial pure-Rust local dispersion trend fitting with DESeq2's
  `dispGeneEst >= 10 * minDisp` local-fit rule, base-mean weights,
  all-near-minimum floor behavior, builder dispatch, and a small offline
  local-trend reference check.
- Deterministic dispersion prior variance estimation, including R-compatible
  MAD scaling, `trigamma((m - p) / 2)` sampling-variance subtraction for
  residual degrees of freedom above 3, saturated designs, low-df histogram/KL
  matching, and the `0.25` floor.
- Initial MAP dispersion fitting with DESeq2's `dispInit` rule,
  `log(dispFit)` prior means, prior-aware line search, grid fallback,
  optional observation weights, `dispMAP`, `dispOutlier`, and final
  dispersion values for the linear-mu and GLM-mu branches.
- Limited native Wald pipeline for the current linear-mu no-weight and GLM-mu
  optionally weighted, deterministic-prior MAP dispersion subsets with
  parametric, local, or mean dispersion trends.
- Default coefficient-level Wald statistic and standard-Normal p-value.
- DESeq2-style Wald t p-values with residual, scalar, or per-gene degrees of
  freedom.
- Log2-scale beta covariance matrices exposed in `DeseqFit` for implemented GLM
  fits and primitive numeric Wald linear contrasts using `c' beta` and
  `sqrt(c' Sigma c)`.
- Result-row assembly with BH adjustment for precomputed primitive numeric
  Wald contrasts.
- Primitive coefficient-name, positive/negative coefficient-list, and common
  factor-level contrast resolution against design coefficient names, with
  stable result-table names and comparison labels for named contrast specs.
- Native linear-mu and GLM-mu Wald contrast entry points reuse the implemented
  MAP dispersion paths, then run numeric, named/list, or caller-supplied
  factor-level contrast result assembly.
- Numeric/expanded DESeq2-style `contrastAllZero` handling for primitive Wald
  contrasts: selected samples are inferred from `modelMatrix %*%
  contrastBinary`, and eligible rows are assigned LFC/stat zero and p-value
  one before result-table adjustment.
- Character/factor-level DESeq2-style `contrastAllZero` handling for primitive
  factor-level Wald contrasts when the caller supplies sample levels.
- Selected-coefficient Wald LFC-threshold alternatives from `results()`:
  `greaterAbs`, `greaterAbs2014`, `greaterAbsUPSHOT` without t p-values,
  `lessAbs`, `greater`, and `less`, with threshold and alternative metadata
  carried by high-level Wald result tables.
- Selected-coefficient Wald result rows with BH-adjusted p-values.
- DESeq2-style result-column metadata descriptions for the reported
  coefficient, primitive contrast, or LRT model comparison, with public
  effect/test description label helpers and table-level scalar metadata for
  wrapper and exporter parity.
- Typed DESeq2-shaped result-table view for implemented primitive rows, with
  row names, numeric/logical columns, and column metadata for wrapper/file
  output, plus regular and `results(tidy = TRUE)`-style TSV result export,
  `mcols(res)`-style result-column metadata export, and table-level result
  metadata export.
- Supplied-dispersion fixed-dispersion Wald pipeline for one coefficient and
  primitive numeric contrasts.
- Supplied-dispersion fixed-dispersion LRT pipeline for full vs reduced design
  matrices, with full and reduced log-likelihood, beta convergence, and
  iteration diagnostics in `DeseqFit`.
- Limited native-dispersion LRT pipelines for the implemented linear-mu and
  GLM-mu MAP dispersion branches, with the full design used for dispersion
  estimation before full-vs-reduced testing.
- All-zero row expansion for the supplied-dispersion Wald/LRT pipelines and
  limited native Wald/LRT paths, using missing numeric outputs for skipped
  all-zero rows.
- Cook's distance matrix for the supplied-dispersion and limited native Wald
  pipelines.
- `maxCooks` over samples in model-matrix cells with at least three replicates.
- Cook's cutoff p-value masking with BH recomputation for result rows.
- Explicit primitive helper for DESeq2's two-group low-count Cook's heuristic:
  rows above cutoff are spared when at least three counts are larger than the
  count in the sample with maximum Cook's distance. This helper is not applied
  automatically because the Rust core does not own R formula/colData metadata.
- Primitive Cook's outlier count replacement transform from `replaceOutliers`:
  trimmed normalized means, size-factor or normalization-factor rescaling,
  integer truncation, replaceable-sample masks, and per-gene `replace` flags.
- Cook's replacement-refit planning metadata: replacement-count base metadata,
  `nrefit`, `refitReplace`, `newAllZero`, and DESeq2-style post-refit
  `maxCooks` masking over nonreplaceable samples.
- Limited Cook's replacement-refit path for the implemented GLM-mu native Wald
  and LRT branches, merging only `refitReplace` rows and preserving original
  size factors.
- Limited Cook's replacement-refit path for GLM-mu native Wald contrasts,
  including primitive numeric, named/list, and caller-supplied factor-level
  contrast routes.
- Base-mean independent filtering with filtered BH adjustment metadata and an
  R `stats::lowess`-shaped rejection-curve smoother for DESeq2's default
  threshold-selection grid, plus a paired `filterNumRej`-shaped
  theta/rejection-count view, a paired `lo.fit`-shaped lowess view, scalar
  metadata entries, and TSV exporters for those metadata tables.

## Reference Points

Implemented normalization behavior is based on DESeq2's public R API and source
comments around `estimateSizeFactorsForMatrix` in `R/core.R`, plus
`estimateSizeFactors.DESeqDataSet` and `normalizationFactors` documentation in
`R/methods.R`. Normalized counts follow `counts.DESeqDataSet` in
`R/methods.R`: `normalizationFactors` are used first and size factors are the
fallback. Fixed-dispersion GLM offsets follow `getSizeOrNormFactors` in
`R/core.R` and `fitNbinomGLMs`. Full-rank model-matrix checks follow
`checkFullRank` in `R/core.R`, where `qr(modelMatrix)$rank < ncol(modelMatrix)`
stops GLM fitting. Base means and variances follow
`getBaseMeansAndVariances` in `R/core.R`: without weights these are ordinary
row summaries of normalized counts; with observation weights, DESeq2 first
multiplies normalized counts by the weights assay and then computes ordinary
row summaries. Cook's
diagnostics follow `calculateCooksDistance`,
`robustMethodOfMomentsDisp`, `nOrMoreInCell`, and `recordMaxCooks` in
`R/core.R`. The primitive outlier replacement transform follows
`replaceOutliers`: replacement candidates are built from trimmed means of
normalized counts and rescaled by size factors or normalization factors, while
only replaceable outlier cells are written into the transformed count matrix.
The refit-planning helper follows `refitWithoutOutliers` bookkeeping for
`replace`, `newAllZero`, `refitReplace`, and the post-refit `maxCooks` rule,
and the limited GLM-mu native Wald replacement-refit path reruns the implemented
dispersion/MAP/Wald stages on replacement counts with original size factors.
The limited GLM-mu native LRT replacement-refit path uses the same replacement
bookkeeping and reruns the implemented dispersion/MAP/LRT stages with original
size factors before merging refit rows. This is not yet a full Bioconductor
`refitWithoutOutliers` implementation.
The limited native Wald pipeline follows the high-level DESeq2
stage order in `DESeq`: size factors, dispersion estimation, then
`nbinomWaldTest`, but only for the currently implemented linear-mu no-weight
and GLM-mu optionally weighted dispersion branches. Cook's cutoff masking
follows the default `results()` path in
`R/results.R`, where `maxCooks > qf(.99, p, m - p)` makes the p-value missing
before p-value adjustment. Independent filtering follows `pvalueAdjustment`
and `filtered_p` in `R/results.R`: `baseMean` is the default filter statistic,
candidate cutoffs are filter quantiles, and BH adjustment is recomputed after
filtering. The local ignored clone at `external/DESeq2` is for inspection only
and is not vendored into this package.

Wald t p-values follow the `useT=TRUE` branch in `nbinomWaldTest`: supplied
or residual degrees of freedom are used with `pt`, and non-positive df values
become missing p-values. Thresholded selected-coefficient Wald tests follow
the `lfcThreshold`/`altHypothesis` branch in `R/results.R`, operating on log2
fold changes and log2 standard errors. LRT behavior follows `nbinomLRT` in
`R/core.R`: full and reduced models are fit on non-all-zero rows, the statistic
is `2 * (logLik_full - logLik_reduced)`, and p-values use a chi-square upper
tail with degrees of freedom `ncol(full) - ncol(reduced)`. DESeq2 stores
`fullBetaConv`, `reducedBetaConv`, full `betaIter`, and full-model `deviance =
-2 * logLike`; `rsdeseq2` exposes the matching full-model fields plus
reduced-model convergence and iteration diagnostics in `DeseqFit`.

The current gene-wise dispersion foundation follows the early
`estimateDispersionsGeneEst` anchors in `R/core.R`: `roughDispEstimate`,
`momentsDispEstimate`, `linearModelMu`, and the two-pass fallback grid shape in
`fitDispGridWrapper`. It implements the Cox-Reid log determinant objective term
for unweighted fits, matches the unweighted GLM-mu Cox-Reid gene-wise fixture,
has low-level weighted objective/derivative/curvature variants, and follows the
Armijo line-search control flow from `fitDisp` for no-prior estimates. The MAP
prior term follows
`fitDisp`/`fitDispGridWrapper`: a normal prior kernel on `log(alpha)` with
mean `log(dispFit)` and variance `dispPriorVar`, plus the corresponding
first and second derivatives. Parametric dispersion trends follow
`parametricDispersionFit` in `R/core.R`, including the
`dispGeneEst > 100 * minDisp` fit rule and the robust residual screen before
fitting the Gamma identity-link model. Mean dispersion trends follow
`estimateDispersionsFit(fitType="mean")`: the shared `100 * minDisp` viability
gate is checked first, then the constant trend is the trimmed mean of estimates
above `10 * minDisp`. The deterministic prior-variance branch
follows `estimateDispersionsPriorVar`: MAD-squared log residual variance,
optional `trigamma((m - p) / 2)` subtraction, low-df histogram/KL matching,
and the `0.25` floor.
The MAP stage follows `estimateDispersionsMAP(type="DESeq2")`, including the
initialization rule, prior-aware optimizer, optional low-level observation
weights, grid fallback for non-converged rows, bounding of `dispMAP`, and
high-side dispersion outlier replacement. The high-level native linear-mu
pipeline intentionally remains no-weight because DESeq2 switches away from the
linear-mu branch when weights are present.

The R reference generator writes both full DESeq2 outputs and narrower
fixed-dispersion references. The latter use `DESeq2:::fitNbinomGLMs` with
caller-supplied dispersions, default `1e-6` beta ridge, `useQR=FALSE`, and
`useOptim=FALSE`; these are the current apples-to-apples references for the
Rust fixed-dispersion Wald/LRT implementation. That default reference set now
includes fixed-dispersion beta, SE, Wald/LRT statistics, log likelihoods,
fitted means, hat diagonals, and Cook's distances. It also writes weighted base
metadata using `getBaseMeansAndVariances` and `getAndCheckWeights`, weighted
fixed-dispersion Wald/LRT references with `weightsFail` rows expanded as
missing, unweighted GLM-mu Cox-Reid mean-trend MAP/Wald/LRT and mean-trend
MAP/Wald/LRT intermediates, plus weighted GLM-mu Cox-Reid mean-trend MAP and
mean-trend and local-trend dispersion/MAP/Wald/LRT intermediate references for
the current native branch. The matched GLM-mu Wald/LRT branches also write compact
DESeq2-shaped result-table fixtures and assert result-row beta, SE, statistic,
p-value, adjusted p-value, dispersion, and convergence parity.
The unweighted GLM-mu local-trend fixture covers MAP, Wald, LRT, and compact
Wald/LRT result rows for the tiny-data edge case where only one row is usable
for local fitting; Rust now follows the corresponding constant local fit shape
instead of failing the trend evaluation.
The weighted GLM-mu local-trend fixture covers the same MAP/Wald/LRT/result-row
surface with DESeq2's `weightsFail` expansion semantics.
The unweighted GLM-mu Cox-Reid local-trend fixture currently covers the MAP
dispersion intermediate; Wald/LRT result-row coverage for that exact
combination remains future work.

The fixed-dispersion IRLS path now includes a bounded pure-Rust optim fallback
for routed rows. The current checked reference fixtures still use
`useOptim=FALSE`; the reference generator also has an optional
`forceOptim=TRUE` fixture and skip-safe Rust comparison hook for validating the
bounded fallback where DESeq2 is installed locally.

## Missing

- Mature wrapper-facing interface around the Rust core.
- R formula parsing in Rust.
- Full DESeq2 dispersion estimation, including broader weighted dispersion
  edge-case parity, exact `locfit` local-trend numerical identity, glmGamPoi trend type,
  and production-ready end-to-end dispersion parity.
- High-level propagation of observation-weight `weights_fail` flags through
  complete DESeq2-like builder and future wrapper workflows.
- Direct weighted low-level `fitNbinomGLMs` parity for rows that DESeq2 marks
  `weightsFail` but still returns ridge-stabilized coefficients for when that
  internal function is called directly.
- Broader DESeq2 parity fixtures for unstable or non-converged rows routed
  through the bounded optim fallback.
- Full DESeq2 `results(contrast=...)` colData/formula-aware factor-level
  semantics, complete coefficient-name cleanup, and contrast-aware Cook's/refit
  edge cases.
- Automatic formula-aware application of the two-group low-count Cook's
  heuristic from high-level wrappers.
- Full Cook's outlier replacement behavior for beta priors, Bioconductor assay
  preservation, and all remaining DESeq2 edge cases.
- General Wald and LRT tests with native dispersion estimation beyond the
  current limited linear-mu/GLM-mu branches and without generated DESeq2
  references for all native LRT intermediates.
- Full Bioconductor result-object metadata. Lightweight `DeseqFit` and
  `DeseqResults` metadata views exist for implemented diagnostics, result
  columns, test type, p-value adjustment method, Wald threshold settings,
  comparison-aware column descriptions, effect/test description labels,
  table-level scalar metadata, typed result-table output, regular and tidy TSV
  result export, result-metadata export, and independent-filtering metadata
  including paired `filterNumRej`- and `lo.fit`-shaped views plus scalar
  metadata entries.
  Shared internal all-zero row expansion helpers exist for compact GLM outputs
  and full-length masked vectors.
- High-level R-style contrast handling beyond primitive coefficient-name
  resolution.
- Full high-level VST object workflow, exact DESeq2 `splinefun` behavior for
  local VST, rlog, and lfcShrink-compatible hooks. Mean-fit, parametric, local,
  fast-subset selection/subsetting, automatic GLM-mu VST, blind automatic VST,
  and fit-state VST helpers are available for normalized-count matrices and
  implemented fitted trends. The transform namespace exports rlog as an
  explicit unsupported marker until the regularized-log workflow is implemented.

## Known Deviations

The Rust core accepts primitive matrices and explicit options. Formula
semantics and model-matrix generation can be handled by a future wrapper.

Gene/sample normalization factors are supported for normalized-count metadata,
supplied-dispersion fixed Wald/LRT pipelines, and the current native linear-mu
dispersion path. The native subset follows DESeq2's
`linearModelMuNormalized` and `momentsDispEstimate` branches: fitted raw means
use `linearModelMu(counts(normalized=TRUE), X) * getSizeOrNormFactors(object)`,
and moments starts use `mean(1 / colMeans(normalizationFactors))` on the
non-all-zero row subset.

The size-factor implementation follows the current DESeq2
`estimateSizeFactorsForMatrix` shape for implemented options: geometric means
are represented as log geometric means, sample locations are computed on
log-ratio values and exponentiated, and `poscounts` sums positive log counts but
divides by the total number of samples. This is documented because a phrase
like "mean over positive counts" can otherwise be interpreted as dividing only
by the number of positive samples.

Independent-filter threshold selection follows the DESeq2 `results()` path with
filtered BH columns, rejection counts at each theta, and
`stats::lowess(numRej ~ theta, f=1/5)`-shaped smoothing for the default
50-point theta grid and dense custom theta grids where R's `delta`
interpolation shortcut skips fitted points. The implementation is tested
against R-generated fixtures for both cases.

The current gene-wise dispersion estimators have the Cox-Reid objective and
Armijo line-search shape, including DESeq2's weighted Cox-Reid
`weightThreshold` subset rule. The linear-mu branch is wired through trend,
MAP, and the limited native Wald path. The newer GLM-mu branch performs
mean/dispersion alternation using fixed-dispersion IRLS, can consume
preprocessed builder observation weights, and is wired through
parametric/local/mean trend fitting, MAP shrinkage, and the limited native
Wald path.
DESeq2 switches away from linear-mu fitting when observation weights are
present. `rsdeseq2` now has the log-prior objective, first and second
derivatives, weighted low-level objective variants, parametric and mean trend
foundations, deterministic prior-variance branches including low residual
degrees of freedom, and initial MAP stages, but not exact seeded Monte
Carlo/loess identity for the low-df prior-variance branch or complete
high-level weighted dispersion parity.
Parity tests
should therefore compare rough/moments starts, objective values, prior
objective values, trend coefficients, `dispFit`, `dispPriorVar`, `dispMAP`,
`dispOutlier`, line-search diagnostics including `last_dlp`/`last_d2lp`, and
grid fallback behavior before full result-table comparison.

Rust `controlGenes` can be provided as zero-based indices or a logical mask.
DESeq2 numeric control genes are one-based because they are R indices.

## Planned Parity Thresholds

Initial stage parity should use exact or near-machine-precision thresholds for
size factors, normalized counts, base means, and BH adjusted p-values.

Future GLM and dispersion parity should compare intermediate fields first:

- `dispGeneEst`
- `dispGeneIter` presence and skipped-row shape; exact counts are diagnostic
  because equivalent Armijo paths can take different numbers of steps.
- `dispFit`
- final dispersion
- beta estimates
- beta standard errors
- Wald or LRT statistics
- p-values and adjusted p-values
- Cook's distances
- `maxCooks`
- independent-filtering flags

Strict mode should stay deterministic. Fast mode can later allow controlled
floating-point drift after strict-mode behavior is established.
