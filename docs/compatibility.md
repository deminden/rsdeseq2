# Compatibility

For a higher-level comparison against original DESeq2, see
[`deseq2-gap-analysis.md`](deseq2-gap-analysis.md).

## Numerical Parity Snapshot

The strongest matches are stage-level primitives with generated
DESeq2 fixtures or direct formula checks. They support matched validation and
benchmarking, but they are not a claim of full
`DESeq()` workflow parity.

| area | numerical status | evidence |
| --- | --- | --- |
| Size factors and normalized counts | Matches DESeq2 `ratio` and `poscounts` behavior for covered fixtures, including supplied geometric means and control genes. | Hand tests, generated DESeq2 normalization fixtures, and real/synthetic CLI benchmark diffs. |
| Base row metadata | Matches `baseMean`, `baseVar`, `allZero`, normalization-factor metadata, and weighted base metadata for implemented inputs. | Generated DESeq2 metadata fixtures and unit tests. |
| Negative-binomial likelihood/deviance | Matches DESeq2's `mu`/dispersion parameterization and `-2 * logLike` convention. | Hand-formula and fixed-dispersion GLM fixture checks. |
| Fixed-dispersion GLM | Matches implemented `fitNbinomGLMs` fields for supplied dispersions: betas, SEs/covariance, fitted means, hats, log likelihood, Wald/LRT, weighted paths, and forced optim fallback fixtures. Optim fallback rows use an R-compatible callback solution with an analytic stability check that does not use reference data at runtime, while preserving DESeq2's pre-optimizer IRLS hats for downstream Cook's decisions. | Optional DESeq2-internal Wald/LRT/Cook's/optim reference tests, deterministic Rust tests, and the versioned 100-row high-error benchmark: v0.2.5 median `6.064e-10` (v0.2.4 `1.464e-4`; 241,402x lower) and maximum `1.526e-3` (v0.2.4 `3.094e-3`; 50.68% lower). |
| Beta prior primitives | Implements DESeq2-shaped beta-prior variance and fixed-dispersion refit math, including Hmisc-style weighted quantiles, log2-to-natural ridge conversion, primitive one-factor and additive expanded design construction for categorical factors, numeric covariates, primitive pairwise interactions, formula-only higher-order interactions, formula term subtraction, primitive numeric transforms, raw and orthogonal polynomial formula transforms, nested additive parenthesized groups, formula offset extraction and integration, owned model-frame formula inputs with factor reference inference/override, builder-stored model-frame formula design helpers, primitive expanded-design refits, grouped collapse, Wald result assembly/workflows from collapsed prior fits with size-factor, normalization-factor, and observation-weight inputs, and primitive expanded beta-prior Wald Cook's replacement refits for selected coefficients and numeric contrasts with size-factor or normalization-factor offsets. One-factor, additive, supported-formula, and owned-model-frame helpers build the expanded design internally for replacement-refit workflows, and the native CLI exposes primitive supplied-dispersion expanded, one-factor, and additive categorical beta-prior Wald paths with normalization-factor and replacement-sidecar coverage. | Source-matched formula tests plus DESeq2 `estimateBetaPriorVar`, beta-prior refit, combined estimated-prior refit fixture checks, Rust expanded-model workflow tests, model-frame wrapper tests, and CLI smoke tests; splines, arbitrary R expressions, and full R-compatible formula parsing are unsupported. |
| Dispersion trend and MAP pieces | Matches or closely approximates parametric/mean trend fixtures, pure-Rust locfit-compatible local trend fixtures including a single-usable-row edge case and real-data fitted-value parity, prior variance, MAP shrinkage, unweighted GLM-mu Cox-Reid mean MAP/Wald/LRT and local MAP/result rows, unweighted GLM-mu mean and local MAP/Wald/LRT, weighted GLM-mu Cox-Reid mean MAP/Wald/LRT and local MAP/result rows, weighted GLM-mu local MAP/Wald/LRT, and weighted GLM-mu deterministic references. | Generated DESeq2 trend/prior/MAP/GLM-mu fixtures and finite-difference objective tests. |
| Results, Cook's, filtering | Matches implemented result-table assembly, including DESeq2-shaped result rows and BH-adjusted p-values for the matched GLM-mu Wald/LRT fixture branches, R-cleaned coefficient/list/factor-level aliases, R-cleaned factor-level result names, explicit and builder-stored owned-model-frame character contrast resolution for supported fixed/native Wald/LRT routes, resolved numeric contrast metadata for implemented contrast result and replacement/refit paths, top-level formula-built Wald designs and LRT full/reduced designs for supported model-frame formulas, Cook's distance/masking/replacement planning, scalar replacement/refit metadata summaries, selected formula-built and explicit-design replacement-refit paths, and independent-filtering lowess fixtures. | Unit tests plus generated Cook's, GLM-mu result-row, and independent-filtering fixtures. |
| Transform primitives | Matches closed-form `normTransform`, mean VST, parametric VST, deterministic unweighted fast-subset selection, implemented local numerical-integration helpers, and the low-level rlog sample-effect ridge-GLM primitive with explicit dispersions plus rlog sample-prior estimation from normalized counts. The weighted fast-VST selector uses Rust's weighted `baseMean` and is not claimed to match DESeq2 subset indices. Convenience rlog helpers compose prior estimation with size-factor or normalization-factor fitting when earlier-stage summaries are supplied; fit-state and builder-level design-aware/blind GLM-mu rlog dispatch are available after MAP dispersions are present and skip/re-expand all-zero rows. Rlog output metadata records shape, prior variance, offset mode, design mode, and retained fit diagnostics for the builder path. | Formula tests, stage-level dispatch tests, builder rlog tests, all-zero rlog expansion tests, and CLI rlog tests; the full Bioconductor object workflow is unsupported. |

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
  and the supported native linear-mu dispersion/Wald subset, plus genes x samples
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
- Primitive expanded-model beta-prior coefficient averaging, covariance
  propagation, and averaged numerator/denominator contrast-vector construction
  for already-built expanded coefficient matrices.
- Primitive expanded-model fit collapse into a standard `NbinomGlmFit` result
  with collapsed betas, standard errors, covariance, and reported design.
- Primitive expanded-design beta-prior refit workflow for already-built
  matrices: expanded MLE fit, prior variance estimation, expanded prior refit,
  and collapse to the standard result type.
- DESeq2-shaped Wald coefficient and numeric-contrast result assembly from
  collapsed expanded-model fits.
- DESeq2-shaped Wald coefficient and numeric-contrast result assembly directly
  from expanded beta-prior refit outputs.
- Primitive expanded beta-prior fit-and-Wald-results workflow helpers for
  selected coefficients and numeric contrasts, including size-factor,
  normalization-factor, and optional observation-weight inputs.
- Primitive size-factor and normalization-factor expanded beta-prior Wald Cook's replacement-refit
  helpers for selected coefficients and numeric contrasts. Cook's distances
  are calculated on the collapsed reported design; replacement-count refits
  reuse the original offsets and supplied dispersions.
- Primitive one-factor expanded design construction from sample labels, with
  expanded intercept-plus-level columns, treatment-style reported design, and
  coefficient groups.
- One-factor expanded beta-prior fit-and-Wald-results helpers that build the
  design internally, then run coefficient or numeric-contrast workflows with
  size-factor, normalization-factor, and optional observation-weight inputs.
- One-factor expanded beta-prior Cook's replacement-refit helpers that build
  the design internally, then run selected-coefficient or numeric-contrast
  replacement workflows with size-factor or normalization-factor offsets.
- Additive expanded beta-prior Cook's replacement-refit helpers that build the
  design internally, then run selected-coefficient or numeric-contrast
  replacement workflows with size-factor or normalization-factor offsets.
- Formula-driven expanded beta-prior Cook's replacement-refit helpers for the
  supported primitive formula subset, including `offset(numeric)` routing into
  normalization-factor replacement workflows.
- Primitive expanded beta-prior replacement refits preserve optional
  observation weights when rerunning the replacement-count GLM.
- Primitive additive expanded design construction for categorical factors,
  numeric covariates, factor-by-factor, factor-by-numeric, and
  numeric-by-numeric interactions, with one expanded indicator per factor
  level, treatment-style reported columns for non-reference factor levels,
  unchanged numeric columns in both design matrices, primitive pairwise
  interaction products, and matching coefficient groups.
- Primitive formula-to-expanded-design parsing for `1` intercept-only,
  intercept-preserving `+`, `:`, `/`, `*` shorthand, and `0`/`-1` intercept
  removal/restoration terms in formula order, including lower-order-omitted
  pairwise interactions and higher-order interaction/nesting/star-expansion
  terms, primitive `- term` subtraction, additive parenthesized groups with
  signed `+`/`-` terms such as `(condition + batch - batch)`, plain numeric
  identity transforms
  such as `I(dose)` and `I(x=dose)` and signed identities such as `I(-dose)`,
  simple scalar arithmetic transforms such as `I(dose + 1)`, `I(dose / 2)`,
  and `I(1 - dose)`, simple two-covariate
  arithmetic transforms such as `I(dose + time)`, integer numeric power
  transforms such as `I(dose^2)`, raw and default orthogonal polynomial
  transforms `poly(numeric, degree, raw=TRUE)` and `poly(numeric, degree)`,
  including positional or named `x`, `degree`, `raw`, and `simple` arguments,
  factor reference transforms such as `relevel(condition, ref="B")` with
  positional or named `x`/`ref` arguments, including nested
  `relevel(factor(condition), ref="B")`, factor identity transforms such as
  `factor(condition)`, `as.factor(condition)`, `ordered(condition)`, and
  `as.ordered(condition)` for already-supplied factor metadata, simple
  `factor(condition, levels=c(...))` and
  `ordered(condition, levels=c(...))` level reordering with positional or named
  `x`/`levels` arguments, conservative
  `factor(condition, levels=c(...), labels=c(...))` and
  `ordered(condition, levels=c(...), labels=c(...))` relabeling with
  positional or named `x`/`levels`/`labels` arguments, and
  `droplevels(condition)` over already-supplied or simple formula-derived
  factor metadata. `as.factor(...)` and `as.ordered(...)` intentionally reject
  `levels`/`labels` arguments. Nested factor-valued transforms in the supported subset,
  such as `factor(relevel(...))`, `ordered(relevel(...))`,
  `relevel(ordered(...))`, and `droplevels(ordered(...))`, preserve
  formula-local references after relabeling or dropped-level pruning. Quoted
  level and label strings in this supported subset can
  contain spaces, hyphens, parentheses, transform-like text, equals signs,
  escaped quotes, and escaped backslashes. Supported factor-valued transform
  arguments use the same exact-first, unambiguous R-cleaned aliases as
  primitive formula variables, so non-syntactic supplied factor names can be
  addressed consistently inside `factor(...)`, `ordered(...)`,
  `relevel(...)`, and `droplevels(...)`. Supported `relevel(..., ref=...)` references also resolve
  exact levels first and then unambiguous R-cleaned level aliases, storing the
  canonical supplied level after resolution. Supported formulas can also include
  common numeric function transforms `log(numeric)`, `log2(numeric)`,
  `log10(numeric)`, `log1p(numeric)`, `sqrt(numeric)`, and numeric identity
  coercions `as.numeric(numeric)`, `as.double(numeric)`, and
  `as.integer(numeric)`, including named `x` arguments, scalar `log()` `base`
  arguments, `scale(numeric)` with
  R-style named `x` plus named or positional `center`/`scale` arguments, and
  `offset(numeric)` plus single-vector supported transform offsets such as
  `offset(log2(numeric))` or `offset(I(numeric + other_numeric))` into
  per-sample log-offset vectors in the supported formula subset.
  `as.numeric(factor)`, `as.double(factor)`, and `as.integer(factor)`
  materialize 1-based factor codes from declared level order when available,
  otherwise first observed level order. Supported
  numeric transform arguments use the same exact-first, unambiguous R-cleaned
  aliases as primitive formula variables, so expressions such as
  `log2(dose.value)`, `I(dose.value + 1)`, `scale(dose.value)`, and
  `poly(dose.value, ...)` resolve canonical supplied numeric metadata without
  silently accepting alias collisions. Public formula model-frame helpers expose
  the same numeric alias resolution for wrapper and object metadata.
  Derived factor and numeric transform names are ordinary formula variables
  after expansion, so exact or R-cleaned collisions with supplied or previously
  derived variables are rejected instead of silently changing coefficient
  lookup.
  Primitive `.`
  expansion adds all supplied factors
  followed by all supplied numeric covariates as a main effect and inside
  supported `*`, `:`, and `/` shorthand terms, including matching subtraction
  of those expanded terms. Repeated supported formula terms are simplified
  rather than duplicated. Backtick-quoted model-frame column names are accepted
  in supported main effects, interactions, offsets, and numeric transforms.
  Formula variable lookup resolves exact supplied metadata names first, then
  unambiguous R-cleaned aliases, so non-syntactic factor and numeric names can
  be used either quoted or in their cleaned form without silently picking
  collisions.
  R-like unary `+` / `-` sign sequences are accepted in supported term lists
  and additive parenthesized groups, including formula-order `0`/`-1`
  intercept removal and `1`, `-0`, and `- -1` restoration. Signed additive
  groups are reduced inside supported `*`, `:`, `/`, `%in%`, and power
  expressions before operator products are formed. Nested additive
  parenthesized groups are expanded through supported `+`, `*`, `:`, `/`,
  `%in%`, subtraction syntax, and positive-integer formula powers such as
  `(condition + batch + dose)^2`; signed power bases such as
  `(condition + batch + dose - dose)^2` are normalized before interaction
  expansion. Metadata-aware dot terms such as `.`, `condition:.`,
  `condition*.`, `.^2`, and `(.)^2` expand over the supplied factor and
  numeric model-frame columns.
- Additive-factor expanded beta-prior fit-and-Wald-results helpers that build
  the design internally, then run coefficient or numeric-contrast workflows
  with size-factor, normalization-factor, and optional observation-weight
  inputs.
- Formula-driven expanded beta-prior fit-and-Wald-results helpers for the
  supported primitive formula subset, including `offset(numeric)` conversion
  into size-factor or normalization-factor GLM offsets, using the same
  coefficient and numeric-contrast result assembly paths.
- Top-level stored-model-frame Wald and LRT formula fit/result helpers,
  including `results(contrast=...)` and Cook's replacement-refit helpers, apply
  supported full-formula offsets by exponentiating the parsed per-sample log
  offsets and multiplying them into size-factor-derived or caller-supplied
  normalization-factor GLM offsets before the native GLM-mu fit. Reduced LRT
  formulas still reject formula offsets.
- DESeq2-style observation-weight preprocessing helper with row-max
  normalization, weighted design-rank checks, thresholded Cox-Reid sub-design
  checks, and `weights_fail` flags.
- DESeq2-style full-rank design matrix checks for supplied-dispersion Wald/LRT
  pipelines and native dispersion stages.
- DESeq2-shaped fit diagnostics expose stored fitted dispersion trend type
  labels and implemented dispersion-stage aliases, including `dispGeneEst`,
  `dispFit`, `dispersion`, `dispIter`, and `dispOutlier`, plus stable present
  column names, typed data-frame assembly, and TSV export for table exporters.
  CLI parity sidecars can also export full-model beta and beta standard-error
  matrices for original and replacement-refit branches. Full/reduced
  fitted-mean and hat-diagonal matrices remain explicit `DeseqFit` fields
  rather than `mcols`-style vector columns.
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
  `dispGeneEst > 100 * minDisp` viability check and the constant trimmed-mean
  fit over estimates above `10 * minDisp`.
- `normTransform` and mean-fit, parametric, and local numerical-integration VST
  primitives for normalized-count matrices, including dispatch from an
  already-fitted `DispersionTrendFit` and DESeq2-shaped local VST size-factor
  summary helpers. Factor-aware VST dispatch helpers accept ordinary size
  factors or normalization factors directly for the local branch. `DeseqFit`
  retains the implemented fitted trend object and can produce normalized
  counts, `normTransform`, and VST output for its source count matrix, with a
  short `vst` method alias for the transform. It can also apply a caller-supplied
  fitted trend to its full count matrix, which supports the fast-VST split
  between full-data normalization and subset-fitted trend estimation.
  Fast-VST support includes the default `nsub=1000`, explicit-size subset
  helpers, aligned subsetting of counts, normalized counts, normalization
  factors and observation weights, named row preservation, subset trend
  fitting, full-matrix transformation, and diagnostic metadata. For inputs
  without observation weights, the subset selector implements DESeq2 1.52.0's
  ordinary normalized-count row-mean threshold, ordering, and rounded-position
  rule. For weighted inputs, Rust uses its stored weighted `baseMean` for
  eligibility and ordering, whereas DESeq2 `vst()` uses ordinary normalized
  row means. Observation weights are nevertheless subset in the selected row
  order and passed to the weighted dispersion-trend fit. Exact weighted
  subset-index parity is therefore not claimed.

  The automatic Rust helper uses the fast subset when enough rows pass its
  eligibility rule and otherwise fits the trend on all rows. This fallback is
  Rust-specific: DESeq2 1.52.0 `vst()` stops when fewer than `nsub` rows have
  ordinary mean normalized count above 5. Metadata reports the selected path,
  eligible-row count, transform and trend-fit shapes, trend type, and optional
  subset indices. The blind helper uses an intercept-only design with the same
  Rust path selection.
- Pure-Rust locfit-compatible local dispersion trend fitting with DESeq2's
  `dispGeneEst >= 10 * minDisp` local-fit rule, base-mean weights,
  all-near-minimum floor behavior, builder dispatch, offline local-trend
  reference checks, and real-data fitted-value parity at median relative error
  `3.74e-13`, p99 `5.85e-12`, max `1.47e-11`. Against the superseded
  polynomial smoother, the fitted-shape local fixture's maximum absolute error
  changed from `3.29e-04` to `9.62e-05`, and the mixed-threshold fixture
  improved from `3.87e-02` to `1.98e-03`; the small out-of-fit prediction
  fixture changed from `5.76e-05` to `2.56e-04` while remaining inside the
  existing `2e-3` tolerance. Downstream GLM-mu local MAP/Wald/LRT
  fixture metrics were already at machine precision and did not move.
- Deterministic dispersion prior variance estimation, including R-compatible
  MAD scaling, `trigamma((m - p) / 2)` sampling-variance subtraction for
  residual degrees of freedom above 3, saturated designs, low-df histogram/KL
  matching, and the `0.25` floor.
- Initial MAP dispersion fitting with DESeq2's `dispInit` rule,
  `log(dispFit)` prior means, prior-aware line search, grid fallback,
  optional observation weights, `dispMAP`, `dispOutlier`, and final
  dispersion values for the linear-mu and GLM-mu branches.
- Limited native Wald pipeline for the supported linear-mu no-weight and GLM-mu
  optionally weighted, deterministic-prior MAP dispersion subsets with
  parametric, local, or mean dispersion trends.
- Default coefficient-level Wald statistic and standard-Normal p-value.
- Top-level Wald result helpers and CLI commands can report the selected design
  coefficient by index or coefficient name.
- DESeq2-style Wald t p-values with residual, scalar, or per-gene degrees of
  freedom for selected coefficients and primitive numeric contrasts, including
  thresholded t-tail alternatives covered by passing cases from DESeq2's
  `test_nbinomWald.R`.
- Log2-scale beta covariance matrices exposed in `DeseqFit` for implemented GLM
  fits and primitive numeric Wald linear contrasts using `c' beta` and
  `sqrt(c' Sigma c)`.
- Result-row assembly with BH adjustment for precomputed primitive numeric
  Wald contrasts.
- Primitive coefficient-name, positive/negative coefficient-list, and common
  factor-level contrast resolution against design coefficient names, with
  stable result-table names and comparison labels for named contrast specs.
  Coefficient-name and list contrasts resolve exact names first, then
  R-cleaned aliases and intercept aliases; top-level named coefficient and
  list-contrast result helpers share the same alias resolution, including
  formula-built designs with quoted/non-syntactic model-frame column names
  and the corresponding limited Cook's replacement-refit contrast-request
  routes. Supplied-dispersion fixed and limited GLM-mu replacement-refit
  named/list/numeric contrast routes can also infer the same two-level factor
  Cook's low-count heuristic when stored formula metadata proves the numeric
  contrast is exactly the non-reference factor comparison. Builder-created fit
  states retain the formula/model-frame metadata used to build them so wrapper and
  already-fitted-object routes can recover factor levels and references for
  later Wald/LRT `results(contrast=...)` assembly, while explicit
  fitted-object factor-level helpers accept caller-supplied sample levels when
  no stored model frame is available. The
  unified fitted-object dispatcher selects Wald or LRT output from the stored
  test state. R-wrapper already-fitted object routes also infer character
  contrast factor levels from `colData`, preserving declared factor levels and
  treating plain character sample annotations as categorical metadata when
  explicit factor-level metadata is absent.
  `FormulaModelFrame` exposes validation, sample-count, and resolved
  factor-reference accessors so wrapper/object code can inspect the same
  explicit-reference, declared-level, and first-observed fallback rules used by
  formula design construction. `DeseqBuilder::try_model_frame` provides a
  checked wrapper-ingestion path, and explicit model-frame contrast helpers use
  the same validation before storing metadata, so invalid model frames are
  rejected before later formula, contrast, or Cook's/refit routing depends on
  them. Cook's replacement/refit helpers that infer factor-level behavior from
  selected coefficients or numeric contrasts also validate stored model-frame
  metadata before deciding whether factor-aware low-count handling applies.
  Explicit model-frame references are resolved exact first and then by
  unambiguous R-cleaned aliases against declared levels and observed sample
  levels, returning the canonical raw level for downstream character contrasts
  and Cook's low-count conditions. Formula design construction
  also returns a formula-local model frame beside offsets, so supported derived factor labels
  such as `relevel(...)`, `factor(..., levels=...)`,
  `factor(..., levels=..., labels=...)`, and `droplevels(...)` can
  drive character contrast metadata and all-zero handling through top-level
  formula result, fixed-dispersion model-frame result, and supported Cook's
  replacement routes. Supported formula-built Wald/LRT replacement-refit
  routes also retain formula-local numeric metadata for selected coefficients,
  coefficient-list contrasts, and numeric-vector contrasts created by supported
  numeric transforms such as `as.double(factor)`, with the same exact-first
  R-cleaned coefficient alias lookup where names are involved. Fitted objects
  from supported formula and model-frame routes retain this metadata for later
  object-style result requests.
  Non-reference factor-level
  comparisons can infer a shared reference from coefficient names such as
  `B_vs_A` and `C_vs_A`. Factor-level coefficient candidates include R-style
  whole-name and component-wise `make.names` cleanup for non-syntactic and
  reserved factor or level names, including stored formula-transform labels
  such as `factor(..., levels=...)`, `relevel(...)`, and `droplevels(...)`,
  with ambiguity errors for treatment-style and expanded/no-intercept alias
  collisions.
  Coefficient-list contrast weights follow DESeq2 `listValues` sign
  validation, and list comparison labels follow DESeq2's two-sided and
  one-sided naming shape. Interaction coefficient aliases can be cleaned
  component-wise for coefficient-name and list contrasts while preserving the
  `:` separator. Resolved Rust character contrast metadata keeps the requested
  factor name together with numerator and denominator levels so downstream
  object and wrapper routes can keep `contrastAllZero` bookkeeping
  self-contained.
- Native linear-mu and GLM-mu Wald contrast entry points reuse the implemented
  MAP dispersion paths, then run numeric, named/list, or caller-supplied
  factor-level contrast result assembly. Compatibility-named parametric-only
  Wald helpers expose the same contrast results while pinning the
  dispersion trend to the parametric branch.
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
- Primitive contrast result builders validate finite effect estimates and
  standard errors before assembling Wald or LRT result tables.
- Typed DESeq2-shaped result-table view for implemented primitive rows, with
  row names, numeric/logical columns, and column metadata for wrapper/file
  output, plus regular and `results(tidy = TRUE)`-style TSV result export,
  `mcols(res)`-style result-column metadata export, and table-level result
  metadata export.
- Supplied-dispersion fixed-dispersion Wald pipeline for one coefficient and
  primitive numeric contrasts.
- Supplied-dispersion fixed-dispersion LRT pipeline for full vs reduced design
  matrices, with full and reduced log-likelihood, beta convergence, and
  iteration diagnostics in `DeseqFit`, plus primitive numeric, named/list, and
  caller-supplied factor-level effect-size contrast result tables.
- Limited native-dispersion LRT pipelines for the implemented linear-mu and
  GLM-mu MAP dispersion branches, with the full design used for dispersion
  estimation before full-vs-reduced testing, and top-level result helpers that
  can report the full-design coefficient by index or coefficient name.
- Native linear-mu and GLM-mu LRT result helpers can report primitive numeric
  or named/list full-model effect-size contrasts, plus caller-supplied
  factor-level contrasts with sample-level all-zero cleanup, while preserving
  the full-vs-reduced LRT statistic and p-values. Compatibility-named
  parametric-only LRT helpers expose the same contrast results. Fit-only
  top-level LRT helpers mirror the default, named-coefficient,
  numeric-contrast, named-contrast, and factor-level contrast result routes.
- Numeric/expanded DESeq2-style `contrastAllZero` handling for LRT contrast
  result tables zeroes only the reported `log2FoldChange`; the LRT statistic,
  p-value, and adjusted p-value remain the model-comparison outputs. The same
  LFC-only cleanup is available for character/factor-level LRT contrasts when
  the caller supplies sample levels.
- R-wrapper primitive and already-fitted-object `results()` routes support
  DESeq2-style character triplet, coefficient-list/listValues, and numeric
  `contrast = ...` forms for Wald and LRT fits. Character triplets use supplied
  `reference` metadata, then explicit wrapper `factorReferences`, then the
  first declared factor level before falling back to observed sample order from
  primitive `factorLevels` / factor-valued `colData`; when counts are
  available, character contrast all-zero handling can also use the matching
  `colData` factor values as per-sample levels.
  Already-fitted object routes can read nested count assays from plain lists or
  bracket-addressable assay containers.
  Named `colData` rows are aligned to named count columns before wrapper
  contrast metadata uses them. Wrapper `factorLevels` and `colData` factor
  lookups resolve exact metadata names first, then unambiguous R-cleaned aliases
  for non-syntactic field names. Wrapper character contrast all-zero handling
  resolves requested numerator, denominator, and explicit reference levels
  exact first, then through unambiguous R-cleaned aliases against observed
  sample levels.
  Declared factor levels may define an unused treatment reference for both
  model-matrix construction and character contrast metadata, matching R factor
  model-matrix behavior. Formula model-frame character contrasts, including
  supplied-dispersion model-frame routes and stored model-frame formula helpers,
  resolve exact factor and level names first, then unambiguous R-cleaned
  aliases for factor names, numerator/denominator levels, and explicit or
  stored reference levels. The public model-frame factor-reference metadata
  helper exposes the same exact-first cleaned-alias resolution for
  wrapper/object code.
  Explicit-reference factor contrasts can resolve either treatment-style
  `factor_level_vs_reference` columns or expanded/no-intercept `factorLevel`
  columns; `contrast` takes precedence over `name`, matching DESeq2.
- All-zero row expansion for the supplied-dispersion Wald/LRT pipelines and
  limited native Wald/LRT paths, using missing numeric outputs for skipped
  all-zero rows.
- Cook's distance matrix for the supplied-dispersion and limited native Wald
  pipelines.
- Cook's diagnostic TSV exports for the distance assay, row-level `maxCooks`
  and robust dispersion, and sample-level `samplesForCooks` eligibility.
- `maxCooks` over samples in model-matrix cells with at least three replicates.
- Cook's cutoff p-value masking with BH recomputation for result rows.
- Explicit primitive helper for DESeq2's two-group low-count Cook's heuristic:
  rows above cutoff are spared when at least three counts are larger than the
  count in the sample with maximum Cook's distance. Supplied-dispersion fixed
  Wald/LRT and limited GLM-mu Wald/LRT factor-level result routes apply
  this automatically when the caller-supplied sample levels contain exactly
  the requested numerator and denominator levels. The same condition is also used
  after limited GLM-mu factor-level replacement refits. Selected-coefficient
  Wald/LRT routes, including the limited replacement-refit paths, also apply it
  when stored formula model-frame metadata proves the selected coefficient or
  numeric contrast is the non-reference comparison, or its reverse, for a
  single two-level factor; stored references are resolved against observed
  levels exact first, then by unambiguous R-cleaned aliases. Direct
  factor-level character all-zero checks use the same exact-first unambiguous
  R-cleaned alias rule when matching requested numerator and denominator levels
  against observed sample labels.
- Primitive Cook's outlier count replacement transform from `replaceOutliers`:
  trimmed normalized means, size-factor or normalization-factor rescaling,
  integer truncation, replaceable-sample masks, and per-gene `replace` flags.
- Cook's replacement-refit planning metadata: replacement-count base metadata,
  `nrefit`, `refitReplace`, `newAllZero`, and DESeq2-style post-refit
  `maxCooks` masking over nonreplaceable samples.
- Cook's replacement-refit plans expose a compact scalar metadata view for
  `nRefit`, refit-row counts, new-all-zero rows, outlier/replaced cell counts,
  replaceable sample counts, and the refit-branch decision, plus TSV export for
  that key/value metadata.
- Cook's replacement/refit plans can export replacement-count assays,
  candidate replacement counts, outlier-cell logical masks, and row-level
  replacement/refit metadata for wrapper and parity-table consumers.
- Limited Cook's replacement-refit path for supplied-dispersion fixed Wald/LRT
  coefficient, primitive contrast, named/list contrast, and caller-supplied
  factor-level contrast rows, plus the implemented GLM-mu native Wald and LRT
  branches, merging only `refitReplace` rows and preserving original size
  factors.
- Limited Cook's replacement-refit path for primitive expanded beta-prior Wald
  selected-coefficient and numeric-contrast workflows using size-factor or
  normalization-factor offsets.
- Limited Cook's replacement-refit path for GLM-mu native Wald contrasts,
  including primitive numeric, named/list, and caller-supplied factor-level
  contrast routes.
- Test-selected top-level Cook's replacement-refit helpers can return a typed
  Wald-or-LRT enum, so `test=Lrt` can reuse a builder-stored reduced design for
  default, named-coefficient, numeric-contrast, named-contrast, and factor-level
  contrast replacement-refit workflows. The enum exposes shared accessors for
  selected test type, original fit/results, replacement refit plan, optional
  refit fit/results, and final merged results.
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
and the limited GLM-mu native Wald replacement-refit path estimates gene-wise
dispersions on replacement counts, reuses the original dispersion function and
dispersion prior variance for MAP shrinkage, and reruns Wald testing with
original size factors.
The limited GLM-mu native LRT replacement-refit path uses the same replacement
bookkeeping and replacement-count size-factor preservation before merging refit
rows. Full Bioconductor `refitWithoutOutliers` behavior is not implemented.
The limited native Wald pipeline follows the high-level DESeq2
stage order in `DESeq`: size factors, dispersion estimation, then
`nbinomWaldTest`, but only for the supported linear-mu no-weight
and GLM-mu optionally weighted dispersion branches. Cook's cutoff masking
follows the default `results()` path in
`R/results.R`, where `maxCooks > qf(.99, p, m - p)` makes the p-value missing
before p-value adjustment. Independent filtering follows `pvalueAdjustment`
and `filtered_p` in `R/results.R`: `baseMean` is the default filter statistic,
candidate cutoffs are filter quantiles, and BH adjustment is recomputed after
filtering. An optional ignored DESeq2 source checkout may be placed at
`external/DESeq2` for inspection; it is not vendored into this package.

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
reduced-model convergence, iteration, fitted-mean, and hat-diagonal diagnostics
in `DeseqFit`. The `mcols`-style diagnostic view intentionally stays
vector-shaped: matrix-valued full/reduced fitted means and hat diagonals are
available on the fit object and are not exported as row metadata columns.

The gene-wise dispersion foundation follows these early
`estimateDispersionsGeneEst` definitions in `R/core.R`: `roughDispEstimate`,
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
check is applied first, then the constant trend is the trimmed mean of estimates
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
fixed-dispersion references. Some legacy fixed-dispersion references use
`DESeq2:::fitNbinomGLMs` with caller-supplied dispersions, default `1e-6` beta
ridge, `useQR=FALSE`, and `useOptim=FALSE`; Rust keeps
`IrlsSolver::NormalEquations` for those checks while the default builder path
uses the DESeq2-style augmented QR solver. That reference set includes
fixed-dispersion beta, SE, Wald/LRT statistics, log likelihoods,
fitted means, hat diagonals, and Cook's distances. It also writes weighted base
metadata using `getBaseMeansAndVariances` and `getAndCheckWeights`, weighted
fixed-dispersion Wald/LRT references with `weightsFail` rows expanded as
missing, including reduced-model fitted means and hat diagonals for weighted
LRT, unweighted GLM-mu Cox-Reid mean-trend MAP/Wald/LRT and mean-trend
MAP/Wald/LRT intermediates, plus weighted GLM-mu Cox-Reid mean-trend MAP and
mean-trend and local-trend dispersion/MAP/Wald/LRT intermediate references for
the supported native branch. The matched GLM-mu Wald/LRT branches also write compact
DESeq2-shaped result-table fixtures and assert result-row beta, SE, statistic,
p-value, adjusted p-value, dispersion, and convergence parity.
The unweighted GLM-mu local-trend fixture covers MAP, Wald, LRT, and compact
Wald/LRT result rows for the tiny-data edge case where only one row is usable
for local fitting; Rust follows the corresponding constant local fit shape
instead of failing the trend evaluation.
The weighted GLM-mu local-trend fixture covers the same MAP/Wald/LRT results
with DESeq2's `weightsFail` expansion semantics.
The GLM-mu Cox-Reid local-trend fixtures cover the unweighted and
weighted MAP/Wald/LRT intermediates plus compact result rows, including
Wald fitted means and hat diagonals, LRT full/reduced likelihoods, deviances,
convergence, and `weightsFail` row expansion.

The fixed-dispersion IRLS path includes a bounded limited-memory BFGS-style
pure-Rust optim fallback for routed rows. Its primary candidate follows
DESeq2's objective-only `optim(..., method="L-BFGS-B")` shape with
independently implemented R-compatible negative-binomial arithmetic, R-style
finite differences, fused predictor accumulation, and DESeq2's raw natural-log
QR backup values passed into the log2 optimizer. For fallback rows, DESeq2 keeps the pre-optimizer IRLS hat
diagonals for downstream Cook's decisions, so Rust preserves those hats too
while using the selected refit beta, SE/covariance, fitted means, row log
likelihoods, and convergence diagnostics.

Final/refit fits also compute an analytic endpoint. A stability check that does
not use reference data selects it only when the analytic solver succeeds, the
compatible objective exceeds the analytic objective by more than
`10 * factr * f64::EPSILON * max(abs(compatible_objective), abs(analytic_objective), 1)`,
and the endpoints differ by more than ten finite-difference steps. Gene-wise dispersion uses the
analytic route directly; the callback-compatible route is limited to
final/refit fitting.

The fallback uses `rcompat-lbfgsb` 0.2.1. Its precision was checked against the
512-case synthetic stress bundle generated by R 4.6.1 (`stats::optim`) with
OpenBLAS 0.3.32, one BLAS thread, and 17-digit round-trip serialization. At the
scan's practical thresholds (maximum parameter error `5e-3`, absolute objective
error `1e-5` or relative objective error `1e-8`), 0.1.6 matched both endpoint
and objective in 493/512 cases, matched the objective in 507/512 cases, and
matched R's function/gradient counts in 311/512. Version 0.2.1 matched all 512
endpoints, objective values, and counts.
With all endpoint and objective tolerances set to zero, the result improved
from 0/512 exact cases on 0.1.6 to 512/512 on 0.2.1. These results apply to the
recorded R version, platform, and objective fixtures; they do not establish
bit-exact behavior for arbitrary platforms or callbacks.

In a dependency-only 65,580-gene kidney replay, the exact optimizer gain did
not materially improve end-to-end Wald parity. From 0.1.6 to 0.2.1, the LFC
median/p99/maximum changed from
`3.11e-14` / `3.79e-12` / `7.70e-4` to
`3.15e-14` / `3.79e-12` / `7.70e-4`; the Wald-statistic values changed from
`5.27e-12` / `3.70e-11` / `1.25e-3` to
`5.33e-12` / `3.70e-11` / `1.25e-3`. The `lfcSE` maximum changed from
`1.95e-6` to `5.10e-6`, and the BH-adjusted-p maximum from `4.50e-5` to
`7.29e-5`. Fallback rows are a small part of the workflow, so other numerical
stages dominate the remaining tails. The complete table is stored in
[`data/lbfgsb_real_data_precision.tsv`](data/lbfgsb_real_data_precision.tsv).
For an analytic-gradient-only variant, median/p99 errors changed from
`3.15e-14` / `3.79e-12` to `2.46e-14` / `3.03e-12` for LFC, from
`1.77e-12` / `1.88e-10` to `1.38e-12` / `1.50e-10` for SE, and from
`5.33e-12` / `3.70e-11` to `4.18e-12` / `2.96e-11` for the Wald statistic
(about 22%/20% lower). Maximum errors were unchanged. Version 0.2.5 uses
the callback-compatible solution plus the stability check described above.

On the fixed 100 rows with the largest v0.2.4 errors from a 69,045-gene heart
contrast, the measured v0.2.5 median absolute error is `6.064e-10`, mean error
is `1.261e-4`, and maximum error is `1.526e-3`. The corresponding v0.2.4
measurements are `1.464e-4`, `3.794e-4`, and `3.094e-3`: using unrounded
values, v0.2.5 is 241,402x lower at the median, 3.009x lower at the mean, and
50.68% lower at the maximum. Of the 100 rows, 89 improved and 78 improved by
at least 10x. The scorer requires all 100 genes to be present in the supplied
diagnostics and does not infer errors for missing rows.

The large median change comes from compensated accumulation of per-gene log
counts during size-factor estimation. It retains low-order information in a
portable binary64 implementation, closely tracks R's extended-precision row
means, and keeps near-boundary dispersion fits on the same route as DESeq2.

Only 26/65,580 kidney rows (`0.040%`) and 305/535,178 fitted rows across eight
real-data contrasts (`0.057%`) used L-BFGS-B. On eight non-replaced kidney tail
rows, median dispersion drift of `4.77e-08` corresponded to median optimizer
beta-target drift of `9.82e-05`, illustrating why trajectory-compatible
arithmetic and a stability check matter on flat objectives. Three
whole-workflow runs showed no clear
speed change: the 0.1.6 and 0.2.1 ranges overlapped (`25.87–28.09 s` versus
`25.76–28.51 s`). A direct finite-difference/analytic-gradient comparison also
overlapped (`28.46–28.66 s` versus `28.18–28.77 s`), with a descriptive 0.8%
median reduction. See
[`data/lbfgsb_real_data_route_summary.tsv`](data/lbfgsb_real_data_route_summary.tsv).
The v0.2.5 100-row validation run measured 128.625 s and 2,947,976 KiB peak
RSS, while the v0.2.4 run measured 120.719 s and 2,945,168 KiB. One run per
version measured absolute differences of +7.906 s and +2,808 KiB (+6.55% and
+0.095%); it does not establish a whole-workflow runtime or memory change.

The publication-style Wald validation set covers high-error, Cook's
replacement, and missingness-sensitive full-blocked contrasts across three
tissues. Its 278,257 result rows have zero missing-row or finite/NA-pattern
mismatches across `baseMean`,
`log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. The largest observed
absolute differences in those checks were `5.700e-4` for p-value and
`8.748e-4` for adjusted p-value; remaining differences are finite numeric
drift rather than row swaps or missingness cascades.

`scripts/analyze_parity_mismatch_sources.py` joins real-data diagnostics with
offline optimizer fixture summaries to separate initial beta differences from
later changes to statistics and adjusted p-values. Across the selected
high-error contrasts, optimizer-refit rows had the largest numeric errors.
Smaller errors arose from MAP dispersion line-search sensitivity, and the
largest adjusted-p-value errors reflected Benjamini-Hochberg propagation from
upstream p-values.

## Missing

- Mature wrapper-facing interface around the Rust core.
- Complete formula parsing in Rust, including arbitrary R expressions, splines,
  and full R-compatible formula semantics beyond the supported native formula
  subset.
- Full DESeq2 dispersion estimation, including broader weighted dispersion
  edge-case parity, broader synthetic `locfit` edge-case parity, glmGamPoi trend
  type, and end-to-end dispersion parity across more designs.
- High-level propagation of observation-weight `weights_fail` flags through
  complete DESeq2-like builder and unimplemented wrapper workflows.
- Direct weighted low-level `fitNbinomGLMs` parity for rows that DESeq2 marks
  `weightsFail` but still returns ridge-stabilized coefficients when that
  internal function is called directly.
- Broader DESeq2 parity fixtures for unstable or non-converged rows routed
  through the bounded optim fallback.
- DESeq2-style `results(contrast=...)` semantics are implemented for primitive
  Rust, CLI, and R wrapper already-fitted result routes: character triplets,
  coefficient-list/listValues, numeric contrasts, factor-level reference
  inference from fitted metadata where available, explicit fitted-object
  factor-level routes with caller-supplied sample levels, contrast-specific
  all-zero behavior for character and two-sided numeric/list contrasts,
  `contrast` precedence over `name`, and fitted-object dispatch from stored
  Wald or LRT state. High-level Bioconductor object integration and additional
  contrast-aware Cook's/refit edge cases remain unsupported.
- Full Cook's outlier replacement behavior for high-level Bioconductor assay
  attachment, high-level wrapper-object metadata, and all remaining DESeq2 edge
  cases.
- General Wald and LRT tests with native dispersion estimation beyond the
  supported limited linear-mu/GLM-mu branches and without generated DESeq2
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
- High-level R-style contrast handling for arbitrary unfitted `DESeqDataSet`
  objects beyond the implemented primitive and already-fitted-object
  `results()` routes.
- Full high-level VST object workflow, exact DESeq2 `splinefun` behavior for
  local VST, full high-level rlog object workflow, and lfcShrink-compatible
  hooks. Mean-fit, parametric, local, fast-subset selection/subsetting,
  automatic GLM-mu VST, blind automatic VST, and fit-state VST helpers are
  available for normalized-count matrices and implemented fitted trends.
  Low-level rlog sample-effect ridge fitting is available when callers supply
  dispersions, with helpers to estimate the sample-effect prior variance from
  normalized counts, `baseMean`, and `dispFit`, or to estimate that prior and
  fit in one call for size-factor and normalization-factor inputs. `DeseqFit`
  can run the same rlog path after fitted trend dispersions and final
  dispersions are available in stored fit state. Builder-level `rlog_glm_mu`
  and `blind_rlog_glm_mu` helpers compose the implemented GLM-mu MAP
  dispersion path with fit-state rlog, with `*_with_fit` variants that retain
  the dispersion fit state for diagnostics. The native CLI exposes this path
  through `rlog` for blind and design-aware GLM-mu MAP dispersion workflows.
- Top-level builder Wald helpers route through the implemented GLM-mu path.
  Top-level LRT helpers can store a reduced design on the builder and route
  default, named-coefficient, numeric-contrast, coefficient-name contrast, and
  factor-level contrast result requests, including the typed top-level
  replacement-refit helpers, through the implemented GLM-mu LRT path.

## Known Deviations

The Rust core accepts primitive matrices and explicit options. Formula
semantics and model-matrix generation require wrapper support.

Gene/sample normalization factors are supported for normalized-count metadata,
supplied-dispersion fixed Wald/LRT pipelines, and the supported native linear-mu
dispersion path. The native subset follows DESeq2's
`linearModelMuNormalized` and `momentsDispEstimate` branches: fitted raw means
use `linearModelMu(counts(normalized=TRUE), X) * getSizeOrNormFactors(object)`,
and moments starts use `mean(1 / colMeans(normalizationFactors))` on the
non-all-zero row subset.

The size-factor implementation follows DESeq2 1.52.0
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

The implemented gene-wise dispersion estimators have the Cox-Reid objective and
Armijo line-search shape, including DESeq2's weighted Cox-Reid
`weightThreshold` subset rule. The linear-mu branch supports trend fitting,
MAP, and the limited native Wald path. The GLM-mu branch performs
mean/dispersion alternation using fixed-dispersion IRLS, can consume
preprocessed builder observation weights, and supports parametric, local, and
mean trend fitting, MAP shrinkage, and the limited native Wald path.
DESeq2 switches away from linear-mu fitting when observation weights are
present. `rsdeseq2` implements the log-prior objective, first and second
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

The 48-row focused fixture generated by
`scripts/generate_se_covariance_hard_fixtures.R` isolates a no-optimizer
standard-error tail. When DESeq2's final MAP dispersions are supplied to the
fixed-dispersion GLM, `betaSE` matches at machine precision (`1.44e-15` maximum
absolute difference). A compact regression test injects selected fixture
dispersions and checks beta, `betaSE`, fitted means, and hat diagonals at tight
tolerances. The remaining tail is therefore MAP dispersion line-search
sensitivity. Initial
MAP derivatives already match tightly (`2.96e-14` max absolute difference),
accepted-step counts match on all focused rows. Source-order negative-binomial
likelihood accumulation, a Stirling log-gamma branch for large positive
arguments, and R's positive-argument digamma calculation reduce the median
focused `dispMAP` difference to `4.72e-16` and the maximum to `2.31e-07`. The
MAP posterior uses plain source-order addition for DESeq2-like Armijo boundary
behavior; this drops the focused line-search iteration-count mismatches from
13/48 to 4/48 while keeping accepted-step counts matched for all focused rows.
The remaining differences are in final near-boundary Armijo history and final
`log(alpha)`, not covariance code.

Rust `controlGenes` can be provided as zero-based indices or a logical mask.
DESeq2 numeric control genes are one-based because they are R indices.

## Parity Checks

Each test defines explicit absolute and relative tolerances for its numeric
fields; there is no repository-wide tolerance. The normalization benchmark
reported zero finite/NA mismatches and a maximum relative error of `9.887e-15`
across 612.7 million normalized cells.

GLM and dispersion comparisons include these intermediate fields:

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
