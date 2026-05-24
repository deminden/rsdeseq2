# Algorithms

## Model Reference

DESeq2 models counts with a negative-binomial GLM:

```text
K_ij ~ NB(mu_ij, alpha_i)
mu_ij = s_j * q_ij
log(q_ij) = x_j * beta_i
```

With gene/sample normalization factors, DESeq2 replaces `s_j` with
`NF_ij`:

```text
mu_ij = NF_ij * q_ij
```

The Rust core will use natural-log scale internally for GLM fitting. DESeq2
reports log2 fold changes, so future result-building code must convert
coefficient estimates explicitly.

## Size Factors

The initial implementation supports two size-factor methods.

`ratio` computes a geometric mean for each gene across all samples. Genes with
any zero count receive a zero geometric mean and are skipped. For each sample,
DESeq2 applies the location function to log count/geometric-mean ratios and then
exponentiates. With the default median, even-sized contributing sets use the
geometric midpoint of the two middle ratios on the original scale.

`poscounts` follows DESeq2's implementation in
`estimateSizeFactorsForMatrix`: zero log counts are replaced with zero before
row means are computed, and all-zero genes are skipped. This is equivalent to
using the nth root of the product of positive counts, where n is the total
number of samples. The sample-wise location is again computed on the log-ratio
scale and exponentiated.

When supplied geometric means are used, size factors are stabilized to have
geometric mean 1, matching DESeq2's frozen normalization behavior.

## Normalized Counts

When only size factors are present:

```text
normalized_count_ij = count_ij / size_factor_j
```

When gene/sample normalization factors are supplied, they preempt size factors:

```text
normalized_count_ij = count_ij / normalization_factor_ij
```

This follows DESeq2 `R/methods.R`: `counts(dds, normalized=TRUE)` returns
`counts / normalizationFactors(dds)` when that assay is present. The same
normalization-factor matrix is used as the fixed-dispersion GLM offset source,
following `getSizeOrNormFactors` in `R/core.R` and `fitNbinomGLMs`.
`io::write_count_matrix_tsv()` and `io::write_normalized_counts_tsv()` export
raw or normalized genes x samples matrices with a first `gene` column and
sample-name columns, matching the practical shape of DESeq2 count matrix
exports. `io::write_normalization_factors_tsv()` exports the same genes x
samples shape for finite positive normalization-factor matrices.

## Base Mean

For each gene:

```text
baseMean_i = mean_j(normalized_count_ij)
```

When DESeq2 observation weights are present for early row metadata, normalized
counts are multiplied by weights before the row summary:

```text
weighted_normalized_ij = normalized_count_ij * weight_ij
baseMean_i = mean_j(weighted_normalized_ij)
```

This mirrors `getBaseMeansAndVariances` in DESeq2. It is not a weighted-mean
estimator; the weights are applied to the normalized-count matrix first and an
ordinary row mean is then computed.

This is tested independently so future dispersion and result-table code can
compare against DESeq2 intermediate values.

## Design Matrix Rank

DESeq2 checks GLM model matrices with `qr(modelMatrix)$rank < ncol(modelMatrix)`
and stops when the design is not full rank. The Rust core mirrors that guard
for supplied-dispersion Wald/LRT pipelines and the current native linear-mu
dispersion path. `DesignMatrix` exposes deterministic rank helpers using
partial-pivot elimination with a fixed tolerance, and rank errors report
zero-only coefficient columns when they can be identified.

## Base Variance And All-Zero Rows

DESeq2's early row metadata includes:

```text
baseMean = rowMeans(normalized_counts)
baseVar  = rowVars(normalized_counts)
allZero  = rowSums(raw_counts) == 0
```

For weighted row metadata:

```text
baseVar = rowVars(normalized_counts * weights)
```

`rsdeseq2` mirrors this in `DeseqFit`. `baseVar` uses sample variance. For one
sample, variance is undefined in R, so the Rust value is `NaN`.
`io::write_base_metadata_tsv()` exports the three implemented row metadata
columns with DESeq2 names (`baseMean`, `baseVar`, `allZero`) and R-style
`NA`/`TRUE`/`FALSE` values for parity logs and wrapper-facing tables.

## Multiple Testing

Benjamini-Hochberg adjustment ranks non-missing p-values, applies
`p * m / rank`, walks backward with a cumulative minimum, and caps values at
1. Missing p-values remain missing.

## Roadmap

Next numerical stages should be implemented in this order:

1. Negative-binomial log-likelihood.
2. Fixed-dispersion GLM beta fitting with IRLS.
3. Wald/LRT statistics and inspectable GLM log-likelihood plus beta
   convergence/iteration metadata for fixed dispersions.
4. Gene-wise dispersion estimation.
5. Parametric dispersion trend, prior variance, and MAP shrinkage.
6. Limited native Wald path for the current linear-mu and GLM-mu MAP subsets.
7. Full Wald and LRT parity.

## Negative-Binomial Likelihood

DESeq2 evaluates likelihoods with R's `dnbinom` using mean `mu` and size
`1 / dispersion`:

```text
log p(y | mu, alpha) =
  lgamma(y + 1/alpha) - lgamma(1/alpha) - lgamma(y + 1)
  + (1/alpha) * log((1/alpha) / (1/alpha + mu))
  + y * log(mu / (1/alpha + mu))
```

The Rust helpers in `glm::nb` implement this parameterization and row-sum
log-likelihoods matching DESeq2's `nbinomLogLike`.

## Variance-Stabilizing Transform

DESeq2's `normTransform` is available as the direct
`log2(normalized_count + 1)` transform for already-normalized count matrices.
Like VST and rlog, it is a secondary visualization transform and is not used by
the differential-expression GLM.

The implemented transformation subset includes the two closed-form DESeq2 VST
branches plus the local numerical-integration branch for already-normalized
counts.

For `fitType="mean"`:

```text
vst(q, alpha) = (2 * asinh(sqrt(alpha * q)) - log(alpha) - log(4)) / log(2)
```

where `q` is a normalized count and `alpha` is the mean dispersion.

For `fitType="parametric"` with fitted trend coefficients `asymptDisp` and
`extraPois`:

```text
vst(q) =
  log(
    (1 + extraPois + 2 * asymptDisp * q
       + 2 * sqrt(asymptDisp * q * (1 + extraPois + asymptDisp * q)))
    / (4 * asymptDisp)
  ) / log(2)
```

For `fitType="local"`, Rust evaluates the fitted local dispersion trend on the
same `sinh(seq(asinh(0), asinh(max(q)), length.out=1000))[-1]` grid shape used
by DESeq2, integrates `1 / sqrt(dispersion(q) * q^2 + xim * q)` with a
trapezoid rule, and rescales the integral using type-7 95% and 99.9% quantiles
of row means so high counts follow `log2(q)`. The current interpolation is
deterministic linear interpolation over the integrated grid. The `xim` term can
be computed from size factors with `mean(1 / sizeFactors)`, or from
normalization factors by DESeq2's local-VST approximation
`sf = exp(colMeans(log(normalizationFactors)))` followed by `mean(1 / sf)`.
Convenience dispatch helpers accept ordinary size factors or a
normalization-factor matrix directly and compute this local variance term
internally; mean and parametric branches ignore the term after factor
validation.

For large counts all three transforms converge toward `log2(q)`. Remaining VST
parity work includes frozen dispersion-function reuse, exact DESeq2 `splinefun`
behavior for the local branch, broader object metadata, and remaining
high-level wrapper semantics. A lower-level
`vst_with_dispersion_trend` and its factor-aware variants mirror DESeq2's
`getVarianceStabilizedData` branch selection once the caller already has a
fitted `DispersionTrendFit`. `DeseqFit` also retains the implemented fitted
trend object and exposes fit-level normalized-count, `normTransform`, and VST
helpers, including a short `vst` alias, so callers can transform the same count
matrix without reconstructing the dispatch inputs by hand. A fit-level
`variance_stabilizing_transform_with_trend` helper also accepts an external
fitted trend, which is the reusable primitive for applying a fast-subset trend
to the full normalized count matrix.

The fast `vst()` wrapper in DESeq2 estimates the dispersion trend on a
deterministic subset before applying the fitted transform to all rows. The
implemented Rust helper for this subset keeps rows with mean normalized count
above 5, orders them by base mean, and selects
`round(seq(1, n, length.out=nsub))` positions using R-style half-to-even
rounding. The exported `DEFAULT_FAST_VST_NSUB` constant is `1000`, matching the
DESeq2 default, while explicit-size helpers remain available for small
fixtures and diagnostics. `fast_vst_eligible_count` and the fit-level
eligibility helper expose how many rows pass the same `baseMean > 5` and finite
input checks before a default-size subset is requested. Companion helpers
return selected normalized-count rows and aligned gene/sample matrix rows, such
as normalization factors or observation weights, in that same deterministic
order for the subset trend fit. A row-aligned
`FastVstSubset` bundle combines the raw count subset, normalized-count subset,
optional normalization factors, optional observation weights, and original row
indices so downstream trend fitting can use one shared subset rule. The bundle
also exposes a compact metadata view with subset shape, original row indices,
and normalization-factor/observation-weight presence flags. `DeseqFit`
exposes the same bundle from its stored `baseMean`, normalization factors, and
preprocessed observation weights, which is the intended entry point for later
automatic fast-VST trend re-estimation.
`DeseqBuilder::fit_fast_vst_dispersion_trend_glm_mu` now fits the selected
GLM-mu dispersion trend on that deterministic subset while preserving full-data
size factors and subset normalization factors. The paired `fast_vst_glm_mu`
helper applies the subset-fitted trend back to the full normalized count matrix,
including the normalization-factor dispatch when factors are present, and
returns a named output containing the transformed matrix, subset fit, and row
bundle for diagnostics. `FastVstGlmMuOutput::metadata()` summarizes the
transformed matrix shape, deterministic subset shape, original subset row
indices, subset trend-fit shape, and stable trend-fit type label
(`parametric`, `local`, or `mean`) for direct fast-VST benchmark logs.
The public fast-VST builder validates `nsub > 0` before other branch checks so
invalid subset requests report the same error shape as the lower-level subset
helpers and automatic VST.
Observation-weighted fast-VST trend fitting remains
explicit future work because preprocessed weight rows must not be treated as
raw input weights a second time.
`DeseqBuilder::vst_glm_mu_auto` now performs the DESeq2-shaped automatic
choice: use the deterministic fast subset when enough rows are eligible and
observation weights are absent, or fit the selected Rust trend on all rows when
the default subset is too large for the dataset or weights are present. The
returned `VstGlmMuOutput` records whether the trend came from the fast subset
or full data, the requested `nsub`, and the eligible-row count;
`VstTrendSource` also exposes accessor helpers for these fields and for the
fast-subset decision. Full-data trend metadata records whether that path was
chosen because too few rows were eligible or because observation weights were
present. The source and full-data reason enums expose stable labels:
`fastSubset`, `fullData`, `insufficientEligibleRows`, and
`observationWeightsPresent`. `VstGlmMuOutput::metadata()` packages those
labels with `nsub`, eligible-row count, transform dimensions, trend-fit row
and sample counts, stable trend-fit type label, optional fast-subset row count,
and optional original fast-subset row indices for wrappers and benchmark logs.
`blind_vst_glm_mu_auto` uses the same automatic choice with a named
intercept-only design, matching the implemented part of DESeq2's `blind=TRUE`
workflow without invoking any external runtime.
`CountMatrix::select_rows` provides the matching raw-count subset while
preserving gene and sample names, which mirrors the object subset passed to the
fast-VST dispersion fit. Full fast-VST parity still requires object metadata
and exact local interpolation semantics.

## Gene-Wise Dispersion Objective

The current gene-wise dispersion foundation follows DESeq2's fixed-mean
dispersion optimizer shape. For one gene, the alpha-dependent likelihood kernel
used for scoring is:

```text
sum_j [
  lgamma(K_j + 1 / alpha)
  - lgamma(1 / alpha)
  - K_j * log(mu_j + 1 / alpha)
  - (1 / alpha) * log(1 + mu_j * alpha)
]
```

Terms independent of `alpha` are omitted because they do not affect the
optimizer. When observation weights are supplied to the low-level dispersion
objective, they multiply the per-sample likelihood terms. With Cox-Reid
correction, DESeq2 uses those weights to choose the `weightThreshold`
sample subset, then computes the determinant on the unweighted NB working
variance diagonal:

```text
-0.5 * log(det(X' W X))
selected = { j | obs_weight_j > weightThreshold }
W_jj = (1 / mu_j + alpha)^-1 for j in selected
```

The Rust objective can also add DESeq2's normal prior kernel on the
log-dispersion scale:

```text
-0.5 * (log(alpha) - log(alpha_prior_mean))^2 / prior_var
```

The first derivative with respect to `log(alpha)` adds:

```text
-(log(alpha) - log(alpha_prior_mean)) / prior_var
```

The second derivative with respect to `log(alpha)` adds:

```text
-1 / prior_var
```

For the NB likelihood kernel and Cox-Reid adjustment, the Rust implementation
exposes both first and second derivatives on the log-alpha scale. The
second-derivative functions follow the chain-rule shape used by DESeq2's
`d2log_posterior`:

```text
d2 objective / d log(alpha)^2 =
  d2 objective / d alpha^2 * alpha^2
  + d objective / d log(alpha)
```

The default current gene-wise estimator still runs without a prior, matching
DESeq2's gene-wise estimate stage. Prior-aware objective, derivative,
curvature, line-search, and grid functions are used by the MAP dispersion
stage; the low-level versions also accept normalized observation weights.

## Parametric Dispersion Trend

The implemented trend foundation follows DESeq2's
`parametricDispersionFit` form:

```text
dispersion(mean) = asymptDisp + extraPois / mean
```

The fitting path starts with DESeq2's coefficients `(0.1, 1.0)`, keeps rows
selected by:

```text
dispGeneEst > 100 * minDisp
```

and applies DESeq2's robust residual screen during the outer loop:

```text
1e-4 < dispGeneEst / fittedDisp < 15
```

For the retained rows, the Rust implementation fits the Gamma GLM with identity
link for:

```text
dispGeneEst ~ 1 + I(1 / baseMean)
```

using the equivalent IRLS weighted least-squares update with weights
`1 / fittedDisp^2`. The outer loop stops when the squared log coefficient
change is below `1e-6`, matching DESeq2's control flow. If all gene-wise
estimates are within two orders of magnitude of `minDisp`, the trend fit
returns an explicit error, matching DESeq2's guidance to use gene-wise
estimates directly. Offline DESeq2 fixtures check both fitted training-row
values and predictions at means outside the original fit rows.

The mean trend type is also implemented. It follows
`estimateDispersionsFit(fitType="mean")`: first require at least one row with
`dispGeneEst > 100 * minDisp`, then compute one constant fitted dispersion from
the trimmed mean of finite estimates with `dispGeneEst > 10 * minDisp`
(`trim = 0.001` by default). The fitted value is expanded back to every finite
positive-base-mean row, with missing rows represented as `NaN`. An offline
DESeq2 fixture checks the separate viability and mean-inclusion thresholds plus
the constant fitted value.

The local trend type has an initial pure-Rust implementation. It follows
DESeq2's local-fit data contract by fitting on `log(dispGeneEst)` versus
`log(baseMean)`, retaining rows with `dispGeneEst >= 10 * minDisp`, weighting
rows by `baseMean`, and returning `minDisp` when all estimates are near the
minimum. Instead of calling R's `locfit`, Rust evaluates a deterministic
adaptive local polynomial smoother with the DESeq2-style default span. A small
offline DESeq2 local-trend fixture checks the fitted shape, with a second
fixture covering mixed rows above and below the `10 * minDisp` fit threshold
and another fixture checking predictions at means outside the original fit
rows. Fitted rows are sorted once on log mean; each prediction uses the
adaptive nearest-neighbor window directly rather than allocating and sorting a
full distance vector per gene. A separate offline fixture covers DESeq2's
all-near-minimum local floor branch, where the helper returns the
minimum-dispersion vector directly rather than a prediction function. Exact
`locfit` numerical identity and glmGamPoi trend support remain future work.

## Dispersion Prior Variance

DESeq2 estimates the variance of the normal prior on log dispersions from
residuals around the fitted trend:

```text
residual_i = log(dispGeneEst_i) - log(dispFit_i)
use_i = dispGeneEst_i >= 100 * minDisp
varLogDispEsts = mad(residual_i for use_i)^2
```

The MAD uses R's default constant `1.4826`. For model matrices with residual
degrees of freedom greater than three, DESeq2 subtracts the expected sampling
variance of log dispersion estimates:

```text
expVarLogDisp = trigamma((n_samples - n_coefficients) / 2)
dispPriorVar = max(varLogDispEsts - expVarLogDisp, 0.25)
```

For saturated designs where `n_samples == n_coefficients`, no sampling variance
is subtracted and `dispPriorVar = varLogDispEsts`.

For residual degrees of freedom one through three, DESeq2 uses a seeded
Monte-Carlo histogram/KL matching branch because the residual distribution is
asymmetric. `rsdeseq2` implements the same branch shape with deterministic
quasi-random chi-square and normal samples, the same `-10..10` histogram range
with `0.5`-wide bins, a `0..8` candidate variance grid, and local-linear
smoothing over a fine grid. This preserves deterministic Rust runs without
depending on R's random-number generator internals. Offline fixtures check the
shared residual selection, robust variance, expected sampling variance, and
bounded prior-variance range. Exact numerical identity with DESeq2's seeded
Monte Carlo plus R `loess` result remains future work.
The public `estimate_dispersion_prior` entry point is a stage-shaped wrapper
over this implemented prior-variance estimator.

The public `fit_dispersion_trend` helper dispatches the implemented
`fitType` values (`parametric`, `local`, and `mean`) with DESeq2-style default
options. `glmGamPoi` remains explicitly unsupported.

## MAP Dispersions

The initial MAP stage follows DESeq2's `type="DESeq2"` flow after gene-wise
estimates, trend fitting, and prior variance are available. Low-level MAP
fitting accepts optional normalized observation weights; the high-level
linear-mu pipeline remains no-weight because DESeq2 disables `linearMu` when
weights are present. For each non-all-zero gene:

```text
dispInit = if dispGeneEst > 0.1 * dispFit then dispGeneEst else dispFit
if dispInit is missing, use dispFit
prior mean = log(dispFit)
prior variance = dispPriorVar
```

The optimizer maximizes the Cox-Reid-adjusted log posterior with the normal
log-dispersion prior. Rows whose line search reaches `maxit` fall back to the
prior-aware two-pass grid optimizer. The MAP estimate is bounded for numerical
stability:

```text
dispMAP = clamp(dispMAP, minDisp, max(10, n_samples))
```

Final dispersions then apply DESeq2's high-side outlier rule:

```text
dispOutlier = log(dispGeneEst) >
  log(dispFit) + outlierSD * sqrt(varLogDispEsts)
dispersion = if dispOutlier then dispGeneEst else dispMAP
```

All-zero genes are expanded with missing numeric values and false outlier flags.
The low-level MAP optimizer accepts optional normalized observation weights,
but the high-level native linear-mu dispersion pipeline still rejects weights
because DESeq2 switches away from the `linearMu` branch when weights are
present. The current MAP stage does not implement glmGamPoi MAP behavior,
or replacement/refitting around dispersion outliers.

## Intercept-Only Fixed-Dispersion GLM

DESeq2 has an explicit shortcut in `fitNbinomGLMs` for an intercept-only design
with effectively no beta prior. `rsdeseq2` implements this stage separately
from general IRLS, and the public `fit_irls` dispatcher uses this shortcut
when the model matrix is intercept-only and the ridge settings are eligible:

```text
beta_i = log2(mean_j(count_ij / size_factor_j))
mu_ij = size_factor_j * 2^beta_i
w_ij = (1 / mu_ij + alpha_i)^-1
betaSE_i = log2(e) * sqrt(1 / sum_j(w_ij))
hat_ij = w_ij / sum_j(w_ij)
```

For weighted fits, the normalized-count mean and working weights are multiplied
by observation weights, matching the DESeq2 shortcut.

## Fixed-Dispersion IRLS

The initial general GLM path implements DESeq2's standard-design-matrix
`fitBeta` branch. Fitting is done on natural-log scale:

```text
mu = normalization_factor * exp(X beta)
w  = mu / (1 + alpha * mu)
z  = log(mu / normalization_factor) + (count - mu) / mu
beta_new = solve(X' W X + ridge, X' W z)
```

The ridge term is represented as a diagonal matrix. The default scalar ridge is
expanded to all coefficients, and callers can supply one natural-log-scale
ridge value per coefficient, matching DESeq2's `diag(lambda)` shape after R has
converted log2-scale beta-prior values to the natural-log scale.

When observation weights are supplied, they multiply `w` and the
negative-binomial log-likelihood terms used for the convergence deviance,
matching DESeq2's low-level `fitBeta` behavior.

DESeq2 floors fitted means to `minmu` inside the IRLS loop because the working
weights contain `1 / mu`. After the final beta is estimated,
`fitNbinomGLMs` recomputes the returned `mu` matrix from `NF * exp(X beta)`
without that floor and uses that returned matrix for `logLike` and Cook's
distance. `rsdeseq2` follows this split: floored means drive IRLS weights,
hat diagonals, and beta covariance, while unfloored final means are stored in
`NbinomGlmFit.mu` and used for row log likelihoods.

## Beta Prior Variance

`estimate_beta_prior_variance()` implements the primitive variance-estimation
part of DESeq2's beta-prior workflow for already-fitted MLE beta matrices. For
each non-intercept coefficient, it drops non-finite betas and betas with
absolute value at or above the configured finite-beta cutoff, then matches the
selected upper absolute-beta quantile to a zero-centered Normal:

```text
priorVar = (quantile(abs(beta), 1 - upperQuantile) /
            qnorm(1 - upperQuantile / 2))^2
```

The unweighted method uses R type-7 quantiles. The weighted method follows the
Hmisc weighted-quantile algorithm vendored inside DESeq2 and uses row weights
matching DESeq2's `estimateBetaPriorVar` input shape:

```text
weight_i = 1 / (1 / baseMean_i + dispFit_i)
```

Intercept columns named `Intercept` or `(Intercept)` receive the configured
wide prior variance. The full high-level beta-prior refit and expanded
model-matrix workflow are still tracked as future parity work.

`fit_glms_with_beta_prior_variance()` performs the primitive fixed-dispersion
refit from a supplied `betaPriorVar` vector. The vector is expressed on the
log2 beta scale, matching DESeq2, then converted to the natural-log ridge used
by the Rust IRLS solver:

```text
lambda_log2 = 1 / betaPriorVar
lambda_natural_log = lambda_log2 / log(2)^2
```

`fit_glms_with_estimated_beta_prior_variance()` combines the two primitive
steps for fixed-dispersion matrices: first fit MLE betas, estimate the
coefficient prior variances, and then refit with the converted per-coefficient
ridge. Expanded-model prior averaging remains separate future work.

## Observation-Weight Preprocessing

DESeq2's user-facing weighted workflow normalizes each gene's weights before
fitting:

```text
weights_i = weights_i / max(weights_i)
```

Then it checks whether the weighted design still permits parameter estimation.
For full-rank model matrices, the Rust helper follows the two DESeq2 checks:

```text
rank(weights_i * X) == ncol(X)
rank(X[weights_i > weightThreshold, nonzero_columns]) == ncol(nonzero_columns)
```

For rank-deficient designs, it follows DESeq2's fallback shape and checks that
no design column is entirely zero after weighting. Rows that fail are returned
as `weights_fail`; higher-level pipelines can treat them like DESeq2 marks
`mcols(dds)$allZero = TRUE` for failed-weight rows. The current helper is
deterministic and primitive-matrix based.

Builder stages expose this preprocessing in the fit state. With a supplied
design, raw observation weights are first used for `baseMean` and `baseVar`,
matching `getBaseMeansAndVariances`; weights are then row-normalized,
`weights_fail` is stored, and the working `all_zero` flags are computed as raw
all-zero rows OR failed-weight rows. Supplied-dispersion Wald/LRT pipelines
pass the compacted normalized weights into the weighted fixed-dispersion GLM
kernels. The GLM-mu gene-wise dispersion, MAP, and native Wald path also
passes compacted normalized weights through the fixed-dispersion mean fit,
fixed-mean dispersion objective, MAP objective, and final Wald GLM. The
high-level native linear-mu dispersion pipeline still rejects observation
weights because DESeq2 switches away from `linearMu` when weights are present.

The default Rust solver preserves the original normal-equations path for
current `useQR=FALSE` references. `IrlsSolver::Qr` follows DESeq2's `useQR=TRUE`
update by solving the augmented least-squares system:

```text
A = [sqrt(W) X; sqrt(ridge)]
b = [sqrt(W) z; 0]
beta_new = solve_qr(A, b)
```

Post-fit hat diagonals and beta covariance still use the DESeq2 formula based
on `(X' W X + ridge)^-1`. Builder pipelines copy the per-gene covariance rows
into `DeseqFit.beta_covariance` so contrast diagnostics can be validated from
the top-level fit state.

The convergence check follows DESeq2:

```text
abs(dev - dev_old) / (abs(dev) + 0.1) < betaTol
```

After IRLS, Rust applies the same row-routing predicate DESeq2 uses before its
optim backup: non-finite beta rows, non-positive coefficient variances, and,
when `IrlsOptions.use_optim` is enabled, non-converged rows are fallback
candidates. Routed rows are refit with a bounded pure-Rust optimizer on the
log2 coefficient scale, using the same `[-30, 30]` coefficient bounds and
normal-prior penalty shape as DESeq2's backup path. The refit stores optimized
betas even when the optimizer does not declare convergence, recomputes fitted
means, standard errors, coefficient covariance, and row log likelihoods, and
sets `beta_converged` from the optimizer convergence flag. `force_optim` sends
every row through the bounded refit after IRLS.

Returned beta estimates and standard errors are converted to log2 scale with
`log2(e)`.

## Wald Statistics

The implemented default coefficient-level Wald helper mirrors DESeq2's
`useT=FALSE` path:

```text
stat = beta / betaSE
pvalue = 2 * pnorm(abs(stat), lower.tail = FALSE)
```

When configured for DESeq2's `useT=TRUE` path, p-values instead use:

```text
pvalue = 2 * pt(abs(stat), df = df_i, lower.tail = FALSE)
```

The Rust options support residual degrees of freedom
`n_samples - n_coefficients`, one scalar df recycled over genes, or one df per
gene. Non-positive or non-finite df values produce missing p-values while
preserving the Wald statistic. In builder pipelines, all-zero rows have missing
t degrees of freedom after DESeq2-style row expansion.

## Wald Linear Contrasts

Implemented GLM fits store each gene's coefficient covariance matrix on the
log2 scale, and builder pipelines expose the matrix through
`DeseqFit.beta_covariance`. A primitive contrast helper accepts an explicit
numeric contrast vector `c` and computes:

```text
contrastEstimate = c' beta
contrastSE       = sqrt(c' Sigma c)
contrastStat     = contrastEstimate / contrastSE
```

This mirrors the numeric core of DESeq2's `getContrast` path after R has
resolved a contrast into numeric coefficients. The precomputed numeric contrast
output can be assembled into DESeq2-shaped result rows with BH adjustment. The
supplied-dispersion Wald builder and the native linear-mu/GLM-mu Wald builders
can run a primitive numeric contrast through GLM fitting, Cook's cutoff
masking, and independent filtering.

The primitive contrast resolver also supports named forms when a
`DesignMatrix` has coefficient names:

```text
name:              condition_B_vs_A -> unit vector for that coefficient
list contrast:     c("coef_a", "coef_b") vs c("coef_c") -> +1,+1,-1
factor level:      condition, B, A -> condition_B_vs_A
expanded factors:  condition, B, A -> conditionB - conditionA
```

Explicit list values are supported for the positive and negative lists,
matching DESeq2's `listValues` shape after names have been resolved. The
factor-level helper covers common DESeq2 coefficient-name shapes:
`factor_level_vs_reference`, the reversed reference comparison with sign flip,
two non-reference levels when an explicit reference level is supplied, and
expanded no-intercept coefficients such as `conditionB - conditionA`. It also
uses a minimal R-like `make.names` pass for common non-syntactic level names.
When named contrast specifications are run through the supplied-dispersion or
native linear-mu/GLM-mu Wald builders, the result-table metadata preserves a
stable result name and comparison label for coefficient, list, and factor-level
contrast requests.

This is still a coefficient-name resolver over an already-built design matrix.
Full DESeq2 `results(contrast=...)` compatibility still needs colData/formula
ownership, all coefficient cleanup rules, and contrast-aware Cook's/refit
edge-case handling.

For primitive numeric contrasts with both positive and negative coefficients,
the supplied-dispersion and native linear-mu/GLM-mu Wald contrast pipelines
mirror DESeq2's `contrastAllZeroNumeric` behavior. Non-zero contrast
coefficients are converted to one, samples are selected by non-zero
`modelMatrix %*% contrastBinary`, and genes with zero raw counts across those
selected samples get
`log2FoldChange = 0`, `stat = 0`, and `pvalue = 1`, unless the gene is already
an all-zero row. One-sided numeric contrasts leave rows unchanged, matching
DESeq2.

For factor-level contrasts, a separate primitive helper mirrors DESeq2's
`contrastAllZeroCharacter` rule when the caller supplies one sample level per
column. Samples whose level is either the numerator or denominator are selected;
genes with zero raw counts across all selected samples get the same
`log2FoldChange = 0`, `stat = 0`, and `pvalue = 1` override in the
factor-level Wald contrast builder. R still owns formula parsing and colData
factor validation.

## Wald LFC Thresholds

The same selected-coefficient path supports DESeq2's LFC-threshold alternatives
from `results(lfcThreshold=..., altHypothesis=...)` for primitive matrix
workflows:

```text
greaterAbs:     |beta| > T
greaterAbs2014: older 2014 two-sided threshold test
greaterAbsUPSHOT: UP-SHOT threshold test, normal p-values only
lessAbs:        |beta| < T
greater:        beta > T
less:           beta < -T
```

These formulas operate on log2-scale beta estimates and log2-scale standard
errors, matching DESeq2 result columns. `lessAbs` requires a positive
threshold, and `greaterAbsUPSHOT` with t p-values is explicitly unsupported,
matching DESeq2's current guard.

The helper operates on an `NbinomGlmFit` and a selected coefficient index. The
builder pipelines wrap it with result-table assembly, Cook's filtering, and
independent filtering. Primitive numeric contrasts have their own result-row
builder and supplied-dispersion Wald pipeline, but full Bioconductor result
metadata remains future work.

## Result Rows

The initial result assembly mirrors the simple non-contrast DESeq2 `results()`
path, with Cook's p-value masking and base-mean independent filtering available
in the builder pipelines:

```text
baseMean
log2FoldChange
lfcSE
stat
pvalue
padj = p.adjust(pvalue, method = "BH")
```

Optional diagnostic fields include dispersion, beta convergence, Cook's
distances, `maxCooks`, Cook's outlier flags, and independent-filtering flags
when those stages are run. `DeseqResults::column_names()` reports the core
DESeq2 result columns plus whichever optional diagnostic columns are present in
the table. `DeseqResults::deseq2_metadata()` reports a primitive metadata view:
test type, reported coefficient or contrast, Wald threshold and alternative,
DESeq2-style comparison-aware column descriptions, p-value adjustment method,
and independent-filtering metadata when present. The table metadata also
exposes scalar key/value entries and the resolved effect-size and
test-statistic description labels: Wald contrast tables prefer the comparison
label when present, while LRT tables keep the reported full-model coefficient
as the effect label and use the model comparison for statistic and p-value
labels. `DeseqResults::data_frame()`
assembles the same rows into a typed data-frame view for wrappers and file
writers: gene identifiers are kept as row names, numeric and logical columns
are separate typed vectors, and each column carries the DESeq2-style metadata
description. `io::write_deseq_results_tsv()` exports the assembled result table
with a first `gene` column, DESeq2-shaped result columns, and R-style `NA`
values for missing numeric or logical entries. `io::write_deseq_results_tidy_tsv()`
exports the same table with the first column named `row`, matching the
`results(tidy = TRUE)` shape.
`io::write_deseq_result_column_metadata_tsv()` exports the corresponding
`mcols(res)`-style column metadata with `name`, `type`, and `description`
columns. `io::write_deseq_result_table_metadata_tsv()` exports table-level
result metadata as `name`/`value` entries.
`DeseqFit::deseq2_mcols_diagnostics()` provides a DESeq2-name-shaped view for
implemented diagnostics such as `dispGeneEst`, `dispGeneIter`, `dispFit`,
`dispersion`, `dispIter`, `dispOutlier`, `betaConv`, `fullBetaConv`,
`reducedBetaConv`, `betaIter`, `deviance`, and `maxCooks`. The diagnostics
view also reports the present column names in stable stage order for wrappers
and parity-table exporters, and can assemble those fields into a typed
diagnostic data-frame view. `io::write_deseq_mcols_diagnostics_tsv()` exports
the same fields with a leading `gene` column and R-style `NA` values.
Full Bioconductor result objects, formula-aware contrast metadata, and complete
wrapper metadata preservation are future work.

## Cook's Distances

The current implementation mirrors DESeq2's `calculateCooksDistance` shape for
the supplied-dispersion Wald path and the limited native linear-mu Wald path:

```text
V_ij = mu_ij + alpha_i^robust * mu_ij^2
PearsonResSq_ij = (K_ij - mu_ij)^2 / V_ij
cooks_ij = PearsonResSq_ij / p * H_ij / (1 - H_ij)^2
```

Here `p` is the number of model-matrix columns and `H` is the fitted hat
diagonal. The Cook's variance uses DESeq2's robust method-of-moments
dispersion with a minimum of `0.04`; it is separate from the fitted/supplied
dispersion. `maxCooks` is recorded only when `n_samples > p` and at least one
model-matrix cell has three or more replicates. The samples included in
`maxCooks` follow the same model-matrix-cell rule as DESeq2's `nOrMoreInCell`.

## Cook's Cutoff

For the current Wald result rows, Cook's cutoff filtering can be enabled,
disabled, or set to an explicit threshold through `CooksCutoff`. The default
cutoff is DESeq2's:

```text
qf(.99, p, m - p)
```

where `p` is the number of model-matrix columns and `m` is the number of
samples. Rows with `maxCooks > cutoff` have `pvalue` set to missing before
Benjamini-Hochberg adjustment is recomputed. This implements the basic
`results()` p-value masking behavior.

DESeq2 applies an additional low-count heuristic for formula designs with one
two-level factor: if a row is above the Cook's cutoff, find the sample with the
maximum Cook's distance, take its raw count, and do not mask the row if at
least three counts in that row are larger. `rsdeseq2` exposes this as
`apply_cooks_cutoff_with_low_count_heuristic`. The helper is explicit because
the Rust core receives primitive matrices and cannot infer R formula/colData
semantics on its own.

The primitive replacement-count transform from DESeq2 `replaceOutliers()` is
also available as a low-level helper. It computes:

```text
trimBaseMean_i = mean(counts_i / factor_i, trim = 0.2)
replacement_ij = as.integer(trimBaseMean_i * factor_ij)
```

where `factor_ij` is either the sample size factor or the gene/sample
normalization factor. Only cells with `cooks_ij > cutoff` in replaceable samples
are changed. Replaceable samples are either supplied explicitly or selected from
model-matrix cells with at least `minReplicates` samples, matching DESeq2's
`whichSamples`/`nOrMoreInCell` shape. The helper returns transformed counts,
candidate replacement counts, per-cell outlier flags, per-gene `replace` flags,
and the replaceable-sample mask.

When `m <= p` (`m` samples and `p` model-matrix columns), replacement is skipped,
matching the early return in DESeq2 `replaceOutliers()`.

A second helper prepares the `refitWithoutOutliers` bookkeeping without running
the refit: it recomputes `baseMean`, `baseVar`, and `allZero` on replacement
counts, records `nrefit = sum(replace, na.rm=TRUE)`, separates `refitReplace`
rows from `newAllZero` rows, and applies DESeq2's post-refit `maxCooks` rule
where original Cook's distances in replaceable sample columns are zeroed before
recording the maximum. If every sample is replaceable, post-refit `maxCooks` is
missing for every gene.

The first end-to-end replacement-refit paths are implemented for the current
GLM-mu native Wald, Wald contrast, and LRT branches. They run the original fit
with Cook's filtering and independent filtering disabled, prepare replacement
counts from the original Cook's distances, refit the implemented GLM-mu
dispersion/MAP path on the replacement counts while preserving the original
size factors, then rerun the relevant Wald, Wald contrast, or LRT test, merge
only `refitReplace` rows into the result table, clear result fields for
`newAllZero` rows, and finally apply the caller's Cook's cutoff and independent
filtering to the merged result rows.

Automatic R formula/colData dispatch for the two-group low-count heuristic,
remaining contrast edge cases, beta-prior behavior, and Bioconductor
assay/metadata preservation are not implemented yet.

## Independent Filtering

The current result path implements DESeq2's default independent-filtering
shape with `baseMean` as the filter statistic:

```text
lowerQuantile = mean(baseMean == 0)
upperQuantile = if lowerQuantile < 0.95 then 0.95 else 1.0
theta = seq(lowerQuantile, upperQuantile, length.out = 50)
cutoff_t = quantile(baseMean, theta_t)
padj_t = BH(pvalue for rows with baseMean >= cutoff_t)
numRej_t = count(padj_t < alpha)
```

If the maximum number of rejections is no greater than 10, the first threshold
is selected, matching DESeq2's guard against over-aggressive filtering when all
genes are null. Otherwise, the threshold is selected using the same
max-smoothed-fit-minus-RMSE rule as DESeq2, with an R `stats::lowess`-shaped
smoother over the rejection curve. The smoother uses `f = 1/5`, R's
floor-based span, sliding contiguous fitting window, tricube distance weights,
local-linear weighted least squares, Tukey robustifying iterations, and the
default `delta = 0.01 * diff(range(theta))` interpolation shortcut. This path is
tested against R-generated fixtures for both the default 50-point DESeq2 theta
grid and a dense custom theta grid where `delta` skips fitted points. The result
metadata records the theta grid, rejection counts, selected theta, selected
filter threshold, lowess fitted values, and alpha. It also exposes paired
`filterNumRej`- and `lo.fit`-shaped views with theta/rejection-count and
theta/smoothed-rejection columns for wrappers and parity exporters. Scalar
metadata is exposed as key/value entries for `filterThreshold`, `filterTheta`,
and `alpha`, omitting optional entries when filtering is disabled. TSV
exporters write those metadata tables using DESeq2-style `theta`/`numRej`,
`x`/`y`, and `name`/`value` column names. Rows below the selected threshold keep
their raw p-value but receive a missing adjusted p-value.

## Wald/LRT Pipelines

The supplied-dispersion Wald pipeline accepts counts, a design matrix, caller
supplied dispersions, and either a selected coefficient or a primitive numeric
contrast. Size factors can either be estimated with the configured size-factor
method or supplied directly by the caller. Gene/sample normalization factors
can also be supplied; when present, they preempt size factors for normalized
counts and GLM offsets. The pipeline performs:

```text
size-factor estimation
normalized counts
baseMean/baseVar/allZero
fixed-dispersion NB GLM fit
all-zero row expansion for GLM outputs
Cook's distance matrix and maxCooks diagnostics
Cook's cutoff p-value masking
baseMean independent filtering
selected-coefficient Wald test
selected-coefficient result rows
```

For primitive numeric contrasts, the same pipeline uses `c' beta` and
`sqrt(c' Sigma c)` after the GLM fit, then assembles contrast result rows
before applying Cook's cutoff masking and independent filtering.

The limited native Wald pipeline reuses the same GLM/result path after running
the implemented dispersion subset:

```text
size-factor estimation
baseMean/baseVar/allZero
linear-mu gene-wise dispersions
or GLM-mu gene-wise dispersions
parametric, local, or mean dispersion trend
deterministic dispersion prior variance
MAP dispersions
fixed-dispersion NB GLM fit using final dispersions
Cook's diagnostics and selected-coefficient Wald results
```

This follows the DESeq2 `DESeq()` stage order at a high level:
`estimateSizeFactors`, `estimateDispersions`, then `nbinomWaldTest`. The
implemented Rust path remains narrower than DESeq2 because it only supports
the current linear-mu no-weight or GLM-mu optionally weighted,
parametric/local/mean-trend, deterministic-prior branches.
`DeseqBuilder::fit()` and `DeseqBuilder::fit_with_results()` now expose the
implemented GLM-mu Wald branch as a limited top-level workflow and report the
last design coefficient by default, matching DESeq2's default result-coefficient
shape. LRT remains available through `DeseqBuilder::fit_lrt()`,
`DeseqBuilder::fit_lrt_with_results()`, and the branch-specific `fit_lrt_*`
entry points because callers must provide a reduced design.
Gene/sample normalization factors are supported in this subset by following
DESeq2's `linearModelMuNormalized`: projected normalized counts are multiplied
by `getSizeOrNormFactors`, and moments starts use
`mean(1 / colMeans(normalizationFactors))`.

The limited native LRT pipeline uses the same implemented dispersion subset,
with the full design driving dispersion estimation. It then fits full and
reduced fixed-dispersion GLMs using the final MAP dispersions and stores the
same full/reduced diagnostics as the supplied-dispersion LRT path. Linear-mu
and GLM-mu branches are available with parametric, local, or mean trends; observation
weights are limited to the GLM-mu branch.

The LRT path accepts full and reduced design matrices, fits both models on the
non-all-zero subset, and computes:

```text
stat_i = 2 * (logLik_full_i - logLik_reduced_i)
pvalue_i = pchisq(stat_i, df = ncol(full) - ncol(reduced), lower.tail = FALSE)
```

`DeseqFit` stores the full-model GLM diagnostics through the standard beta,
standard-error, covariance, convergence, iteration, log-likelihood,
DESeq2-style full deviance (`-2 * logLike`), fitted-mean, and hat-diagonal
fields. LRT pipelines also store reduced-model log likelihood, beta
convergence, and beta iteration vectors for stage-by-stage parity checks.
`DeseqFit::deseq2_mcols_diagnostics()` also exposes the stable fitted
dispersion trend type label when a parametric, local, or mean trend has been
attached, keeping transform and result parity logs on the same label set.

The R reference generator validates this current scope with
`DESeq2:::fitNbinomGLMs` using supplied dispersions, default `1e-6` beta ridge,
`useQR=FALSE`, and `useOptim=FALSE`. Full DESeq2 parity still requires broader
native dispersion validation, optim fallback, and richer
result metadata.

## Gene-Wise Dispersion Foundation

The current dispersion stage implements reusable pieces of DESeq2's
`estimateDispersionsGeneEst` path for both the `linearMu` branch and an
initial GLM-mu branch.

The rough starting estimate follows `roughDispEstimate`:

```text
mu = max(1, linearModelMu(normalized_counts, X))
rough_i = sum_j (((y_ij - mu_ij)^2 - mu_ij) / mu_ij^2) / (m - p)
rough_i = max(rough_i, 0)
```

The moments estimate follows `momentsDispEstimate`:

```text
xim = mean_j(1 / size_factor_j)
moments_i = (baseVar_i - xim * baseMean_i) / baseMean_i^2
```

With gene/sample normalization factors, the DESeq2 branch is:

```text
xim = mean_j(1 / colMeans_i(normalization_factor_ij))
moments_i = (baseVar_i - xim * baseMean_i) / baseMean_i^2
```

All-zero rows are excluded from the normalization-factor column means, matching
the `objectNZ` subset used before DESeq2 calls `momentsDispEstimate`.

The bounded starting value is:

```text
alpha_init_i = clamp(min(rough_i, moments_i), minDisp, maxDisp)
maxDisp = max(10, n_samples)
```

For the initial Rust estimator, fitted raw means are computed from projected
normalized counts and size factors or normalization factors:

```text
mu_raw_ij = max(minmu, linearModelMu(normalized_counts, X)_ij * size_factor_j)
mu_raw_ij = max(minmu, linearModelMu(normalized_counts, X)_ij * normalization_factor_ij)
```

For the GLM-mu branch, non-all-zero rows alternate between fixed-dispersion NB
GLM mean fitting and fixed-mean dispersion optimization:

```text
alpha_hat = clamp(min(rough, moments), minDisp, maxDisp)
fitidx = non-all-zero rows

for iter in seq_len(niter):
  mu[fitidx, ] = fitNbinomGLMs(alpha_hat[fitidx])$mu
  alpha_hat_new[fitidx] = fitDisp(mu[fitidx, ], alpha_hat[fitidx])
  fitidx = abs(log(alpha_hat_new) - log(alpha_hat)) > 0.05
  alpha_hat = alpha_hat_new
```

This mirrors DESeq2's `niter` and `fitidx` shape. The current Rust branch uses
the existing fixed-dispersion IRLS implementation for `fitNbinomGLMs`, applies
the same `minmu` floor, and can consume builder-supplied normalized
observation weights. It feeds the same parametric/local/mean trend, prior-variance,
MAP shrinkage, and native Wald stages as the linear-mu branch.

Then each non-all-zero gene is optimized on the log-alpha scale. By default,
the score includes DESeq2's Cox-Reid adjustment. Observation weights multiply
the likelihood terms; for the Cox-Reid determinant they define the
`weightThreshold` sample subset, matching `fitDisp`:

```text
score(alpha) = ll_alpha_dependent(counts_i, mu_i, alpha)
               - 0.5 * log(det(X_i' W_i X_i))
X_i = X[obs_weight_i > weightThreshold, nonzero_columns]
W_ij = (1 / mu_ij + alpha)^-1 for selected samples
```

The alpha-dependent likelihood kernel follows DESeq2's `log_posterior`
objective and omits constants independent of alpha:

```text
ll = sum_j [
  lgamma(y_ij + 1/alpha)
  - lgamma(1/alpha)
  - y_ij * log(mu_ij + 1/alpha)
  - (1/alpha) * log(1 + mu_ij * alpha)
]
```

The default optimizer follows DESeq2's Armijo line-search shape:

```text
a = log(alpha_init)
lp = score(a)
dlp = d score(a) / d a
kappa = kappa_0

proposal = a + kappa * dlp
accept when -score(proposal) <= -lp - kappa * 1e-4 * dlp^2
on accept:
  a = proposal
  stop if score(a) - lp < dispTol
  update dlp
  kappa = min(kappa * 1.1, kappa_0)
  every five accepts, halve kappa
on reject:
  halve kappa
```

Rows that do not converge and remain above `minDisp * 10` fall back to the
two-pass grid search, matching DESeq2's fallback structure:

```text
grid_1 = seq(log(minDisp), log(maxDisp), length.out = 20)
best_1 = argmax_alpha score(alpha)
grid_2 = seq(best_1 - delta, best_1 + delta, length.out = 20)
dispGeneEst_i = clamp(exp(argmax grid_2), minDisp, maxDisp)
```

All-zero genes are expanded with missing numeric values (`NaN` internally),
following DESeq2's missing-row expansion convention. The implemented MAP path
adds the log-dispersion prior, deterministic prior variance, and weighted
optimizer support. The GLM-mu builder path can pass normalized observation
weights through gene-wise dispersion, MAP, and native Wald/LRT stages; the current
native linear-mu branch is still no-weight. The line-search diagnostics now
record the final first and second derivatives (`last_dlp` and `last_d2lp`) for
DESeq2-style curvature checks, and high-level fit states retain the gene-wise
iteration counts as `dispGeneIter`-compatible diagnostics. The current weighted
GLM-mu mean and local-trend MAP/Wald/LRT branches have DESeq2 golden
references, including compact result-table and BH-adjusted p-value checks for
the matched result rows. The unweighted GLM-mu Cox-Reid local-trend branch has
a MAP dispersion anchor. Broader full parity still requires exact local
`locfit` identity, glmGamPoi trend support, and remaining edge-case coverage.
The local smoother now treats a single usable local-fit row as a constant
log-dispersion trend, matching the generated tiny GLM-mu local fixture instead
of failing on zero tricube neighborhood weight.

The full-model beta and standard error for a selected coefficient are reported
alongside the model-level LRT statistic and p-value, matching the shape of
DESeq2 results after `nbinomLRT`.

This mirrors the shape of DESeq2's `nbinomWaldTest` and `nbinomLRT` when
dispersion estimates already exist. The native Wald path now supplies those
estimates from the limited linear-mu or GLM-mu parametric/local/mean MAP
subsets.
Cook's outlier replacement/refitting, beta priors, exact R `lowess` parity,
and full contrast result-table handling are not implemented. Shared internal
helpers skip all-zero rows during GLM fitting and expand compact outputs back
to full gene order with `NaN` internal numeric values and `None` result-table
values, matching the intent of DESeq2's `buildMatrixWithNARows` and
`buildDataFrameWithNARows` helpers without copying their implementation.
