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
log-likelihoods matching DESeq2's `nbinomLogLike`. Individual log-PMF terms,
weighted terms, and accumulated row log-likelihood sums are checked for finite
values before they are returned to GLM, Wald, or LRT stages.

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
Observation-weighted fast-VST trend fitting reuses the design-preprocessed
weight matrix and subsets those rows alongside counts and normalization data.
`DeseqBuilder::vst_glm_mu_auto` now performs the DESeq2-shaped automatic
choice: use the deterministic fast subset when enough rows are eligible and
carry any preprocessed observation weights through that subset, or fit the
selected Rust trend on all rows when the default subset is too large for the
dataset. The returned `VstGlmMuOutput` records whether the trend came from the
fast subset or full data, the requested `nsub`, and the eligible-row count;
`VstTrendSource` also exposes accessor helpers for these fields and for the
fast-subset decision. Full-data trend metadata records whether that path was
chosen because too few rows were eligible. The source and full-data reason
enums expose stable labels: `fastSubset`, `fullData`, and
`insufficientEligibleRows`. `VstGlmMuOutput::metadata()` packages those labels
with `nsub`, eligible-row count, transform dimensions, trend-fit row and sample
counts, stable trend-fit type label, optional fast-subset row count, and
optional original fast-subset row indices for wrappers and benchmark logs.
`blind_vst_glm_mu_auto` uses the same automatic choice with a named
intercept-only design, matching the implemented part of DESeq2's `blind=TRUE`
workflow without invoking any external runtime.
`CountMatrix::select_rows` provides the matching raw-count subset while
preserving gene and sample names, which mirrors the object subset passed to the
fast-VST dispersion fit. Full fast-VST parity still requires object metadata
and exact local interpolation semantics.

## Regularized Log

The Rust core includes a low-level regularized-log sample-effect primitive.
`rlog_with_size_factors()` and `rlog_with_normalization_factors()` build the
same foundational design shape used by DESeq2 rlog: one intercept plus one
indicator column per sample. The intercept receives a wide log2-scale prior
variance, while every sample effect receives a shared caller-supplied
log2-scale prior variance. The implementation reuses the existing
negative-binomial beta-prior GLM machinery and returns a genes x samples matrix
assembled as:

```text
rlog_ij = intercept_i + sample_effect_ij
```

`estimate_rlog_sample_prior_variance()` implements the rlog-specific prior
estimate from already-normalized counts, `baseMean`, and `dispFit`:

```text
logFoldChange_ij = log2(normalized_count_ij + 0.5) - log2(baseMean_i + 0.5)
weight_i = 1 / (1 / baseMean_i + dispFit_i)
```

The flattened log-fold-change values and repeated row weights are matched to a
zero-centered Normal using DESeq2's default weighted upper-tail quantile rule.
`rlog_with_estimated_prior_and_size_factors()` and
`rlog_with_estimated_prior_and_normalization_factors()` compose that
prior-estimation step with the sample-effect fit when the caller already has
`baseMean`, `dispFit`, and gene-wise dispersions from earlier stages.
`rlog_fit_with_size_factors()` and `rlog_fit_with_normalization_factors()`
retain the fitted intercept-plus-sample-effect GLM beside the transformed
matrix, exposing the intermediate beta surface needed for stricter parity
diagnostics and future frozen-rlog reuse.
`rlog_frozen_with_size_factors()` and
`rlog_frozen_with_normalization_factors()` provide the matching low-level
frozen-intercept transform: caller-supplied log2 intercepts are converted to
gene-specific offsets, only sample-effect coefficients are fit, and the output
is assembled as the frozen intercept plus the fitted sample effect.
`DeseqFit::regularized_log_transform()` and the short `rlog()` alias expose
the same composition from stored fit state after the implemented dispersion
MAP stages have produced `dispFit` and final dispersions.
`DeseqFit::regularized_log_transform_with_frozen_intercept()` and the short
`frozen_rlog()` alias expose frozen-intercept reuse from the same fit state,
using stored final dispersions and offsets. For one-call Rust workflows,
`DeseqBuilder::rlog_glm_mu()` first runs the implemented GLM-mu MAP dispersion
stages and then applies fit-state rlog. `blind_rlog_glm_mu()` uses a named
intercept-only design for the same implemented `blind=TRUE` shape used by the
CLI. `frozen_rlog_glm_mu_with_fit()` and
`blind_frozen_rlog_glm_mu_with_fit()` learn the rlog intercept/prior, then run
a frozen-intercept reuse pass from the same dispersion fit. The `*_with_fit`
variants return `RlogGlmMuOutput` or `FrozenRlogGlmMuOutput`, retaining the MAP
dispersion fit state beside the rlog matrices for wrappers, diagnostics, and
parity checks. Fit-state and builder-level rlog skip all-zero rows during the
sample-effect fit and re-expand them as zero-valued transform rows. Fit-state
frozen rlog uses the supplied frozen intercept for all-zero rows, so ordinary
unfiltered count matrices remain accepted. `RlogOutput` includes the
transformed matrix, fitted per-gene intercepts, estimated sample-effect prior
variance, the offset source, and a compact metadata view for wrappers and
benchmark logs.
This is still not the full high-level DESeq2 `rlog()` object workflow. The
low-level and builder-level frozen-intercept numeric surfaces are present,
while full object dispatch and Bioconductor metadata remain future work.

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
change is below `1e-6`, matching DESeq2's control flow. The Rust path rejects
non-finite weighted least-squares terms, determinant and coefficient-numerator
products, Gamma deviances, coefficient changes, and fitted trend values before
they can propagate into later dispersion stages. If all gene-wise estimates are
within two orders of magnitude of `minDisp`, the trend fit returns an explicit
error, matching DESeq2's guidance to use gene-wise estimates directly. Offline
DESeq2 fixtures check both fitted training-row values and predictions at means
outside the original fit rows.

The mean trend type is also implemented. It follows
`estimateDispersionsFit(fitType="mean")`: first require at least one row with
`dispGeneEst > 100 * minDisp`, then compute one constant fitted dispersion from
the trimmed mean of finite estimates with `dispGeneEst > 10 * minDisp`
(`trim = 0.001` by default). Rust exposes those two row-selection masks
separately, because DESeq2 uses the stricter threshold only as the preliminary
viability gate. The fitted value is expanded back to every finite
positive-base-mean row, with missing rows represented as `NaN`. The trimmed mean
uses stable averaging so very large finite dispersion estimates do not overflow
before division. An offline DESeq2 fixture checks the separate viability and
mean-inclusion thresholds plus the constant fitted value.

The local trend type has a pure-Rust locfit-compatible implementation. It follows
DESeq2's local-fit data contract by fitting on `log(dispGeneEst)` versus
`log(baseMean)`, retaining rows with `dispGeneEst >= 10 * minDisp`, weighting
rows by `baseMean`, and returning `minDisp` when all estimates are near the
minimum. Instead of calling R's `locfit`, Rust uses a pure-Rust compatibility
backend for tricube local polynomial fits. The default DESeq2-sized path uses
the backend's locfit-style Hermite evaluation grid; smaller or custom-span
local trends use the same backend's direct local prediction path. A small
offline DESeq2 local-trend fixture checks the fitted shape, with a second
fixture covering mixed rows above and below the `10 * minDisp` fit threshold
and another fixture checking predictions at means outside the original fit
rows. The real-data local trend fixture currently checks 64,344 finite fitted
values with median relative error `3.74e-13`, p99 `5.85e-12`, and max
`1.47e-11` against DESeq2 1.46.0. Compared with the previous in-repo smoother,
that real-data fixture moved from median relative error `7.99e-03`, p99
`2.00e-01`, and max `4.28e-01` to the near-exact values above. The newest
locfit-compatible backend revision tightened the same real-data fixture from
median relative error `4.04e-10`, p99 `2.80e-09`, and max `3.19e-09` to the
current values. On the committed small local-trend fixtures, the fitted-shape
max absolute error improved from `3.29e-04` to `9.62e-05`, and the
mixed-threshold max absolute error improved
from `3.87e-02` to `1.98e-03`. The small out-of-fit prediction fixture moved
from max absolute error `5.76e-05` to `2.56e-04`; this remains inside its
`2e-3` DESeq2 reference tolerance, but it is not counted as an improvement.
The committed GLM-mu local MAP/Wald/LRT fixtures were already at machine
precision and are unchanged by this smoother swap.
Manually constructed local trends reuse the
builder shape checks for span and polynomial degree, then validate finite log
dispersions, positive finite local weights, and a consistent empty state for the
minimum-dispersion floor branch before prediction. Batch trend evaluation
prevalidates the fitted trend once before expanding finite positive rows and
missing `NaN` rows, so malformed manual trend state is reported even when every
requested mean is missing. A separate offline fixture covers DESeq2's
all-near-minimum local floor branch, where the helper returns the
minimum-dispersion vector directly rather than a prediction function. Broader
synthetic locfit edge cases and glmGamPoi trend support remain future work.

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

The MAP path checks initial-value threshold arithmetic, optimizer-produced
dispersion values, and outlier-threshold arithmetic before those values drive
clamping or final dispersion replacement.

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
wide prior variance.

The expanded-model beta-prior surface now has primitive helpers for the core
coefficient and covariance arithmetic used after fitting an expanded model
matrix. Given a log2-scale expanded beta matrix and groups of expanded columns,
the collapse helper returns per-gene group averages:

```text
beta_collapsed_ij = mean(beta_expanded_i, columns(group_j))
```

For per-gene expanded covariance matrices, the covariance helper applies the
same averaging matrix:

```text
Sigma_collapsed_i = A Sigma_expanded_i A'
```

The companion contrast helper builds averaged numerator-vs-denominator vectors
with `+1 / n` weights over numerator columns and `-1 / m` weights over
denominator columns. A fit-level collapse helper packages the same arithmetic
back into an `NbinomGlmFit`: collapsed betas, propagated covariance, standard
errors recomputed from the collapsed covariance diagonal, the caller-supplied
standard design matrix, and the original fitted means and diagnostics.
Result-table helpers can then build DESeq2-shaped Wald rows for a selected
collapsed coefficient or for a caller-supplied numeric contrast on the
collapsed standard-design coefficient scale.

The expanded refit helper runs the same sequence end to end for primitive
matrices: MLE fit on the expanded design, prior-variance estimation from the
expanded MLE coefficients, ridge refit on the expanded design, and collapse
onto the caller-supplied standard design. Result-table companions can consume
that expanded beta-prior fit output directly for selected coefficients or
numeric contrasts. A fit-and-results workflow helper combines those pieces for
the common Wald coefficient or numeric-contrast path, including size-factor,
normalization-factor, and optional observation-weight inputs. A primitive
one-factor design helper builds the expanded intercept-plus-level-indicator
matrix, the matching treatment-style reported design, and coefficient groups
from caller-supplied sample labels. One-factor fit-and-results helpers now own
that construction and pass the generated expanded design, standard design, and
column groups into the same Wald coefficient or numeric-contrast workflows.
For additive designs, a multi-term helper covers primitive
`~ factor1 + factor2 + numeric1 + factor1:factor2 + factor1:numeric1 +
numeric1:numeric2 + ...` construction: intercept, one expanded indicator per
factor level, treatment-style reported columns for non-reference levels,
numeric covariates included unchanged in both design surfaces, primitive
pairwise interaction columns, and matching collapse groups. Additive
fit-and-results helpers run the same coefficient and numeric-contrast
workflows with size-factor, normalization-factor, and optional
observation-weight inputs. A primitive formula helper now parses a
DESeq2-style subset (`1` intercept-only, `+`, `:`, `/`, `*` shorthand, and
`0`/`-1` intercept removal), including lower-order-omitted pairwise
interactions, higher-order interaction/nesting/star-expansion terms, and
primitive `- term` subtraction for the same supported term subset. Additive
parenthesized groups, including nested additive groups, distribute through
`*`, `:`, `/`, and subtraction in that same primitive subset. Integer numeric
power transforms such as `I(dose^2)` and common numeric function transforms
`log(numeric)`, `log2(numeric)`, `log10(numeric)`, `log1p(numeric)`,
`sqrt(numeric)`, and `scale(numeric)` are materialized as derived numeric
covariates with sanitized coefficient names. `scale(numeric)` supports boolean
or scalar numeric
`center=` and `scale=` arguments and follows R's single-column default:
centered scaling divides by sample standard deviation, while uncentered
scaling divides by root mean square.
Raw polynomial transforms `poly(numeric, degree, raw=TRUE)` are materialized as
`numeric_poly_1` through `numeric_poly_degree`, and the generated additive
group participates in supported main-effect, interaction, shorthand, and
subtraction expansion. Orthogonal `poly()` remains outside the native formula
subset.
`offset(numeric)` terms and single-vector supported transform offsets such as
`offset(log2(numeric))` or `offset(I(numeric + other_numeric))` are extracted
by `expanded_formula_design_with_offsets()` into per-sample log-offset
vectors, with multiple offsets summed sample-wise. Formula-driven
fit-and-results helpers exponentiate those log offsets and multiply them into either
sample-level size factors expanded across genes or supplied gene/sample
normalization factors, then run the same expanded beta-prior Wald coefficient
and numeric-contrast workflows from the parsed design. Splines, arbitrary R
expressions, orthogonal polynomial bases, and complete R-compatible formula
semantics remain future work; the runtime numeric path is pure Rust and expects
callers or wrappers to provide or derive unsupported design surfaces explicitly.

`fit_glms_with_beta_prior_variance()` performs the primitive fixed-dispersion
refit from a supplied `betaPriorVar` vector. Size-factor, normalization-factor,
and optional observation-weight entry points share the same conversion. The
vector is expressed on the log2 beta scale, matching DESeq2, then converted to
the natural-log ridge used by the Rust IRLS solver:

```text
lambda_log2 = 1 / betaPriorVar
lambda_natural_log = lambda_log2 / log(2)^2
```

`fit_glms_with_estimated_beta_prior_variance()` combines the two primitive
steps for fixed-dispersion matrices: first fit MLE betas, estimate the
coefficient prior variances, and then refit with the converted per-coefficient
ridge. Companion helpers cover the same estimated-prior workflow for
normalization-factor offsets and optional observation weights.

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

The default Rust solver follows DESeq2's `useQR=TRUE` update by solving the
augmented least-squares system. `IrlsSolver::NormalEquations` remains available
for older fixed-dispersion `useQR=FALSE` references:

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
candidates. Routed rows are refit with a mature pure-Rust L-BFGS-B port on the
log2 coefficient scale, using the same `[-30, 30]` coefficient bounds,
100-iteration default, five correction pairs, and normal-prior penalty shape as
DESeq2's backup path. The optimizer now uses the same R-compatible wrapper
shape as DESeq2's `optim(..., method="L-BFGS-B")` fallback: objective-only calls
with R-style finite-difference gradients, default `factr=1e7`, `pgtol=0`, and
no extra post-optimizer polish. Failed high-gradient excursions still fall back
to the optimizer start instead of storing an unstable move. The refit stores
optimized betas even when the optimizer does not declare convergence, recomputes fitted
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
masking, and independent filtering. Non-finite contrast estimates, covariance
quadratic forms, or standard errors produce missing contrast statistics for the
affected gene instead of propagating non-finite values downstream.

The primitive contrast resolver also supports named forms when a
`DesignMatrix` has coefficient names:

```text
name:              condition_B_vs_A -> unit vector for that coefficient
list contrast:     c("coef_a", "coef_b") vs c("coef_c") -> +1,+1,-1
factor level:      condition, B, A -> condition_B_vs_A
expanded factors:  condition, B, A -> conditionB - conditionA
```

Explicit list values are supported for the positive and negative lists,
matching DESeq2's `listValues` shape after names have been resolved: the first
value must be greater than zero and the second value must be less than zero. The
factor-level helper covers common DESeq2 coefficient-name shapes:
`factor_level_vs_reference`, the reversed reference comparison with sign flip,
two non-reference levels when an explicit reference level is supplied or can be
inferred from shared coefficient suffixes, and expanded no-intercept
coefficients such as `conditionB - conditionA`. It also uses a minimal R-like
`make.names` pass for common non-syntactic level names.
When named contrast specifications are run through the supplied-dispersion or
native linear-mu/GLM-mu Wald builders, the result-table metadata preserves a
stable result name and comparison label for coefficient, list, and factor-level
contrast requests. Coefficient-list comparison labels follow DESeq2's
`cleanContrast` shape for two-sided, positive-only, and negative-only list
contrasts.

The core resolver operates over an already-built design matrix, while the R
wrapper can use fitted `factorLevels` or factor-valued `colData` to infer the
reference level for character triplet contrasts. Primitive Rust, CLI, and R
already-fitted result routes now cover DESeq2-style character, list/listValues,
and numeric `results(contrast=...)` semantics, including `contrast` taking
precedence over `name`. Stored `DeseqFit` objects expose a single
`results_contrast_request()` dispatcher that selects Wald or LRT behavior from
the fitted test state while preserving branch-specific contrast semantics and
contrast-all-zero cleanup for character and two-sided numeric/list contrasts.
Remaining work is high-level Bioconductor object mutation/plumbing and broader
contrast-aware Cook's/refit edge-case handling.

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
factor-level Wald contrast builder. The R wrapper can derive the reference from
fitted factor metadata for already-fitted result objects; full high-level
`DESeqDataSet` mutation remains an interface boundary.

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
independent filtering. Primitive numeric and DESeq2-style results contrasts
have result-row builders and supplied/native Wald pipeline routes. Full
Bioconductor result-object metadata remains interface work.

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
label when present, while LRT tables use the reported full-model coefficient
or contrast as the effect label and keep the model comparison for statistic
and p-value labels. `DeseqResults::data_frame()`
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
Cook's replacement/refit plans expose a scalar metadata summary with refit
counts, new-all-zero rows, outlier/replaced cell counts, replaceable sample
counts, and the refit-branch decision; `io::write_cooks_replacement_metadata_tsv()`
exports those entries as a `name`/`value` table.
Cook's diagnostic outputs can also be exported directly:
`io::write_cooks_distance_tsv()` writes the per-gene, per-sample Cook's
distance assay with `NA` for missing values,
`io::write_cooks_row_metadata_tsv()` writes `maxCooks` and robust dispersion,
and `io::write_cooks_sample_metadata_tsv()` writes the sample eligibility mask
used to record `maxCooks`.
The same plan can export replacement-count assays through
`io::write_cooks_replaced_counts_tsv()`,
`io::write_cooks_candidate_replacement_counts_tsv()`, and
`io::write_cooks_outlier_cells_tsv()`, plus row-level replacement/refit
metadata through `io::write_cooks_replacement_row_metadata_tsv()`.
`DeseqFit::deseq2_mcols_diagnostics()` provides a DESeq2-name-shaped view for
implemented diagnostics such as `dispGeneEst`, `dispGeneIter`, `dispFit`,
`dispMAP`, `dispersion`, `dispIter`, `dispOutlier`, `betaConv`, `fullBetaConv`,
`reducedBetaConv`, `betaIter`, `deviance`, and `maxCooks`. For Rust-side
fallback debugging, the same view also includes `rustBetaOptimIter` and
`rustBetaOptimStartObjective`/`rustBetaOptimObjective` plus
`rustBetaOptimGradientNorm` when GLM fitting state is present. The diagnostics
view reports the present column names in stable stage order for wrappers and
parity-table exporters, and can assemble those fields into a typed diagnostic
data-frame view.
`io::write_deseq_mcols_diagnostics_tsv()` exports the same fields with a
leading `gene` column and R-style `NA` values. Matrix diagnostics such as
fitted means and full/reduced hat diagonals remain on `DeseqFit`; they are not
flattened into the `mcols(dds)`-style row-metadata view.
The native Wald and LRT CLI commands can export the same row diagnostics with
`--fit-diagnostics-output`; replacement/refit runs can also export the refit
branch with `--refit-diagnostics-output`. They can also export full-model beta
and beta standard-error matrices with `--fit-beta-output`,
`--fit-beta-se-output`, `--fit-beta-optim-start-output`,
`--refit-beta-output`, `--refit-beta-se-output`, and
`--refit-beta-optim-start-output`. These sidecars are intended for
stage-by-stage parity checks, especially when result-table differences need to
be traced back to `dispGeneEst`, `dispFit`, `dispMAP`, final dispersion
estimates, beta estimates, beta optimizer starts, beta standard errors,
iteration/convergence fields, deviance, or Cook's summaries.
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
counts from the original Cook's distances, estimate gene-wise dispersions on
the replacement counts, reuse the original dispersion function and dispersion
prior variance for MAP shrinkage, then rerun the relevant Wald, Wald contrast,
or LRT test, merge only `refitReplace` rows into the result table, clear result
fields for `newAllZero` rows, and finally apply the caller's Cook's cutoff and
independent filtering to the merged result rows. Top-level GLM-mu result
helpers expose this limited replacement-refit path for the default Wald/LRT
coefficient and
primitive numeric, named/list, and caller-supplied factor-level Wald contrasts.

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
shape. `fit_with_results_name()` exposes the selected-coefficient path by
design coefficient name, using the same exact-first R-cleaned aliases as named
primitive contrasts. `fit_contrast()`,
`fit_contrast_spec()`, and
`fit_factor_level_contrast()` expose fit-only top-level Wald contrast helpers;
`fit_with_results_contrast()`, `fit_with_results_contrast_spec()`, and
`fit_with_results_factor_level_contrast()` expose the same branch with result
rows for primitive numeric, named/list, and caller-supplied factor-level Wald
contrasts. Compatibility-named parametric Wald helpers expose the same
contrast shapes while pinning the dispersion trend to the parametric branch.
Builder-created `DeseqFit` objects retain optional formula/model-frame
metadata, including factor levels and references, for later wrapper or
already-fitted-object Wald/LRT `results(contrast=...)` assembly through the
unified fitted-object contrast dispatcher.
The corresponding `*_with_cooks_replacement()` helpers expose the limited
replacement-refit workflow when callers supply replacement options; formula
contrast-request replacement helpers use the same exact-first R-cleaned
coefficient aliases as the non-replacement result routes. Stored formula
model-frame metadata also lets selected-coefficient Wald replacement refits
reuse the two-group low-count Cook's gate when the selected coefficient is the
non-reference comparison for one two-level factor. LRT
remains available through `DeseqBuilder::fit_lrt()`,
`DeseqBuilder::fit_lrt_name()`,
`DeseqBuilder::fit_lrt_contrast()`,
`DeseqBuilder::fit_lrt_contrast_spec()`,
`DeseqBuilder::fit_lrt_factor_level_contrast()`,
`DeseqBuilder::fit_lrt_with_results()`,
`DeseqBuilder::fit_lrt_with_results_name()`,
`DeseqBuilder::fit_lrt_with_results_contrast()`,
`DeseqBuilder::fit_lrt_with_results_contrast_spec()`,
`DeseqBuilder::fit_lrt_with_results_factor_level_contrast()`,
`DeseqBuilder::fit_lrt_with_results_with_cooks_replacement()`,
`DeseqBuilder::fit_lrt_with_results_name_with_cooks_replacement()`,
`DeseqBuilder::fit_lrt_with_results_factor_level_contrast_with_cooks_replacement()`,
and the
branch-specific `fit_lrt_*` entry points because callers must provide a
reduced design. LRT contrast result tables report the requested full-model
effect size, but keep the full-vs-reduced LRT statistic and p-value; when
DESeq2's numeric or character/factor-level `contrastAllZero` rule applies,
only the reported `log2FoldChange` is zeroed for LRT rows.
Gene/sample normalization factors are supported in this subset by following
DESeq2's `linearModelMuNormalized`: projected normalized counts are multiplied
by `getSizeOrNormFactors`, and moments starts use
`mean(1 / colMeans(normalizationFactors))`.

The limited native LRT pipeline uses the same implemented dispersion subset,
with the full design driving dispersion estimation. It then fits full and
reduced fixed-dispersion GLMs using the final MAP dispersions and stores the
same full/reduced diagnostics as the supplied-dispersion LRT path. Linear-mu
and GLM-mu branches are available with parametric, local, or mean trends; observation
weights are limited to the GLM-mu branch. Both native branches can report
numeric, named/list, or caller-supplied factor-level full-model effect sizes
without changing the likelihood-ratio test itself. The compatibility-named
parametric LRT entry points provide the same contrast routes while pinning the
dispersion trend to the parametric branch regardless of the builder's
configured `fit_type`.

The supplied-dispersion LRT path mirrors the implemented contrast result
surface for numeric, named/list, and caller-supplied factor-level full-model
effect sizes. Factor-level LRT contrasts use caller-provided sample levels for
DESeq2-style character `contrastAllZero` cleanup, again zeroing only the
reported LFC while preserving the likelihood-ratio statistic and p-values.

The LRT path accepts full and reduced design matrices, fits both models on the
non-all-zero subset, and computes:

```text
stat_i = 2 * (logLik_full_i - logLik_reduced_i)
pvalue_i = pchisq(stat_i, df = ncol(full) - ncol(reduced), lower.tail = FALSE)
```

For LRT result tables, selected coefficients and primitive numeric or named
full-model contrasts affect only the reported `log2FoldChange` and `lfcSE`
columns. The `stat`, `pvalue`, and `padj` columns remain tied to the same
full-vs-reduced likelihood-ratio test.

`DeseqFit` stores the full-model GLM diagnostics through the standard beta,
standard-error, covariance, convergence, iteration, log-likelihood,
DESeq2-style full deviance (`-2 * logLike`), fitted-mean, and hat-diagonal
fields. LRT pipelines also store reduced-model log likelihood, beta
convergence, beta iteration, fitted-mean, and hat-diagonal fields for
stage-by-stage parity checks. The full/reduced fitted-mean and hat-diagonal
matrices stay as matrix-valued `DeseqFit` state rather than `mcols(dds)` row
metadata columns.
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

The shared bounded optimizer and the GLM beta fallback reject non-finite norm,
dot-product, movement, directional-derivative, and objective-building
accumulations. Overflow in optimizer control quantities now
produces a non-converged optimizer result for the affected row; overflow while
building the beta objective, gradient, Hessian, or ridge penalty is converted
into the same large finite penalty used for invalid fitted means. Fitted-mean reconstruction
uses checked row-wise linear predictors so overflow is reported on the affected
sample before exponentiation.
The intercept-only fixed-dispersion shortcut also checks normalized mean,
weighted mean, reconstructed fitted-mean, and intercept-information sums, and
uses the stable `(1 / mu + alpha)^-1` working-weight formula.
Weighted beta-prior quantile matching checks total, normalized, merged, and
cumulative weights before interpolation so extreme observation-derived weights
cannot silently overflow the shrinkage-scale estimate.
Observation-weight preprocessing checks weighted design values and Cox-Reid
subset column sums for finite arithmetic, while the rank-deficient zero-column
path avoids unnecessary products when testing whether a column is represented.
Gene-wise dispersion Cox-Reid matrix construction checks design cross-products,
weighted derivative matrix terms, selected-column sums, and trace accumulations
before determinant, derivative, or second-derivative values are used by the
dispersion optimizer.
Rough and moments dispersion starts check residual squares, inverse
size-factor summaries, normalization-factor column means, and per-gene moment
terms. The dispersion likelihood objective, derivative, and second derivative
also check per-sample terms, observation-weight products, and row-level
accumulations before line-search or grid-search decisions consume them.
Ratio and poscounts size-factor estimation check per-gene geometric-mean log
sums before median ratio construction, and frozen size-factor stabilization
checks the final log-size-factor sum before rescaling to geometric mean one.
Base mean and base variance metadata use checked row sums, checked
observation-weight products, and checked centered sum-of-squares before those
values are stored on the fitted state.
Result assembly validates base means before building Wald or LRT result rows,
so invalid upstream metadata is reported before adjusted p-values or
diagnostic columns are emitted.
Precomputed Wald/LRT outputs are also checked for finite statistics, bounded
p-values, valid t degrees of freedom, finite contrast estimates and standard
errors, and aligned reduced-model convergence flags before result tables are
built.
Later mutable result-table operations reuse the same p-value contract before
recomputing BH adjustments, and Cook's filtering rejects non-finite `maxCooks`
diagnostics before masking p-values.
DESeq2-shaped diagnostic exports validate aligned column lengths and gene names
before writing row-wise metadata.
Size-factor and base-mean exports validate name alignment and finite value
domains before writing scalar sample or gene metadata.
Result-table exports validate manually assembled rows for finite effects,
bounded p-values, positive dispersions, and non-negative Cook's diagnostics
before serializing TSV output.
Independent-filtering metadata exports validate theta domains, alpha, selected
indices, lowess length, and finite rejection-curve values before writing.
Wald p-value helpers clamp finite tail probabilities to `[0, 1]` and preserve
missing p-values for non-finite tail calculations.
Cook's default cutoff validates F-distribution degrees of freedom and the
resulting quantile before it can mask result p-values.
LRT result construction checks log-likelihood subtraction and deviance scaling
before chi-square p-values are computed, preserving missing values for rows
whose likelihood difference cannot be represented finitely.
Finite LRT chi-square upper-tail probabilities are clamped to `[0, 1]`, and
reduced-model convergence flags are validated before being copied into LRT
diagnostics.
Stored full-model deviance diagnostics use the same finite guard around
`-2 * logLike`, turning non-representable values into the existing missing
numeric representation.
The public negative-binomial `-2 * logLik` helper also checks the final
deviance scaling after the row log-likelihood has been accumulated.
Cook's-distance robust dispersion and outlier replacement helpers check
trimmed means, residual squares, variance factors, and Cook's-distance
arithmetic so extreme finite values are reported before replacement or refit
planning consumes them.
The local VST path checks inverse size-factor summaries, normalization-factor
geometric-mean summaries, row-mean accumulation, variance-curve products, and
integration/interpolation arithmetic before producing transformed values.
Independent-filtering lowess checks residual, robust-weight, local-linear
center/variance, interpolation, fitted-value, RMSE, and threshold accumulations
before selecting the final filter threshold.
MAP shrinkage checks initial-threshold arithmetic, optimizer output, and
outlier-threshold arithmetic before clamping or replacing final dispersions.
Numeric `contrastAllZero` checks design-score accumulation before selecting
samples for the all-zero contrast rule.
The low-residual-df dispersion-prior branch also checks local smoother weighted
regression sums and prediction arithmetic, so extreme finite KL grids fail
explicitly instead of choosing a spurious prior-variance minimum.

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
the matched result rows. The unweighted and weighted GLM-mu Cox-Reid
local-trend branches have MAP/Wald/LRT anchors, including fitted means, hat
diagonals, log-likelihoods, deviances, convergence, and compact result rows.
Broader full parity still requires glmGamPoi trend support and remaining
weighted-dispersion and synthetic smoother edge-case coverage.
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
Cook's outlier replacement/refitting, beta priors, exact R `lowess` parity, and
contrast result-table handling are available on the implemented workflow
surfaces described above, with remaining gaps concentrated in high-level
Bioconductor object attachment and edge-case coverage. Shared internal helpers
skip all-zero rows during GLM fitting and expand compact outputs back to full
gene order with `NaN` internal numeric values and `None` result-table values,
matching the intent of DESeq2's `buildMatrixWithNARows` and
`buildDataFrameWithNARows` helpers without copying their implementation.
