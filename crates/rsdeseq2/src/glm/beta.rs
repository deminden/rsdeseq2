use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::NbinomGlmFit;
use crate::glm::irls::{
    IrlsOptions, fit_fixed_dispersion_irls,
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights, fit_irls,
};
use crate::glm::nb::nbinom_log_likelihood_matrix;
use crate::matrix::RowMajorMatrix;
use statrs::distribution::{ContinuousCDF, Normal};

/// Method used to estimate DESeq2's beta prior variance.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BetaPriorVarianceMethod {
    /// Match the upper absolute-beta quantile without mean/dispersion weights.
    Quantile,
    /// Match the upper absolute-beta quantile using DESeq2's
    /// `1 / (1 / baseMean + dispFit)` weights.
    Weighted,
}

/// Options for DESeq2-style beta prior variance estimation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorVarianceOptions {
    /// Estimation method. DESeq2 defaults to `weighted`.
    pub method: BetaPriorVarianceMethod,
    /// Upper tail probability matched against a zero-centered Normal.
    pub upper_quantile: f64,
    /// Wide prior variance used for intercepts and columns with no finite betas.
    pub wide_prior_variance: f64,
    /// Absolute beta cutoff used to discard near-infinite MLEs.
    pub finite_beta_cutoff: f64,
}

impl Default for BetaPriorVarianceOptions {
    fn default() -> Self {
        Self {
            method: BetaPriorVarianceMethod::Weighted,
            upper_quantile: 0.05,
            wide_prior_variance: 1e6,
            finite_beta_cutoff: 10.0,
        }
    }
}

/// Two-stage fixed-dispersion GLM output for a DESeq2-style beta prior refit.
#[derive(Clone, Debug, PartialEq)]
pub struct BetaPriorGlmFit {
    /// First pass MLE GLM fit used to estimate beta prior variances.
    pub mle_fit: NbinomGlmFit,
    /// Refit GLM using beta-prior variance as a per-coefficient ridge.
    pub prior_fit: NbinomGlmFit,
    /// Log2-scale beta prior variances, one per design coefficient.
    pub beta_prior_variance: Vec<f64>,
}

/// Expanded-model fixed-dispersion GLM output for a DESeq2-style beta-prior refit.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedModelBetaPriorGlmFit {
    /// First pass MLE GLM fit on the expanded design.
    pub expanded_mle_fit: NbinomGlmFit,
    /// Refit GLM on the expanded design using beta-prior variance as ridge.
    pub expanded_prior_fit: NbinomGlmFit,
    /// Prior fit collapsed onto the caller-supplied standard design matrix.
    pub prior_fit: NbinomGlmFit,
    /// Log2-scale beta prior variances, one per expanded-design coefficient.
    pub beta_prior_variance: Vec<f64>,
}

/// Expanded and reported design matrices for a beta-prior expanded-model refit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExpandedModelBetaPriorDesignInput<'a> {
    /// Design used for expanded-model fitting.
    pub expanded_design: &'a DesignMatrix,
    /// Reported standard design used after grouped coefficient collapse.
    pub standard_design: &'a DesignMatrix,
    /// Expanded coefficient columns averaged into each standard coefficient.
    pub coefficient_groups: &'a [Vec<usize>],
}

/// Options for the two-stage beta-prior fixed-dispersion refit helper.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BetaPriorRefitOptions {
    /// Options used for the first MLE GLM fit and the final prior refit.
    ///
    /// The prior refit replaces any ridge settings with values derived from
    /// the estimated beta prior variances.
    pub fit_options: IrlsOptions,
    /// Options used when estimating beta prior variances from the MLE fit.
    pub variance_options: BetaPriorVarianceOptions,
}

/// Size-factor offsets and optional observation weights for beta-prior refits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorSizeFactorWeightInput<'a> {
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
}

/// Normalization-factor offsets and optional observation weights for beta-prior refits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorNormalizationFactorWeightInput<'a> {
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
}

/// Estimate fixed-dispersion beta coefficients with DESeq2-style GLM dispatch.
///
/// This is a public beta-estimation convenience entry point for callers that
/// already have size factors and per-gene dispersions. Intercept-only designs
/// use the closed-form DESeq2 shortcut through `fit_irls`; other designs use
/// the general fixed-dispersion IRLS implementation.
pub fn estimate_beta(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_irls(counts, design, size_factors, dispersions, options)
}

/// Convert log2-scale beta prior variances to natural-log IRLS ridge values.
///
/// DESeq2 computes `lambda = 1 / betaPriorVar` on the log2 beta scale and then
/// divides by `log(2)^2` before fitting on the natural-log scale. This helper
/// exposes that conversion for primitive Rust GLM refits.
pub fn beta_prior_variance_to_ridge_lambda(
    beta_prior_variance: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    beta_prior_variance
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, variance)| {
            validate_positive_finite(variance, "beta prior variance", idx)?;
            let inv_ln2 = std::f64::consts::LOG2_E;
            Ok(variance.recip() * inv_ln2 * inv_ln2)
        })
        .collect()
}

// Expanded-model coefficient, covariance and contrast collapse.
include!("beta/expanded.rs");

// Fixed and estimated beta-prior fitting workflows.
include!("beta/fitting.rs");

// Prior-variance estimation and upper-quantile matching.
include!("beta/prior.rs");

// Intercept-only fits, normalization and working weights.
include!("beta/intercept.rs");

// Expanded-fit options, validation and covariance handling.
include!("beta/expanded_validation.rs");

// Prior input checks, stable arithmetic and quantiles.
include!("beta/prior_math.rs");

// Fit and normalization-factor validation.
include!("beta/validation.rs");

// Private numerical regression tests.
include!("beta/tests.rs");
