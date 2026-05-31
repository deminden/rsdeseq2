//! Negative-binomial GLM fitting primitives.

pub mod beta;
pub mod dispersion_fit;
pub mod fallback;
pub mod irls;
pub mod lrt;
pub mod nb;
pub mod wald;
pub mod weights;

pub use beta::{
    average_expanded_model_coefficients, average_expanded_model_covariances,
    beta_prior_variance_to_ridge_lambda, collapse_expanded_model_fit, estimate_beta,
    estimate_beta_prior_variance, expanded_model_group_contrast,
    fit_expanded_glms_with_estimated_beta_prior_variance,
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors,
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights,
    fit_expanded_glms_with_estimated_beta_prior_variance_and_weights,
    fit_glms_with_beta_prior_variance, fit_glms_with_beta_prior_variance_and_normalization_factors,
    fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights,
    fit_glms_with_beta_prior_variance_and_weights, fit_glms_with_estimated_beta_prior_variance,
    fit_glms_with_estimated_beta_prior_variance_and_normalization_factors,
    fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights,
    fit_glms_with_estimated_beta_prior_variance_and_weights, fit_intercept_only_fixed_dispersion,
    fit_intercept_only_fixed_dispersion_with_normalization_factors,
    fit_intercept_only_fixed_dispersion_with_weights, match_upper_quantile_for_variance,
    match_weighted_upper_quantile_for_variance, BetaPriorGlmFit,
    BetaPriorNormalizationFactorWeightInput, BetaPriorRefitOptions, BetaPriorSizeFactorWeightInput,
    BetaPriorVarianceMethod, BetaPriorVarianceOptions, ExpandedModelBetaPriorDesignInput,
    ExpandedModelBetaPriorGlmFit,
};
pub use dispersion_fit::fit_with_dispersion;
pub use fallback::{optim_fallback_rows, OptimFallbackRows};
pub use irls::{
    fit_fixed_dispersion_irls, fit_fixed_dispersion_irls_with_normalization_factors,
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights, fit_irls, IrlsOptions, IrlsSolver,
};
pub use lrt::lrt_test;
pub use nb::{
    nbinom_log_likelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
    nbinom_log_pmf, nbinom_negative_twice_log_likelihood,
};
pub use wald::{
    two_sided_normal_pvalue, two_sided_t_pvalue, wald_stat_and_pvalue,
    wald_stat_and_pvalue_with_options, wald_test, wald_test_coefficient,
    wald_test_coefficient_with_options, wald_test_contrast, wald_test_contrast_with_options,
    WaldAlternative, WaldContrastOutput, WaldDegreesOfFreedom, WaldPvalueType, WaldTestOptions,
};
pub use weights::{
    preprocess_observation_weights, preprocess_observation_weights_with_options,
    ObservationWeightOptions, ObservationWeights,
};

use crate::design::DesignMatrix;
use crate::matrix::RowMajorMatrix;

/// Negative-binomial GLM fit output matching DESeq2 low-level result fields.
#[derive(Clone, Debug)]
pub struct NbinomGlmFit {
    /// Per-gene log likelihood.
    pub log_like: Vec<f64>,
    /// Per-gene beta convergence flags.
    pub beta_converged: Vec<bool>,
    /// Beta estimates on log2 scale, matching DESeq2 returned `betaMatrix`.
    pub beta: RowMajorMatrix<f64>,
    /// Beta standard errors on log2 scale.
    pub beta_se: RowMajorMatrix<f64>,
    /// Fallback optimizer starting beta values on log2 scale.
    ///
    /// Rows/coefficients that did not enter the fallback optimizer are `NaN`.
    pub beta_optim_start: RowMajorMatrix<f64>,
    /// Per-gene beta covariance matrices on log2 scale.
    ///
    /// Stored as genes x `(n_terms * n_terms)`, with each gene row containing a
    /// row-major coefficient covariance matrix.
    pub beta_covariance: Option<RowMajorMatrix<f64>>,
    /// Fitted mean matrix.
    pub mu: RowMajorMatrix<f64>,
    /// Per-gene beta iteration counts.
    pub beta_iter: Vec<usize>,
    /// Rust fallback-optimizer iterations for rows routed after IRLS.
    ///
    /// Rows that did not enter the fallback optimizer are `NaN` so diagnostic
    /// TSV exports can preserve the full gene shape without implying an
    /// optimizer result.
    pub beta_optim_iter: Vec<f64>,
    /// Rust fallback-optimizer objective at the starting parameter vector.
    ///
    /// Rows that did not enter the fallback optimizer are `NaN`.
    pub beta_optim_start_objective: Vec<f64>,
    /// Final Rust fallback-optimizer objective for rows routed after IRLS.
    ///
    /// Rows that did not enter the fallback optimizer are `NaN`.
    pub beta_optim_objective: Vec<f64>,
    /// Projected gradient norm at the final Rust fallback-optimizer parameters.
    ///
    /// Rows that did not enter the fallback optimizer are `NaN`.
    pub beta_optim_gradient_norm: Vec<f64>,
    /// Model matrix used for fitting.
    pub model_matrix: DesignMatrix,
    /// Number of model terms.
    pub n_terms: usize,
    /// Hat diagonal matrix.
    pub hat_diagonal: RowMajorMatrix<f64>,
}

impl PartialEq for NbinomGlmFit {
    fn eq(&self, other: &Self) -> bool {
        nan_equal_vec(&self.log_like, &other.log_like)
            && self.beta_converged == other.beta_converged
            && nan_equal_matrix(&self.beta, &other.beta)
            && nan_equal_matrix(&self.beta_se, &other.beta_se)
            && nan_equal_matrix(&self.beta_optim_start, &other.beta_optim_start)
            && nan_equal_optional_matrix(&self.beta_covariance, &other.beta_covariance)
            && nan_equal_matrix(&self.mu, &other.mu)
            && self.beta_iter == other.beta_iter
            && nan_equal_vec(&self.beta_optim_iter, &other.beta_optim_iter)
            && nan_equal_vec(
                &self.beta_optim_start_objective,
                &other.beta_optim_start_objective,
            )
            && nan_equal_vec(&self.beta_optim_objective, &other.beta_optim_objective)
            && nan_equal_vec(
                &self.beta_optim_gradient_norm,
                &other.beta_optim_gradient_norm,
            )
            && self.model_matrix == other.model_matrix
            && self.n_terms == other.n_terms
            && nan_equal_matrix(&self.hat_diagonal, &other.hat_diagonal)
    }
}

fn nan_equal_optional_matrix(
    left: &Option<RowMajorMatrix<f64>>,
    right: &Option<RowMajorMatrix<f64>>,
) -> bool {
    match (left, right) {
        (Some(left), Some(right)) => nan_equal_matrix(left, right),
        (None, None) => true,
        _ => false,
    }
}

fn nan_equal_matrix(left: &RowMajorMatrix<f64>, right: &RowMajorMatrix<f64>) -> bool {
    left.n_rows() == right.n_rows()
        && left.n_cols() == right.n_cols()
        && nan_equal_vec(left.as_slice(), right.as_slice())
}

fn nan_equal_vec(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .copied()
            .zip(right.iter().copied())
            .all(nan_equal_f64)
}

fn nan_equal_f64((left, right): (f64, f64)) -> bool {
    left == right || (left.is_nan() && right.is_nan())
}

/// Future Wald-test output.
#[derive(Clone, Debug, PartialEq)]
pub struct WaldOutput {
    /// Per-gene Wald statistics for the selected coefficient.
    pub stat: Vec<Option<f64>>,
    /// Per-gene p-values.
    pub pvalue: Vec<Option<f64>>,
    /// Per-gene t degrees of freedom when t p-values are requested.
    pub degrees_of_freedom: Option<Vec<Option<f64>>>,
}

/// Future likelihood-ratio-test output.
#[derive(Clone, Debug, PartialEq)]
pub struct LrtOutput {
    /// Per-gene deviance differences.
    pub deviance: Vec<Option<f64>>,
    /// Per-gene p-values.
    pub pvalue: Vec<Option<f64>>,
    /// Chi-square degrees of freedom.
    pub degrees_of_freedom: usize,
    /// Reduced-model beta convergence flags.
    pub reduced_converged: Vec<bool>,
}
