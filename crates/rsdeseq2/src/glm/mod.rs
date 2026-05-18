//! Future negative-binomial GLM fitting.

pub mod beta;
pub mod dispersion_fit;
pub mod fallback;
pub mod irls;
pub mod lrt;
pub mod nb;
pub mod wald;
pub mod weights;

pub use beta::{
    fit_intercept_only_fixed_dispersion,
    fit_intercept_only_fixed_dispersion_with_normalization_factors,
    fit_intercept_only_fixed_dispersion_with_weights,
};
pub use irls::{
    fit_fixed_dispersion_irls, fit_fixed_dispersion_irls_with_normalization_factors,
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights, IrlsOptions, IrlsSolver,
};
pub use lrt::lrt_test;
pub use nb::{
    nbinom_log_likelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
    nbinom_log_pmf, nbinom_negative_twice_log_likelihood,
};
pub use wald::{
    two_sided_normal_pvalue, two_sided_t_pvalue, wald_stat_and_pvalue,
    wald_stat_and_pvalue_with_options, wald_test_coefficient, wald_test_coefficient_with_options,
    wald_test_contrast, wald_test_contrast_with_options, WaldAlternative, WaldContrastOutput,
    WaldDegreesOfFreedom, WaldPvalueType, WaldTestOptions,
};
pub use weights::{
    preprocess_observation_weights, preprocess_observation_weights_with_options,
    ObservationWeightOptions, ObservationWeights,
};

use crate::design::DesignMatrix;
use crate::matrix::RowMajorMatrix;

/// Negative-binomial GLM fit output matching DESeq2 low-level result fields.
#[derive(Clone, Debug, PartialEq)]
pub struct NbinomGlmFit {
    /// Per-gene log likelihood.
    pub log_like: Vec<f64>,
    /// Per-gene beta convergence flags.
    pub beta_converged: Vec<bool>,
    /// Beta estimates on log2 scale, matching DESeq2 returned `betaMatrix`.
    pub beta: RowMajorMatrix<f64>,
    /// Beta standard errors on log2 scale.
    pub beta_se: RowMajorMatrix<f64>,
    /// Per-gene beta covariance matrices on log2 scale.
    ///
    /// Stored as genes x `(n_terms * n_terms)`, with each gene row containing a
    /// row-major coefficient covariance matrix.
    pub beta_covariance: Option<RowMajorMatrix<f64>>,
    /// Fitted mean matrix.
    pub mu: RowMajorMatrix<f64>,
    /// Per-gene beta iteration counts.
    pub beta_iter: Vec<usize>,
    /// Model matrix used for fitting.
    pub model_matrix: DesignMatrix,
    /// Number of model terms.
    pub n_terms: usize,
    /// Hat diagonal matrix.
    pub hat_diagonal: RowMajorMatrix<f64>,
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
