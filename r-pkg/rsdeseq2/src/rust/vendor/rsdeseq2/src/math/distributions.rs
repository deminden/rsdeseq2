use crate::core::CountMatrix;
use crate::errors::DeseqError;
use crate::glm::nb::{
    nbinom_log_likelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
    nbinom_log_pmf, nbinom_negative_twice_log_likelihood,
};
use crate::matrix::RowMajorMatrix;

/// Convenience namespace for DESeq2-parameterized negative-binomial helpers.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct NegativeBinomialHelpers;

/// Return a stateless helper namespace for negative-binomial calculations.
pub fn negative_binomial_helpers() -> NegativeBinomialHelpers {
    NegativeBinomialHelpers
}

impl NegativeBinomialHelpers {
    /// Negative-binomial log PMF using DESeq2's `mu` and dispersion shape.
    pub fn log_pmf(self, count: u32, mu: f64, dispersion: f64) -> Result<f64, DeseqError> {
        negative_binomial_log_pmf(count, mu, dispersion)
    }

    /// Row log likelihood for one gene.
    pub fn log_likelihood(
        self,
        counts: &[u32],
        mu: &[f64],
        dispersion: f64,
    ) -> Result<f64, DeseqError> {
        negative_binomial_log_likelihood(counts, mu, dispersion)
    }

    /// Weighted row log likelihood for one gene.
    pub fn log_likelihood_weighted(
        self,
        counts: &[u32],
        mu: &[f64],
        dispersion: f64,
        weights: Option<&[f64]>,
    ) -> Result<f64, DeseqError> {
        negative_binomial_log_likelihood_weighted(counts, mu, dispersion, weights)
    }

    /// DESeq2-style `-2 * logLik` for one gene.
    pub fn negative_twice_log_likelihood(
        self,
        counts: &[u32],
        mu: &[f64],
        dispersion: f64,
    ) -> Result<f64, DeseqError> {
        negative_binomial_negative_twice_log_likelihood(counts, mu, dispersion)
    }
}

/// Negative-binomial log PMF using DESeq2's `mu` and dispersion parameterization.
pub fn negative_binomial_log_pmf(count: u32, mu: f64, dispersion: f64) -> Result<f64, DeseqError> {
    nbinom_log_pmf(count, mu, dispersion)
}

/// Row log likelihood for one gene.
pub fn negative_binomial_log_likelihood(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
) -> Result<f64, DeseqError> {
    nbinom_log_likelihood(counts, mu, dispersion)
}

/// Weighted row log likelihood for one gene.
pub fn negative_binomial_log_likelihood_weighted(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    nbinom_log_likelihood_weighted(counts, mu, dispersion, weights)
}

/// Row-wise log likelihoods for a count matrix.
pub fn negative_binomial_log_likelihood_matrix(
    counts: &CountMatrix,
    mu: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<Vec<f64>, DeseqError> {
    nbinom_log_likelihood_matrix(counts, mu, dispersions, weights)
}

/// DESeq2-style `-2 * logLik` for one gene.
pub fn negative_binomial_negative_twice_log_likelihood(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
) -> Result<f64, DeseqError> {
    nbinom_negative_twice_log_likelihood(counts, mu, dispersion)
}
