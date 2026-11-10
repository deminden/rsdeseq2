//! Shared component-family dispatch for standalone prior fitting and summaries.
use super::{
    NormalMixture, NormalMixturePosterior, TruncatedNormalMixture, TruncatedNormalMixturePosterior,
    UniformMixture, UniformMixturePosterior,
};
use crate::DeseqError;

/// Fixed-mixture families with a normal observation likelihood by default.
/// Uniform components may explicitly select a Student-t likelihood.
#[derive(Clone, Debug)]
pub enum FixedMixture {
    Normal(NormalMixture),
    Uniform(UniformMixture),
    TruncatedNormal(TruncatedNormalMixture),
}
/// A conditioned fixed mixture, retaining its component-specific tail semantics.
#[derive(Clone, Debug)]
pub enum FixedMixturePosterior {
    Normal(NormalMixturePosterior),
    Uniform(UniformMixturePosterior),
    TruncatedNormal(TruncatedNormalMixturePosterior),
}
impl FixedMixture {
    /// Component log convolution densities in input component order.
    pub fn log_likelihoods(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<Vec<f64>, DeseqError> {
        match self {
            Self::Normal(p) => p.log_likelihoods(estimate, standard_error),
            Self::Uniform(p) => p.log_likelihoods(estimate, standard_error),
            Self::TruncatedNormal(p) => p.log_likelihoods(estimate, standard_error),
        }
    }
    /// Condition on the same finite observation through every interface.
    pub fn posterior(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<FixedMixturePosterior, DeseqError> {
        Ok(match self {
            Self::Normal(p) => {
                FixedMixturePosterior::Normal(p.posterior(estimate, standard_error)?)
            }
            Self::Uniform(p) => {
                FixedMixturePosterior::Uniform(p.posterior(estimate, standard_error)?)
            }
            Self::TruncatedNormal(p) => {
                FixedMixturePosterior::TruncatedNormal(p.posterior(estimate, standard_error)?)
            }
        })
    }
    pub(super) fn posterior_at(
        &self,
        estimate: f64,
        standard_error: f64,
        row: usize,
        observations: usize,
    ) -> Result<FixedMixturePosterior, DeseqError> {
        if let Self::TruncatedNormal(p) = self {
            Ok(FixedMixturePosterior::TruncatedNormal(p.posterior_at(
                estimate,
                standard_error,
                row,
                observations,
            )?))
        } else {
            self.posterior(estimate, standard_error)
        }
    }
    pub(super) fn cdf_values_at(
        &self,
        estimate: f64,
        standard_error: f64,
        thresholds: &[f64],
        row: usize,
        observations: usize,
    ) -> Result<Vec<f64>, DeseqError> {
        if let Self::Uniform(p) = self {
            return p.cdf_values(estimate, standard_error, thresholds);
        }
        let p = self.posterior_at(estimate, standard_error, row, observations)?;
        Ok(thresholds.iter().map(|&c| p.cdf(c)).collect())
    }
    /// Condition a complete observation vector, retaining any family-specific
    /// batch-shape semantics required by the established numerical contract.
    pub fn posteriors(
        &self,
        estimates: &[f64],
        standard_errors: &[f64],
    ) -> Result<Vec<FixedMixturePosterior>, DeseqError> {
        if let Self::TruncatedNormal(p) = self {
            return Ok(p
                .posteriors(estimates, standard_errors)?
                .into_iter()
                .map(FixedMixturePosterior::TruncatedNormal)
                .collect());
        }
        if estimates.len() != standard_errors.len() {
            return Err(crate::errors::invalid_dimensions(
                "posterior standard errors",
                estimates.len(),
                standard_errors.len(),
            ));
        }
        estimates
            .iter()
            .zip(standard_errors)
            .map(|(&x, &s)| self.posterior(x, s))
            .collect()
    }
}
impl FixedMixturePosterior {
    /// Posterior expectation on the observation scale.
    pub fn mean(&self) -> f64 {
        match self {
            Self::Normal(p) => p.mean(),
            Self::Uniform(p) => p.mean(),
            Self::TruncatedNormal(p) => p.mean(),
        }
    }
    /// Posterior standard deviation.
    pub fn standard_deviation(&self) -> f64 {
        match self {
            Self::Normal(p) => p.standard_deviation(),
            Self::Uniform(p) => p.standard_deviation(),
            Self::TruncatedNormal(p) => p.standard_deviation(),
        }
    }
    /// Normalized component probabilities in prior input order.
    pub fn weights(&self) -> &[f64] {
        match self {
            Self::Normal(p) => p.weights(),
            Self::Uniform(p) => p.weights(),
            Self::TruncatedNormal(p) => p.weights(),
        }
    }
    /// Right-continuous posterior CDF.
    pub fn cdf(&self, threshold: f64) -> f64 {
        match self {
            Self::Normal(p) => p.cdf(threshold),
            Self::Uniform(p) => p.cdf(threshold),
            Self::TruncatedNormal(p) => p.cdf(threshold),
        }
    }
}
