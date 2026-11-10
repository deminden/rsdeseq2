use super::LocationLikelihood;
use crate::math::elementary::language_exp;
use crate::{DeseqError, errors::invalid_dimensions};

/// Fixed finite uniform-mixture prior, including zero-width point masses.
///
/// The finite-arithmetic likelihood and moment rules preserve original ashr
/// compatibility. Numerically undefined observations retain NaN summaries.
#[derive(Clone, Debug)]
pub struct UniformMixture {
    log_weights: Vec<f64>,
    likelihood: LocationLikelihood,
    lower: Vec<f64>,
    upper: Vec<f64>,
}
/// A uniform-mixture prior conditioned on one finite location-scale observation.
#[derive(Clone, Debug)]
pub struct UniformMixturePosterior {
    weights: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    estimate: f64,
    standard_error: f64,
    likelihood: LocationLikelihood,
    mean: f64,
    standard_deviation: f64,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("uniform-mixture posterior: {reason}"),
    }
}
fn observation(x: f64, s: f64) -> Result<(), DeseqError> {
    if !x.is_finite() || !s.is_finite() || s <= 0.0 {
        return Err(invalid(
            "estimate and positive standard error must be finite",
        ));
    }
    Ok(())
}
impl UniformMixture {
    /// Validate finite ordered endpoints and nonnegative, nonzero total weight.
    /// Weights need not sum to one; zero-weight positions are retained.
    pub fn new(weights: Vec<f64>, lower: Vec<f64>, upper: Vec<f64>) -> Result<Self, DeseqError> {
        if weights.is_empty() {
            return Err(invalid("at least one component is required"));
        }
        for (name, len) in [
            ("uniform lower endpoints", lower.len()),
            ("uniform upper endpoints", upper.len()),
        ] {
            if len != weights.len() {
                return Err(invalid_dimensions(name, weights.len(), len));
            }
        }
        if weights.iter().any(|v| !v.is_finite() || *v < 0.0) || !weights.iter().any(|v| *v > 0.0) {
            return Err(invalid(
                "weights must be finite, nonnegative and not all zero",
            ));
        }
        if lower
            .iter()
            .zip(&upper)
            .any(|(&a, &b)| !a.is_finite() || !b.is_finite() || a > b || !(b - a).is_finite())
        {
            return Err(invalid(
                "finite ordered endpoints and representable widths are required",
            ));
        }
        Ok(Self {
            likelihood: LocationLikelihood::Normal,
            log_weights: weights.into_iter().map(f64::ln).collect(),
            lower,
            upper,
        })
    }
    /// Select a Student-t observation likelihood. Infinite degrees of freedom
    /// retains the normal likelihood; finite values must be strictly positive.
    /// Low-df moment undefinedness follows the original component evaluator.
    pub fn with_student_t(self, degrees_of_freedom: f64) -> Result<Self, DeseqError> {
        self.with_likelihood(LocationLikelihood::StudentT(degrees_of_freedom))
    }
    /// Select a validated standalone location-scale likelihood.
    pub fn with_likelihood(mut self, likelihood: LocationLikelihood) -> Result<Self, DeseqError> {
        self.likelihood = likelihood.validated()?;
        Ok(self)
    }
    fn log_cdf(&self, q: f64) -> f64 {
        self.likelihood.log_cdf(q)
    }
    /// Component log convolution densities in prior input order.
    pub fn log_likelihoods(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<Vec<f64>, DeseqError> {
        observation(estimate, standard_error)?;
        Ok(self
            .lower
            .iter()
            .zip(&self.upper)
            .map(|(&a, &b)| {
                if a == b {
                    let q = (estimate - a) / standard_error;
                    let density = self.likelihood.log_density(q);
                    return density - standard_error.ln();
                }
                let left = self.log_cdf((estimate - a) / standard_error);
                let right = self.log_cdf((estimate - b) / standard_error);
                // Preserve the established orientation and ordinary subtraction semantics.
                left + (1.0 - language_exp(right - left)).ln() - (b - a).ln()
            })
            .collect())
    }
    fn posterior_weights(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<Vec<f64>, DeseqError> {
        let likelihoods = self.log_likelihoods(estimate, standard_error)?;
        let mut weights = likelihoods
            .iter()
            .zip(&self.log_weights)
            .map(|(l, w)| l + w)
            .collect::<Vec<_>>();
        let maximum = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if weights.iter().any(|v| v.is_nan()) || !maximum.is_finite() {
            return Ok(vec![f64::NAN; weights.len()]);
        }
        for w in &mut weights {
            *w = language_exp(*w - maximum);
        }
        let total: f64 = weights.iter().sum();
        for w in &mut weights {
            *w /= total;
        }
        Ok(weights)
    }
    /// Evaluate CDF cutoffs without evaluating unused moment expressions.
    pub fn cdf_values(
        &self,
        estimate: f64,
        standard_error: f64,
        thresholds: &[f64],
    ) -> Result<Vec<f64>, DeseqError> {
        let weights = self.posterior_weights(estimate, standard_error)?;
        let view = CdfView {
            lower: &self.lower,
            upper: &self.upper,
            weights: &weights,
            estimate,
            standard_error,
            likelihood: &self.likelihood,
        };
        Ok(thresholds.iter().map(|&c| view.cdf(c)).collect())
    }
    /// Condition on a finite estimate and positive finite standard error.
    /// Undefined endpoint differences propagate to NaN posterior summaries.
    pub fn posterior(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<UniformMixturePosterior, DeseqError> {
        let weights = self.posterior_weights(estimate, standard_error)?;
        if weights.iter().any(|v| v.is_nan()) {
            return Ok(UniformMixturePosterior {
                weights,
                lower: self.lower.clone(),
                upper: self.upper.clone(),
                estimate,
                standard_error,
                likelihood: self.likelihood.clone(),
                mean: f64::NAN,
                standard_deviation: f64::NAN,
            });
        }
        let mut mean = 0.0;
        let mut second = 0.0;
        for ((&a, &b), &weight) in self.lower.iter().zip(&self.upper).zip(&weights) {
            let (component_mean, component_second) =
                if a == b && !matches!(self.likelihood, LocationLikelihood::Restricted(_)) {
                    (a, a * a)
                } else {
                    let alpha = (a - estimate) / standard_error;
                    let beta = (b - estimate) / standard_error;
                    let m = self.likelihood.moment(alpha, beta, 1);
                    let m2 = self.likelihood.moment(alpha, beta, 2);
                    (
                        if matches!(self.likelihood, LocationLikelihood::Normal) {
                            estimate + standard_error * m
                        } else {
                            estimate
                                - standard_error
                                    * self.likelihood.moment(
                                        (estimate - b) / standard_error,
                                        (estimate - a) / standard_error,
                                        1,
                                    )
                        },
                        estimate * estimate
                            + 2.0 * estimate * standard_error * m
                            + standard_error * standard_error * m2,
                    )
                };
            mean += weight * component_mean;
            second += weight * component_second;
        }
        Ok(UniformMixturePosterior {
            weights,
            lower: self.lower.clone(),
            upper: self.upper.clone(),
            estimate,
            standard_error,
            likelihood: self.likelihood.clone(),
            mean,
            standard_deviation: {
                let second = if second < 0.0 { 0.0 } else { second };
                let variance = second - mean * mean;
                (if variance < 0.0 { 0.0 } else { variance }).sqrt()
            },
        })
    }
}
impl UniformMixturePosterior {
    /// Normalized posterior component probabilities in prior order.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    /// Posterior expectation on the original observation scale.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Posterior SD calculated from the first and second moments.
    pub fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }
    /// Right-continuous posterior CDF, including point-mass boundaries.
    pub fn cdf(&self, threshold: f64) -> f64 {
        CdfView {
            lower: &self.lower,
            upper: &self.upper,
            weights: &self.weights,
            estimate: self.estimate,
            standard_error: self.standard_error,
            likelihood: &self.likelihood,
        }
        .cdf(threshold)
    }
}
struct CdfView<'a> {
    lower: &'a [f64],
    upper: &'a [f64],
    weights: &'a [f64],
    estimate: f64,
    standard_error: f64,
    likelihood: &'a LocationLikelihood,
}
impl CdfView<'_> {
    fn log_cdf(&self, x: f64) -> f64 {
        self.likelihood.log_cdf(x)
    }
    pub fn cdf(&self, threshold: f64) -> f64 {
        if threshold.is_nan() || self.weights.iter().any(|v| v.is_nan()) {
            return f64::NAN;
        }
        self.lower
            .iter()
            .zip(self.upper)
            .zip(self.weights)
            .map(|((&a, &b), &w)| {
                let probability = if threshold < a {
                    0.0
                } else if threshold >= b {
                    1.0
                } else {
                    let left = self.log_cdf((self.estimate - a) / self.standard_error);
                    let middle = self.log_cdf((self.estimate - threshold) / self.standard_error);
                    let right = self.log_cdf((self.estimate - b) / self.standard_error);
                    let value =
                        (language_exp(middle - left) - 1.0) / (language_exp(right - left) - 1.0);
                    if value.is_nan() {
                        (threshold - a) / (b - a)
                    } else {
                        value
                    }
                };
                w * probability
            })
            .sum()
    }
}
