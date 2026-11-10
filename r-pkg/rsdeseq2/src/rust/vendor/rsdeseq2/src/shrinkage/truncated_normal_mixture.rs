//! Finite truncated-normal mixtures for the optional shrinkage migration.
use super::truncated_standard::{log_phi, scaled_mean, scaled_sd};
use crate::probability::normal_upper_tail;
use crate::{DeseqError, errors::invalid_dimensions, math::elementary::language_exp};

/// Fixed truncated-normal components, including half normals and point masses.
///
/// This compatibility primitive preserves the established ashr component
/// likelihood, including its zero-centered convolution for supplied components.
/// It does not yet replace the public optional shrinkage methods.
#[derive(Clone, Debug)]
pub struct TruncatedNormalMixture {
    log_weights: Vec<f64>,
    means: Vec<f64>,
    sds: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

/// One finite normal observation conditioned on a truncated-normal mixture.
#[derive(Clone, Debug)]
pub struct TruncatedNormalMixturePosterior {
    weights: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    means: Vec<f64>,
    sds: Vec<f64>,
    mean: f64,
    standard_deviation: f64,
}

fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("truncated-normal mixture: {reason}"),
    }
}
fn cdf(q: f64) -> f64 {
    normal_upper_tail(-q)
}
impl TruncatedNormalMixture {
    /// Validate nonnegative weights/SDs and ordered, possibly infinite bounds.
    /// Zero weights retain their component positions; weights need not sum to one.
    pub fn new(
        weights: Vec<f64>,
        means: Vec<f64>,
        sds: Vec<f64>,
        lower: Vec<f64>,
        upper: Vec<f64>,
    ) -> Result<Self, DeseqError> {
        if weights.is_empty() {
            return Err(invalid("at least one component is required"));
        }
        for (name, len) in [
            ("component means", means.len()),
            ("component SDs", sds.len()),
            ("lower bounds", lower.len()),
            ("upper bounds", upper.len()),
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
        if means.iter().any(|v| !v.is_finite())
            || sds
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0 || !(v * v).is_finite())
        {
            return Err(invalid(
                "finite means and finite nonnegative SDs/variances are required",
            ));
        }
        if lower.iter().zip(&upper).any(|(&a, &b)| {
            a.is_nan() || b.is_nan() || a > b || a == f64::INFINITY || b == f64::NEG_INFINITY
        }) {
            return Err(invalid(
                "ordered nonmissing bounds with nonempty finite support are required",
            ));
        }
        Ok(Self {
            log_weights: weights.into_iter().map(f64::ln).collect(),
            means,
            sds,
            lower,
            upper,
        })
    }

    /// Established component log convolution densities in prior input order.
    /// Undefined finite-arithmetic endpoint differences remain NaN or infinity.
    pub fn log_likelihoods(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<Vec<f64>, DeseqError> {
        let variance = standard_error * standard_error;
        if !estimate.is_finite()
            || !standard_error.is_finite()
            || standard_error <= 0.0
            || !variance.is_finite()
            || variance == 0.0
        {
            return Err(invalid(
                "finite estimate and positive finite SE/variance are required",
            ));
        }
        let mut result = Vec::with_capacity(self.sds.len());
        for ((&sd, &a), &b) in self.sds.iter().zip(&self.lower).zip(&self.upper) {
            let v = sd * sd + variance;
            let marginal_sd = v.sqrt();
            let density = log_phi(estimate / marginal_sd) - marginal_sd.ln();
            let value = if sd == 0.0 {
                density
            } else {
                let precision = (1.0 / (sd * sd) + 1.0 / variance).sqrt();
                let shift = ((1.0 / marginal_sd) * (sd / standard_error)) * estimate;
                if a == b {
                    // Match the established point-boundary limit, rather than
                    // replacing the supplied component by a different prior.
                    (log_phi(precision * b - shift).exp() / log_phi(b / sd).exp()).ln() + density
                        - v.ln()
                } else {
                    (cdf(precision * b - shift) - cdf(precision * a - shift)).ln() + density
                        - (cdf(b / sd) - cdf(a / sd)).ln()
                }
            };
            result.push(value);
        }
        Ok(result)
    }

    /// Condition on one observation. Undefined component arithmetic propagates
    /// to the posterior instead of silently discarding the affected component.
    pub fn posterior(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<TruncatedNormalMixturePosterior, DeseqError> {
        self.posterior_at(estimate, standard_error, 0, 1)
    }

    /// Condition a complete observation vector. The established component-SD
    /// routine recycles its bound vectors across a transposed observation matrix;
    /// retaining the batch shape is required to reproduce those summaries.
    pub fn posteriors(
        &self,
        estimates: &[f64],
        standard_errors: &[f64],
    ) -> Result<Vec<TruncatedNormalMixturePosterior>, DeseqError> {
        if estimates.len() != standard_errors.len() {
            return Err(invalid_dimensions(
                "posterior standard errors",
                estimates.len(),
                standard_errors.len(),
            ));
        }
        estimates
            .len()
            .checked_mul(self.sds.len())
            .ok_or_else(|| invalid("posterior matrix shape overflows"))?;
        estimates
            .iter()
            .zip(standard_errors)
            .enumerate()
            .map(|(i, (&x, &s))| self.posterior_at(x, s, i, estimates.len()))
            .collect()
    }

    pub(super) fn posterior_at(
        &self,
        estimate: f64,
        standard_error: f64,
        row: usize,
        observations: usize,
    ) -> Result<TruncatedNormalMixturePosterior, DeseqError> {
        let mut weights = self
            .log_likelihoods(estimate, standard_error)?
            .iter()
            .zip(&self.log_weights)
            .map(|(l, w)| l + w)
            .collect::<Vec<_>>();
        let max = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let undefined = weights.iter().any(|v| v.is_nan()) || !max.is_finite();
        if undefined {
            weights.fill(f64::NAN);
        } else {
            for w in &mut weights {
                *w = language_exp(*w - max);
            }
            let sum: f64 = weights.iter().sum();
            for w in &mut weights {
                *w /= sum;
            }
        }
        let mut means = Vec::with_capacity(weights.len());
        let mut sds = Vec::with_capacity(weights.len());
        let mut mean = 0.0;
        let mut second = 0.0;
        let v = standard_error * standard_error;
        for (component, ((((&prior_mean, &sd), &a), &b), &w)) in self
            .means
            .iter()
            .zip(&self.sds)
            .zip(&self.lower)
            .zip(&self.upper)
            .zip(&weights)
            .enumerate()
        {
            // Moment and CDF paths deliberately retain their respective source
            // arithmetic order: algebraic rearrangement changes narrow tails.
            let prior_fraction = 1.0 / ((sd * sd) / v + 1.0);
            let location = prior_fraction * prior_mean + (1.0 - prior_fraction) * estimate;
            let scale = (1.0 / (1.0 / (sd * sd) + 1.0 / v)).sqrt();
            let component_mean = scaled_mean(a, b, location, scale);
            let bounds = (component * observations + row) % self.sds.len();
            let component_sd = scaled_sd(self.lower[bounds], self.upper[bounds], location, scale);
            mean += w * component_mean;
            second += w * (component_sd * component_sd + component_mean * component_mean);
            let fraction = 1.0 / (v / (sd * sd) + 1.0);
            means.push(fraction * estimate + (1.0 - fraction) * prior_mean);
            sds.push(scale);
        }
        let second = if second < 0.0 { 0.0 } else { second };
        let var = second - mean * mean;
        Ok(TruncatedNormalMixturePosterior {
            weights,
            lower: self.lower.clone(),
            upper: self.upper.clone(),
            means,
            sds,
            mean,
            standard_deviation: (if var < 0.0 { 0.0 } else { var }).sqrt(),
        })
    }
}
impl TruncatedNormalMixturePosterior {
    /// Normalized component probabilities in prior order.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    /// Posterior expectation on the observation scale.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Posterior SD from component moments and mixture moments.
    pub fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }
    /// Right-continuous posterior CDF, preserving the source's zero fallback
    /// when a truncated component probability is numerically undefined.
    pub fn cdf(&self, threshold: f64) -> f64 {
        if threshold.is_nan() || self.weights.iter().any(|v| v.is_nan()) {
            return f64::NAN;
        }
        self.weights
            .iter()
            .zip(&self.lower)
            .zip(&self.upper)
            .zip(&self.means)
            .zip(&self.sds)
            .map(|((((&w, &a), &b), &mean), &sd)| {
                let mut p = if threshold < a {
                    0.0
                } else if threshold >= b {
                    1.0
                } else {
                    let left = cdf((a - mean) / sd);
                    (cdf((threshold - mean) / sd) - left) / (cdf((b - mean) / sd) - left)
                };
                if p.is_nan() {
                    p = 0.0;
                }
                w * p
            })
            .sum()
    }
}
