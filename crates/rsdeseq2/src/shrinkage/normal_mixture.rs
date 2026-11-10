use crate::{DeseqError, errors::invalid_dimensions, probability::normal_upper_tail};

/// Fixed finite normal-mixture prior, including point masses (`sd == 0`).
///
/// Nonnegative component weights need not sum to one: each posterior normalizes
/// their likelihood-weighted values. Prior fitting and missing-observation
/// presentation belong to the higher-level shrinkage engine.
#[derive(Clone, Debug)]
pub struct NormalMixture {
    log_weights: Vec<f64>,
    means: Vec<f64>,
    variances: Vec<f64>,
}

/// A fixed normal-mixture prior conditioned on one finite normal observation.
#[derive(Clone, Debug)]
pub struct NormalMixturePosterior {
    weights: Vec<f64>,
    means: Vec<f64>,
    standard_deviations: Vec<f64>,
    mean: f64,
    standard_deviation: f64,
}

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("normal-mixture posterior: {}", reason.into()),
    }
}

impl NormalMixture {
    /// Validate a fixed mixture. Zero-weight components retain their positions.
    pub fn new(weights: Vec<f64>, means: Vec<f64>, sds: Vec<f64>) -> Result<Self, DeseqError> {
        if weights.is_empty() {
            return Err(invalid("at least one component is required"));
        }
        for (name, actual) in [
            ("component means", means.len()),
            ("component SDs", sds.len()),
        ] {
            if actual != weights.len() {
                return Err(invalid_dimensions(name, weights.len(), actual));
            }
        }
        if weights.iter().any(|w| !w.is_finite() || *w < 0.0) || !weights.iter().any(|w| *w > 0.0) {
            return Err(invalid(
                "weights must be finite, nonnegative and not all zero",
            ));
        }
        if means.iter().any(|m| !m.is_finite()) {
            return Err(invalid("component means must be finite"));
        }
        if sds
            .iter()
            .any(|s| !s.is_finite() || *s < 0.0 || !(s * s).is_finite())
        {
            return Err(invalid(
                "component SDs and variances must be finite and nonnegative",
            ));
        }
        Ok(Self {
            log_weights: weights.into_iter().map(f64::ln).collect(),
            means,
            variances: sds.into_iter().map(|s| s * s).collect(),
        })
    }

    /// Condition on an estimate with a strictly positive finite standard error.
    ///
    /// Likelihood weights are normalized in log space to retain far-tail
    /// observations and very small prior weights. No R runtime is used.
    pub fn posterior(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<NormalMixturePosterior, DeseqError> {
        let variance = standard_error * standard_error;
        if !estimate.is_finite()
            || !standard_error.is_finite()
            || standard_error <= 0.0
            || !variance.is_finite()
            || variance == 0.0
        {
            return Err(invalid(
                "estimate and positive standard error/variance must be finite",
            ));
        }
        let mut weights = Vec::with_capacity(self.means.len());
        let mut means = Vec::with_capacity(self.means.len());
        let mut standard_deviations = Vec::with_capacity(self.means.len());
        for ((&mean, &prior_variance), &log_weight) in self
            .means
            .iter()
            .zip(&self.variances)
            .zip(&self.log_weights)
        {
            let marginal_variance = variance + prior_variance;
            if !marginal_variance.is_finite() {
                return Err(invalid("marginal variance overflow"));
            }
            let marginal_sd = marginal_variance.sqrt();
            let z = (estimate - mean) / marginal_sd;
            let log_density =
                -0.5 * z * z - 0.5 * (2.0 * std::f64::consts::PI).ln() - marginal_sd.ln();
            weights.push(log_weight + log_density);
            means.push((variance * mean + estimate * prior_variance) / marginal_variance);
            standard_deviations.push((variance * prior_variance / marginal_variance).sqrt());
        }
        let max = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !max.is_finite() {
            return Err(invalid(
                "observation exceeds the representable log-likelihood range",
            ));
        }
        for weight in &mut weights {
            *weight = (*weight - max).exp();
        }
        let total: f64 = weights.iter().sum();
        for weight in &mut weights {
            *weight /= total;
        }
        let mean = weights.iter().zip(&means).map(|(w, m)| w * m).sum::<f64>();
        let second_moment = weights
            .iter()
            .zip(&means)
            .zip(&standard_deviations)
            .map(|((w, m), s)| w * (s * s + m * m))
            .sum::<f64>();
        if !mean.is_finite()
            || !second_moment.is_finite()
            || standard_deviations.iter().any(|s| !s.is_finite())
        {
            return Err(invalid("posterior moments exceed the representable range"));
        }
        Ok(NormalMixturePosterior {
            weights,
            means,
            standard_deviations,
            mean,
            standard_deviation: (second_moment.max(0.0) - mean * mean).max(0.0).sqrt(),
        })
    }
}

impl NormalMixturePosterior {
    /// Normalized posterior component probabilities, in prior input order.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    /// Posterior expectation on the observation's original scale.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Posterior SD from the mixture's first and second moments.
    pub fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }
    /// Right-continuous posterior CDF. Point masses contribute at their location.
    pub fn cdf(&self, threshold: f64) -> f64 {
        self.probability(threshold, true)
    }
    /// Posterior probability strictly above a threshold, evaluated directly
    /// without subtracting a near-one CDF. This retains small upper tails.
    pub fn survival(&self, threshold: f64) -> f64 {
        self.probability(threshold, false)
    }
    fn probability(&self, threshold: f64, lower: bool) -> f64 {
        if threshold.is_nan() {
            return f64::NAN;
        }
        self.weights
            .iter()
            .zip(&self.means)
            .zip(&self.standard_deviations)
            .map(|((w, mean), sd)| {
                let cdf = if *sd == 0.0 {
                    f64::from(if lower {
                        threshold >= *mean
                    } else {
                        threshold < *mean
                    })
                } else if lower {
                    normal_upper_tail((*mean - threshold) / sd)
                } else {
                    normal_upper_tail((threshold - *mean) / sd)
                };
                w * cdf
            })
            .sum()
    }
}
