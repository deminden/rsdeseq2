use crate::errors::scoped_invalid_options;
use crate::numeric::{NumericProfile, math::ProfileMath, numeric_profile};
const ERROR_SCOPE: &str = "uniform-mixture posterior";
use super::{LocationLikelihood, student_t, truncated_standard, uniform_interval::NarrowInterval};
use crate::math::elementary::language_exp;
use crate::{Error, errors::invalid_dimensions};

/// Fixed finite uniform-mixture prior, including zero-width point masses.
///
/// Host R profiles follow the original ashr likelihood and moment rules.
/// The portable profile evaluates Gaussian and Student-t interval
/// probabilities in the outer tail and integrates narrow intervals locally.
/// Numerically undefined observations retain NaN summaries.
#[derive(Clone, Debug)]
pub struct UniformMixture {
    log_weights: Vec<f64>,
    cdf_reversed: Option<bool>,
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
fn observation(x: f64, s: f64) -> Result<(), Error> {
    if !x.is_finite() || !s.is_finite() || s <= 0.0 {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "estimate and positive standard error must be finite",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum PortableLocation {
    Normal,
    StudentT(f64),
}

impl PortableLocation {
    fn for_likelihood(likelihood: &LocationLikelihood) -> Option<Self> {
        if numeric_profile() != NumericProfile::Portable {
            return None;
        }
        match likelihood {
            LocationLikelihood::Normal => Some(Self::Normal),
            LocationLikelihood::StudentT(df) => Some(Self::StudentT(*df)),
            _ => None,
        }
    }

    fn narrow_interval(self, a: f64, b: f64, estimate: f64, se: f64) -> Option<NarrowInterval> {
        match self {
            Self::Normal => truncated_standard::narrow_interval(a, b, estimate, se),
            Self::StudentT(df) => student_t::narrow_interval(a, b, estimate, se, df),
        }
    }

    fn log_average_density(self, a: f64, b: f64, estimate: f64, se: f64) -> f64 {
        if let Some(interval) = self.narrow_interval(a, b, estimate, se) {
            return interval.log_density(se);
        }
        self.log_interval_mass(a, b, estimate, se) - (b - a).ln_p()
    }

    fn log_interval_mass(self, a: f64, b: f64, estimate: f64, se: f64) -> f64 {
        let alpha = (a - estimate) / se;
        let beta = (b - estimate) / se;
        match self {
            Self::Normal => truncated_standard::log_interval_mass(alpha, beta),
            Self::StudentT(df) => student_t::log_interval_mass(alpha, beta, df),
        }
    }

    fn interval_cdf(self, a: f64, t: f64, b: f64, estimate: f64, se: f64) -> f64 {
        if t == a {
            return 0.0;
        }
        if matches!(self, Self::Normal)
            && let Some(cdf) = truncated_standard::far_tail_cdf(a, t, b, estimate, se)
        {
            return cdf;
        }
        // Ordinary interval masses already include the width. For narrow
        // intervals, defer the width correction until after subtracting the
        // log densities so it survives even when the densities are tiny.
        let log_mass_terms = |end| match self.narrow_interval(a, end, estimate, se) {
            Some(interval) => (interval.log_density(se), (end - a).ln_p()),
            None => (self.log_interval_mass(a, end, estimate, se), 0.0),
        };
        let (partial, partial_width) = log_mass_terms(t);
        let (full, full_width) = log_mass_terms(b);
        ((partial - full) + (partial_width - full_width))
            .exp_p()
            .clamp(0.0, 1.0)
    }

    fn narrow_moments(self, a: f64, b: f64, estimate: f64, se: f64) -> Option<(f64, f64)> {
        if matches!(self, Self::StudentT(df) if df <= 1.0) {
            return None;
        }
        let (mean, second) = self.narrow_interval(a, b, estimate, se)?.moments(a, b);
        let second = if matches!(self, Self::StudentT(df) if df <= 2.0) {
            f64::NAN
        } else {
            second
        };
        Some((mean, second))
    }
}

impl UniformMixture {
    /// Validate finite ordered endpoints and nonnegative, nonzero total weight.
    /// Weights need not sum to one; zero-weight positions are retained.
    pub fn new(weights: Vec<f64>, lower: Vec<f64>, upper: Vec<f64>) -> Result<Self, Error> {
        if weights.is_empty() {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "at least one component is required",
            ));
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
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "weights must be finite, nonnegative and not all zero",
            ));
        }
        if lower
            .iter()
            .zip(&upper)
            .any(|(&a, &b)| !a.is_finite() || !b.is_finite() || a > b || !(b - a).is_finite())
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "finite ordered endpoints and representable widths are required",
            ));
        }
        Ok(Self {
            likelihood: LocationLikelihood::Normal,
            cdf_reversed: None,
            log_weights: weights.into_iter().map(f64::ln_p).collect(),
            lower,
            upper,
        })
    }
    /// Select a Student-t observation likelihood. Infinite degrees of freedom
    /// retains the normal likelihood; finite values must be strictly positive.
    /// Low-df moment undefinedness follows the original component evaluator.
    pub fn with_student_t(self, degrees_of_freedom: f64) -> Result<Self, Error> {
        self.with_likelihood(LocationLikelihood::StudentT(degrees_of_freedom))
    }
    /// Select a validated standalone location-scale likelihood.
    pub fn with_likelihood(mut self, likelihood: LocationLikelihood) -> Result<Self, Error> {
        self.likelihood = likelihood.validated()?;
        Ok(self)
    }
    pub(super) fn has_count_likelihood(&self) -> bool {
        matches!(&self.likelihood,LocationLikelihood::Builtin(b) if b.is_count())
    }
    pub(super) fn orient_for(
        &mut self,
        x: &[f64],
        se: &[f64],
        rows: impl Iterator<Item = usize> + Clone,
    ) -> Result<(), Error> {
        if !self.has_count_likelihood() {
            return Ok(());
        }
        let mut total = 0.0;
        // Source matrices are column-major: accumulate observations within each component.
        for (&a, &b) in self.lower.iter().zip(&self.upper) {
            for i in rows.clone() {
                let l = self.likelihood.at(i, x.len())?;
                let d = l.log_cdf((x[i] - a) / se[i]) - l.log_cdf((x[i] - b) / se[i]);
                if !d.is_nan() {
                    total += d;
                }
            }
        }
        if total.is_nan() {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "count-likelihood CDF orientation is undefined",
            ));
        }
        self.cdf_reversed = Some(total < 0.0);
        Ok(())
    }
    /// Component log convolution densities in prior input order.
    pub fn log_likelihoods(&self, estimate: f64, standard_error: f64) -> Result<Vec<f64>, Error> {
        self.log_likelihoods_with(estimate, standard_error, &self.likelihood)
    }
    pub(super) fn log_likelihoods_at(
        &self,
        x: f64,
        s: f64,
        row: usize,
        n: usize,
    ) -> Result<Vec<f64>, Error> {
        self.log_likelihoods_with(x, s, self.likelihood.at(row, n)?.as_ref())
    }
    fn log_likelihoods_with(
        &self,
        estimate: f64,
        standard_error: f64,
        likelihood: &LocationLikelihood,
    ) -> Result<Vec<f64>, Error> {
        observation(estimate, standard_error)?;
        let portable = PortableLocation::for_likelihood(likelihood);
        let reverse = if matches!(likelihood,LocationLikelihood::Builtin(b) if b.is_count()) {
            match self.cdf_reversed {
                Some(reverse) => reverse,
                None => {
                    let total = self
                        .lower
                        .iter()
                        .zip(&self.upper)
                        .map(|(&a, &b)| {
                            likelihood.log_cdf((estimate - a) / standard_error)
                                - likelihood.log_cdf((estimate - b) / standard_error)
                        })
                        .filter(|d| !d.is_nan())
                        .sum::<f64>();
                    if total.is_nan() {
                        return Err(scoped_invalid_options(
                            ERROR_SCOPE,
                            "count-likelihood CDF orientation is undefined",
                        ));
                    }
                    total < 0.0
                }
            }
        } else {
            false
        };
        Ok(self
            .lower
            .iter()
            .zip(&self.upper)
            .map(|(&a, &b)| {
                if a == b {
                    let q = (estimate - a) / standard_error;
                    let density = likelihood.log_density(q);
                    return density - standard_error.ln_p();
                }
                if let Some(portable) = portable {
                    return portable.log_average_density(a, b, estimate, standard_error);
                }
                let left = likelihood.log_cdf((estimate - a) / standard_error);
                let right = likelihood.log_cdf((estimate - b) / standard_error);
                let (left, right) = if reverse {
                    (right, left)
                } else {
                    (left, right)
                };
                // Preserve the established orientation and ordinary subtraction semantics.
                left + (1.0 - language_exp(right - left)).ln_p() - (b - a).ln_p()
            })
            .collect())
    }
    fn posterior_weights(
        &self,
        estimate: f64,
        standard_error: f64,
        likelihood: &LocationLikelihood,
    ) -> Result<Vec<f64>, Error> {
        let likelihoods = self.log_likelihoods_with(estimate, standard_error, likelihood)?;
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
    ) -> Result<Vec<f64>, Error> {
        self.cdf_values_with(estimate, standard_error, thresholds, &self.likelihood)
    }
    pub(super) fn cdf_values_at(
        &self,
        x: f64,
        s: f64,
        c: &[f64],
        row: usize,
        n: usize,
    ) -> Result<Vec<f64>, Error> {
        self.cdf_values_with(x, s, c, self.likelihood.at(row, n)?.as_ref())
    }
    fn cdf_values_with(
        &self,
        estimate: f64,
        standard_error: f64,
        thresholds: &[f64],
        likelihood: &LocationLikelihood,
    ) -> Result<Vec<f64>, Error> {
        let weights = self.posterior_weights(estimate, standard_error, likelihood)?;
        let view = CdfView {
            lower: &self.lower,
            upper: &self.upper,
            weights: &weights,
            estimate,
            standard_error,
            likelihood,
        };
        Ok(thresholds.iter().map(|&c| view.cdf(c)).collect())
    }
    /// Condition on a finite estimate and positive finite standard error.
    /// Undefined endpoint differences propagate to NaN posterior summaries.
    pub fn posterior(
        &self,
        estimate: f64,
        standard_error: f64,
    ) -> Result<UniformMixturePosterior, Error> {
        self.posterior_with(estimate, standard_error, &self.likelihood)
    }
    pub(super) fn posterior_at(
        &self,
        x: f64,
        s: f64,
        row: usize,
        n: usize,
    ) -> Result<UniformMixturePosterior, Error> {
        self.posterior_with(x, s, self.likelihood.at(row, n)?.as_ref())
    }
    fn posterior_with(
        &self,
        estimate: f64,
        standard_error: f64,
        likelihood: &LocationLikelihood,
    ) -> Result<UniformMixturePosterior, Error> {
        let weights = self.posterior_weights(estimate, standard_error, likelihood)?;
        if weights.iter().any(|v| v.is_nan()) {
            return Ok(UniformMixturePosterior {
                weights,
                lower: self.lower.clone(),
                upper: self.upper.clone(),
                estimate,
                standard_error,
                likelihood: likelihood.clone(),
                mean: f64::NAN,
                standard_deviation: f64::NAN,
            });
        }
        let (mean, variance) = if matches!(
            PortableLocation::for_likelihood(likelihood),
            Some(PortableLocation::Normal)
        ) {
            self.portable_gaussian_moments(&weights, estimate, standard_error)
        } else {
            let (mean, second) = self.raw_moments(&weights, estimate, standard_error, likelihood);
            let second = if second < 0.0 { 0.0 } else { second };
            (mean, second - mean * mean)
        };
        Ok(UniformMixturePosterior {
            weights,
            lower: self.lower.clone(),
            upper: self.upper.clone(),
            estimate,
            standard_error,
            likelihood: likelihood.clone(),
            mean,
            standard_deviation: if variance < 0.0 { 0.0 } else { variance.sqrt() },
        })
    }

    fn raw_moments(
        &self,
        weights: &[f64],
        estimate: f64,
        standard_error: f64,
        likelihood: &LocationLikelihood,
    ) -> (f64, f64) {
        let mut mean = 0.0;
        let mut second = 0.0;
        let portable = PortableLocation::for_likelihood(likelihood);
        for ((&a, &b), &weight) in self.lower.iter().zip(&self.upper).zip(weights) {
            let (component_mean, component_second) = if a == b
                && matches!(
                    likelihood,
                    LocationLikelihood::Normal | LocationLikelihood::StudentT(_)
                ) {
                (a, a * a)
            } else if let Some(moments) = portable
                .and_then(|location| location.narrow_moments(a, b, estimate, standard_error))
            {
                moments
            } else {
                let alpha = (a - estimate) / standard_error;
                let beta = (b - estimate) / standard_error;
                let m = likelihood.moment(alpha, beta, 1);
                let m2 = likelihood.moment(alpha, beta, 2);
                (
                    if matches!(likelihood, LocationLikelihood::Normal) {
                        estimate + standard_error * m
                    } else {
                        estimate
                            - standard_error
                                * likelihood.moment(
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
        (mean, second)
    }

    fn portable_gaussian_moments(
        &self,
        weights: &[f64],
        estimate: f64,
        standard_error: f64,
    ) -> (f64, f64) {
        let mut total = 0.0;
        let mut mean = 0.0;
        let mut variance_sum = 0.0;
        for ((&a, &b), &weight) in self.lower.iter().zip(&self.upper).zip(weights) {
            if weight == 0.0 {
                continue;
            }
            let (component_mean, variance) =
                truncated_standard::portable_mean_variance(a, b, estimate, standard_error);
            variance_sum += weight * variance;
            if total == 0.0 {
                mean = component_mean;
                total = weight;
                continue;
            }
            // Combine within-component variance with the between-component
            // separation around the current mean.
            let next_total = total + weight;
            let fraction = weight / next_total;
            let spread_scale = (total * fraction).sqrt();
            let difference = component_mean - mean;
            let spread = if difference.is_finite() {
                mean += fraction * difference;
                spread_scale * difference
            } else {
                let separation = spread_scale * component_mean - spread_scale * mean;
                mean = (total / next_total) * mean + fraction * component_mean;
                separation
            };
            variance_sum += spread * spread;
            total = next_total;
        }
        (mean, variance_sum / total)
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
        let portable = PortableLocation::for_likelihood(self.likelihood);
        self.lower
            .iter()
            .zip(self.upper)
            .zip(self.weights)
            .map(|((&a, &b), &w)| {
                let probability = if threshold < a {
                    0.0
                } else if threshold >= b {
                    1.0
                } else if let Some(portable) = portable {
                    portable.interval_cdf(a, threshold, b, self.estimate, self.standard_error)
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
