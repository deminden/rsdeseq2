//! Prior-only summaries and component selection for complete shrinkage batches.
use super::{
    MixtureComponents,
    truncated_standard::{scaled_mean, scaled_sd},
};
use crate::{DeseqError, probability::normal_upper_tail};

fn normal_cdf(q: f64, mean: f64, sd: f64) -> f64 {
    if q.is_nan() || mean.is_nan() {
        f64::NAN
    } else if sd == 0.0 {
        f64::from(q >= mean)
    } else {
        normal_upper_tail((mean - q) / sd)
    }
}
impl MixtureComponents {
    pub(super) fn select(&self, keep: &[usize]) -> Self {
        let select = |v: &[f64]| keep.iter().map(|&i| v[i]).collect::<Vec<_>>();
        match self {
            Self::Normal { means, sds } => Self::Normal {
                means: select(means),
                sds: select(sds),
            },
            Self::Uniform { lower, upper } => Self::Uniform {
                lower: select(lower),
                upper: select(upper),
            },
            Self::TruncatedNormal {
                means,
                sds,
                lower,
                upper,
            } => Self::TruncatedNormal {
                means: select(means),
                sds: select(sds),
                lower: select(lower),
                upper: select(upper),
            },
        }
    }
    pub(super) fn moments(&self, weights: &[f64]) -> (f64, f64, Vec<bool>) {
        let mut mean = 0.0;
        let mut second = 0.0;
        let mut zero = Vec::with_capacity(weights.len());
        for (i, &weight) in weights.iter().enumerate() {
            let (m, s) = match self {
                Self::Normal { means, sds } => (means[i], sds[i]),
                Self::Uniform { lower, upper } => (
                    (lower[i] + upper[i]) / 2.0,
                    (upper[i] - lower[i]) / 12.0_f64.sqrt(),
                ),
                Self::TruncatedNormal {
                    means,
                    sds,
                    lower,
                    upper,
                } => (
                    scaled_mean(lower[i], upper[i], means[i], sds[i]),
                    scaled_sd(lower[i], upper[i], means[i], sds[i]),
                ),
            };
            mean += weight * m;
            second += weight * (s * s + m * m);
            zero.push(m == 0.0 && s == 0.0);
        }
        (mean, (second - mean * mean).sqrt(), zero)
    }
    pub(super) fn prior_cdf(&self, weights: &[f64], threshold: f64) -> Result<f64, DeseqError> {
        // The established truncated-prior CDF fails when no component contains
        // the requested cutoff. Preserve that accepted/rejected decision here.
        if let Self::TruncatedNormal { lower, upper, .. } = self
            && !lower
                .iter()
                .zip(upper)
                .any(|(&a, &b)| a <= threshold && threshold <= b)
        {
            return Err(DeseqError::InvalidOptions {
                reason: "truncated-normal prior CDF requires a component containing the cutoff"
                    .into(),
            });
        }
        Ok(weights
            .iter()
            .enumerate()
            .map(|(i, &w)| w * self.component_cdf(i, threshold))
            .sum())
    }
    pub(super) fn component_cdf(&self, i: usize, threshold: f64) -> f64 {
        match self {
            Self::Normal { means, sds } => normal_cdf(threshold, means[i], sds[i]),
            Self::Uniform { lower, upper } => {
                if threshold.is_nan() {
                    f64::NAN
                } else if threshold >= upper[i] {
                    1.0
                } else if threshold <= lower[i] {
                    0.0
                } else {
                    (threshold - lower[i]) / (upper[i] - lower[i])
                }
            }
            Self::TruncatedNormal {
                means,
                sds,
                lower,
                upper,
            } => {
                if threshold < lower[i] {
                    0.0
                } else if threshold > upper[i] {
                    1.0
                } else {
                    let left = normal_cdf(lower[i], means[i], sds[i]);
                    (normal_cdf(threshold, means[i], sds[i]) - left)
                        / (normal_cdf(upper[i], means[i], sds[i]) - left)
                }
            }
        }
    }
}
