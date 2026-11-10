//! Excluded-row raw CDF/likelihood semantics, distinct from summary fallback.
use super::{MixtureComponents, truncated_standard::log_phi};
use crate::math::elementary::language_exp;

pub(super) const NA: f64 = f64::from_bits(0x7ff0_0000_0000_07a2);
pub(super) fn is_na(x: f64) -> bool {
    x.is_nan() && x.to_bits() as u32 == 1954
}
pub(super) fn log_marginal(likelihoods: &[f64], weights: &[f64]) -> f64 {
    if likelihoods.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let max = likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = likelihoods
        .iter()
        .zip(weights)
        .map(|(l, w)| w * language_exp(l - max))
        .sum();
    sum.ln() + max
}
pub(super) fn excluded(
    components: &MixtureComponents,
    weights: &[f64],
    x: f64,
    se: f64,
    thresholds: &[f64],
    missing_row: usize,
    missing_rows: usize,
) -> (Vec<f64>, f64) {
    if x.is_nan() || se.is_nan() {
        let missing = if is_na(x) || is_na(se) { NA } else { f64::NAN };
        // The original N-by-K posterior-probability matrix recycles the K
        // prior weights across its selected M missing rows, column-major.
        let indices = (0..weights.len())
            .map(|j| (j * missing_rows + missing_row) % weights.len())
            .collect::<Vec<_>>();
        let recycled = indices.iter().map(|&j| weights[j]).collect::<Vec<_>>();
        let cdf = thresholds
            .iter()
            .map(|&cutoff| match components {
                MixtureComponents::Normal { .. } => indices
                    .iter()
                    .map(|&j| weights[j] * components.component_cdf(j, cutoff))
                    .sum(),
                MixtureComponents::Uniform { lower, upper } => {
                    if is_na(missing)
                        && lower
                            .iter()
                            .zip(upper)
                            .any(|(&a, &b)| a <= cutoff && cutoff < b)
                    {
                        NA
                    } else {
                        components.prior_cdf(&recycled, cutoff).unwrap_or(f64::NAN)
                    }
                }
                MixtureComponents::TruncatedNormal { lower, upper, .. } => recycled
                    .iter()
                    .zip(lower)
                    .zip(upper)
                    .map(|((&w, &a), &b)| {
                        w * if cutoff < a {
                            0.0
                        } else if cutoff >= b {
                            1.0
                        } else if is_na(missing) {
                            NA
                        } else {
                            0.0
                        }
                    })
                    .sum(),
            })
            .collect();
        return (cdf, missing);
    }
    if se != 0.0 || !x.is_finite() {
        return (vec![f64::NAN; thresholds.len()], f64::NAN);
    }
    match components {
        MixtureComponents::Normal { means, sds } => {
            if sds.contains(&0.0) {
                return (vec![f64::NAN; thresholds.len()], f64::NAN);
            }
            let ll = means
                .iter()
                .zip(sds)
                .map(|(&m, &s)| log_phi((x - m) / s) - s.ln())
                .collect::<Vec<_>>();
            (
                thresholds
                    .iter()
                    .map(|&cutoff| f64::from(cutoff >= x))
                    .collect(),
                log_marginal(&ll, weights),
            )
        }
        MixtureComponents::Uniform { lower, upper } => {
            if lower
                .iter()
                .zip(upper)
                .any(|(&a, &b)| a == b || x == a || x == b)
            {
                return (vec![f64::NAN; thresholds.len()], f64::NAN);
            }
            let ll = lower
                .iter()
                .zip(upper)
                .map(|(&a, &b)| {
                    if a < x && x < b {
                        -(b - a).ln()
                    } else {
                        f64::NEG_INFINITY
                    }
                })
                .collect::<Vec<_>>();
            let log_likelihood = log_marginal(&ll, weights);
            if log_likelihood.is_nan() {
                return (vec![f64::NAN; thresholds.len()], log_likelihood);
            }
            let log_weights = ll
                .iter()
                .zip(weights)
                .map(|(l, w)| l + w.ln())
                .collect::<Vec<_>>();
            let max = log_weights
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let mut posterior = log_weights
                .iter()
                .map(|w| language_exp(w - max))
                .collect::<Vec<_>>();
            let total: f64 = posterior.iter().sum();
            for w in &mut posterior {
                *w /= total;
            }
            let cdf = thresholds
                .iter()
                .map(|&cutoff| {
                    posterior
                        .iter()
                        .zip(lower)
                        .zip(upper)
                        .map(|((&w, &a), &b)| {
                            let p = if cutoff < a {
                                0.0
                            } else if cutoff >= b {
                                1.0
                            } else if cutoff == x {
                                (cutoff - a) / (b - a)
                            } else {
                                f64::from(cutoff > x)
                            };
                            w * p
                        })
                        .sum()
                })
                .collect();
            (cdf, log_likelihood)
        }
        // Non-normal fixed-prior SE-zero compatibility is tracked separately
        // from ordinary excluded-row summaries; these expressions are undefined.
        MixtureComponents::TruncatedNormal { .. } => (vec![f64::NAN; thresholds.len()], f64::NAN),
    }
}
