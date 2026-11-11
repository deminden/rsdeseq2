//! Independent queries of a fitted mixture, including original-input thresholds.
use super::{
    LocationLikelihood, MixtureComponents, MixtureGrid,
    location_missing::{NA, excluded, is_na},
};
use crate::DeseqError;
/// Raw posterior CDF configuration. Observation preprocessing is explicit so a
/// fitted alpha/t model can still answer the wrapper's original normal queries.
#[derive(Clone, Debug)]
pub struct LocationCdfOptions {
    pub likelihood: LocationLikelihood,
    pub alpha: f64,
    pub thresholds: Vec<f64>,
}
impl Default for LocationCdfOptions {
    fn default() -> Self {
        Self {
            likelihood: LocationLikelihood::Normal,
            alpha: 0.0,
            thresholds: vec![],
        }
    }
}
/// Evaluate observation-major raw posterior CDFs with the unpruned prior.
pub fn location_cdf_batch(
    grid: &MixtureGrid,
    weights: &[f64],
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationCdfOptions,
) -> Result<Vec<f64>, DeseqError> {
    let invalid = |s: &str| DeseqError::InvalidOptions {
        reason: format!("location CDF query: {s}"),
    };
    if standard_errors.len() != estimates.len() && standard_errors.len() != 1 {
        return Err(invalid("SEs must have length one or match estimates"));
    }
    if !options.alpha.is_finite() || options.thresholds.iter().any(|v| v.is_nan()) {
        return Err(invalid("finite alpha and nonmissing cutoffs are required"));
    }
    let n = estimates.len();
    options.likelihood.validate_observations(n)?;
    n.checked_mul(weights.len())
        .ok_or_else(|| invalid("posterior shape overflows"))?;
    let mut x = Vec::with_capacity(n);
    let mut se = Vec::with_capacity(n);
    for (i, &estimate) in estimates.iter().enumerate() {
        let s = standard_errors[if standard_errors.len() == 1 { 0 } else { i }];
        if s < 0.0 {
            return Err(invalid("SEs must be nonnegative"));
        }
        x.push(estimate / s.powf(options.alpha));
        se.push(s.powf(1.0 - options.alpha));
    }
    let prior = grid.prior_with_likelihood(weights, &options.likelihood)?;
    if let MixtureComponents::TruncatedNormal { lower, upper, .. } = &grid.components
        && options
            .thresholds
            .iter()
            .any(|&c| !lower.iter().zip(upper).any(|(&a, &b)| a <= c && c <= b))
    {
        return Err(invalid("a truncated component must contain each cutoff"));
    }
    let count = x
        .iter()
        .zip(&se)
        .filter(|(x, s)| x.is_nan() || s.is_nan())
        .count();
    let mut rank = 0;
    let mut cdf = Vec::with_capacity(
        n.checked_mul(options.thresholds.len())
            .ok_or_else(|| invalid("CDF output shape overflows"))?,
    );
    for i in 0..n {
        if x[i].is_finite() && se[i] > 0.0 && se[i].is_finite() {
            cdf.extend(prior.cdf_values_at(x[i], se[i], &options.thresholds, i, n)?);
        } else {
            cdf.extend(
                excluded(
                    &grid.components,
                    weights,
                    x[i],
                    se[i],
                    &options.thresholds,
                    rank,
                    count,
                )
                .0,
            );
            if x[i].is_nan() || se[i].is_nan() {
                rank += 1;
            }
        }
    }
    Ok(cdf)
}
/// Wrapper-compatible thresholded s-values. Missing local rates retain NA versus
/// NaN; signed zero and ties preserve original order. Rates are not clipped.
pub fn threshold_svalues(
    means: &[f64],
    cdf_positive: &[f64],
    cdf_negative: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    if means.len() != cdf_positive.len() || means.len() != cdf_negative.len() {
        return Err(DeseqError::InvalidOptions {
            reason: "threshold s-value vectors must have matching lengths".into(),
        });
    }
    let rates = means
        .iter()
        .zip(cdf_positive)
        .zip(cdf_negative)
        .map(|((&m, &p), &n)| {
            if m.is_nan() {
                NA
            } else if m > 0.0 {
                p
            } else {
                1.0 - n
            }
        })
        .collect::<Vec<_>>();
    let mut order = (0..rates.len()).collect::<Vec<_>>();
    order.sort_by(|&i, &j| match (rates[i].is_nan(), rates[j].is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => rates[i].partial_cmp(&rates[j]).unwrap(),
    });
    let mut cumulative = 0.0;
    let mut result = vec![0.0; rates.len()];
    for (rank, &i) in order.iter().enumerate() {
        if is_na(cumulative) || is_na(rates[i]) {
            cumulative = NA;
        } else {
            cumulative += rates[i];
        }
        result[i] = cumulative / (rank + 1) as f64;
    }
    Ok(result)
}
