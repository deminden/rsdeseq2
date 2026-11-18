//! Independent queries of a fitted mixture, including original-input thresholds.
use crate::errors::scoped_invalid_options;
use crate::numeric::math::ProfileMath;
const ERROR_SCOPE: &str = "location CDF query";
use super::{
    LocationLikelihood, MixtureComponents, MixtureGrid,
    location_missing::{NA, excluded, is_na},
};
use crate::Error;
/// Raw posterior CDF configuration. Observation preprocessing is explicit so a
/// fitted alpha/t model can still answer the wrapper's original normal queries.
#[derive(Clone, Debug)]
pub struct LocationCdfOptions {
    /// Observation likelihood.
    pub likelihood: LocationLikelihood,
    /// ashr's `alpha`.
    pub alpha: f64,
    /// CDF cutoffs.
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
) -> Result<Vec<f64>, Error> {
    if standard_errors.len() != estimates.len() && standard_errors.len() != 1 {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "SEs must have length one or match estimates",
        ));
    }
    if !options.alpha.is_finite() || options.thresholds.iter().any(|v| v.is_nan()) {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "finite alpha and nonmissing cutoffs are required",
        ));
    }
    let n = estimates.len();
    options.likelihood.validate_observations(n)?;
    n.checked_mul(weights.len())
        .ok_or_else(|| scoped_invalid_options(ERROR_SCOPE, "posterior shape overflows"))?;
    let mut x = Vec::with_capacity(n);
    let mut se = Vec::with_capacity(n);
    for (i, &estimate) in estimates.iter().enumerate() {
        let s = standard_errors[if standard_errors.len() == 1 { 0 } else { i }];
        if s < 0.0 {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "SEs must be nonnegative",
            ));
        }
        x.push(estimate / s.powf_p(options.alpha));
        se.push(s.powf_p(1.0 - options.alpha));
    }
    let prior = grid
        .prior_with_likelihood(weights, &options.likelihood)?
        .orient_for(&x, &se, 0..n)?;
    if let MixtureComponents::TruncatedNormal { lower, upper, .. } = &grid.components
        && options
            .thresholds
            .iter()
            .any(|&c| !lower.iter().zip(upper).any(|(&a, &b)| a <= c && c <= b))
    {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "a truncated component must contain each cutoff",
        ));
    }
    let count = x
        .iter()
        .zip(&se)
        .filter(|(x, s)| x.is_nan() || s.is_nan())
        .count();
    let mut rank = 0;
    let mut cdf = Vec::with_capacity(
        n.checked_mul(options.thresholds.len())
            .ok_or_else(|| scoped_invalid_options(ERROR_SCOPE, "CDF output shape overflows"))?,
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
) -> Result<Vec<f64>, Error> {
    if means.len() != cdf_positive.len() || means.len() != cdf_negative.len() {
        return Err(Error::InvalidOptions {
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
        _ => rates[i]
            .partial_cmp(&rates[j])
            .unwrap_or(std::cmp::Ordering::Equal),
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
