use super::intercept::intercept_value_error;
use super::*;
use crate::errors::{ErrorKind, Requirement, invalid_calculation, invalid_parameter};
use crate::matrix::try_alloc_matrix;
// Prior input checks, stable arithmetic and quantiles.

pub(super) fn validate_beta_prior_inputs(
    beta_matrix: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    coefficient_names: Option<&[String]>,
    options: BetaPriorVarianceOptions,
) -> Result<(), DeseqError> {
    if base_mean.len() != beta_matrix.n_rows() {
        return Err(invalid_dimensions(
            "beta prior variance baseMean rows",
            beta_matrix.n_rows(),
            base_mean.len(),
        ));
    }
    if disp_fit.len() != beta_matrix.n_rows() {
        return Err(invalid_dimensions(
            "beta prior variance dispFit rows",
            beta_matrix.n_rows(),
            disp_fit.len(),
        ));
    }
    if let Some(names) = coefficient_names
        && names.len() != beta_matrix.n_cols()
    {
        return Err(invalid_dimensions(
            "beta prior variance coefficient names",
            beta_matrix.n_cols(),
            names.len(),
        ));
    }
    validate_upper_quantile(options.upper_quantile)?;
    if !options.wide_prior_variance.is_finite() || options.wide_prior_variance <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Options,
            "wide beta prior variance",
            Requirement::FiniteAndPositive,
        ));
    }
    if !options.finite_beta_cutoff.is_finite() || options.finite_beta_cutoff <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Options,
            "finite beta cutoff",
            Requirement::FiniteAndPositive,
        ));
    }
    Ok(())
}

pub(super) fn beta_prior_weights(
    base_mean: &[f64],
    disp_fit: &[f64],
    n_rows: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut weights = Vec::with_capacity(n_rows);
    for row in 0..n_rows {
        validate_positive_finite(base_mean[row], "beta prior baseMean", row)?;
        validate_positive_finite(disp_fit[row], "beta prior dispFit", row)?;
        weights.push(checked_mean_dispersion_working_weight(
            base_mean[row],
            disp_fit[row],
            row,
            "beta prior row weight",
        )?);
    }
    Ok(weights)
}

pub(super) fn checked_mean_dispersion_working_weight(
    mean: f64,
    dispersion: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let inv_mean = checked::try_div(1.0, mean).nonfinite_at(context, Some(index))?;
    let denominator = checked::try_add(inv_mean, dispersion).nonfinite_at(context, Some(index))?;
    checked::try_div(1.0, denominator).nonfinite_at(context, Some(index))
}

pub(super) fn checked_intercept_beta_se(sigma: f64, gene: usize) -> Result<f64, DeseqError> {
    checked::try_mul(std::f64::consts::LOG2_E, sigma.sqrt()).map_err(|failure| {
        intercept_value_error(
            ErrorKind::Dispersion,
            gene,
            None,
            "intercept beta standard error",
            "non-finite",
            failure.value,
        )
    })
}

pub(super) fn checked_intercept_beta_covariance(
    sigma: f64,
    gene: usize,
) -> Result<f64, DeseqError> {
    const LOG2_E_SQUARED: f64 = std::f64::consts::LOG2_E * std::f64::consts::LOG2_E;
    checked::try_mul(LOG2_E_SQUARED, sigma).map_err(|failure| {
        intercept_value_error(
            ErrorKind::Dispersion,
            gene,
            None,
            "intercept beta covariance",
            "non-finite",
            failure.value,
        )
    })
}

pub(super) fn validate_upper_quantile(upper_quantile: f64) -> Result<(), DeseqError> {
    if !upper_quantile.is_finite() || upper_quantile <= 0.0 || upper_quantile >= 1.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "upper quantile must be finite and between 0 and 1".to_string(),
        });
    }
    Ok(())
}

pub(super) fn finite_abs_values(values: &[f64]) -> Result<Vec<f64>, DeseqError> {
    let out = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .map(f64::abs)
        .collect::<Vec<_>>();
    if out.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "beta prior variance quantile needs at least one finite beta".to_string(),
        });
    }
    Ok(out)
}

/// R `quantile(type = 7)` on an unsorted copy. Equal neighbours return the
/// neighbour unchanged (`stable_interpolate`). A failed interpolation retains
/// the rejected value. The VST and filtering versions always evaluate
/// `lo + f * (hi - lo)`.
pub(super) fn quantile_type7_stable(
    mut values: Vec<f64>,
    probability: f64,
) -> Result<f64, DeseqError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::InvalidOptions {
            reason: "quantile probability must be finite and between 0 and 1".to_string(),
        });
    }
    values.sort_by(|a, b| a.total_cmp(b));
    if values.len() == 1 {
        return Ok(values[0]);
    }
    let h = (values.len() as f64 - 1.0) * probability + 1.0;
    let lower = h.floor() as usize;
    let fraction = h - lower as f64;
    if lower == 0 {
        Ok(values[0])
    } else if lower >= values.len() {
        Ok(values[values.len() - 1])
    } else {
        stable_interpolate(values[lower - 1], values[lower], fraction).map_err(|failure| {
            invalid_calculation(
                ErrorKind::Options,
                "beta prior quantile",
                None,
                failure.value,
                "beta prior variance quantile produced non-finite quantile",
            )
        })
    }
}

pub(super) fn weighted_abs_quantile(
    betas: &[f64],
    weights: impl Iterator<Item = f64>,
    probability: f64,
) -> Result<f64, DeseqError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted quantile probability must be finite and between 0 and 1".to_string(),
        });
    }
    let mut pairs = try_alloc_matrix(betas.len(), 1, "weighted beta prior pairs")?;
    for (idx, (beta, weight)) in betas.iter().copied().zip(weights).enumerate() {
        if beta.is_finite() && weight.is_finite() && weight > 0.0 {
            pairs.push((beta.abs(), weight));
        } else if !weight.is_finite() || weight < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "beta prior variance weight".to_string(),
                index: Some(idx),
                value: weight,
            });
        }
    }
    if pairs.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite weights"
                .to_string(),
        });
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

    let weight_sum =
        checked::try_sum(pairs.iter().map(|(_, weight)| *weight)).map_err(|failure| {
            invalid_calculation(
                ErrorKind::Options,
                "weighted beta prior total weight",
                None,
                failure.value,
                "weighted beta prior variance quantile needs finite total weight",
            )
        })?;
    if weight_sum <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite total weight"
                .to_string(),
        });
    }

    let norm_scale = checked::try_div(pairs.len() as f64, weight_sum).map_err(|failure| {
        invalid_calculation(
            ErrorKind::Options,
            "weighted beta prior normalization scale",
            None,
            failure.value,
            "weighted beta prior variance quantile produced non-finite normalization scale",
        )
    })?;
    // Merge ties in stable input order while compacting the sorted pairs.
    let mut unique = 0;
    for index in 0..pairs.len() {
        let (value, weight) = pairs[index];
        let normalized_weight = checked::try_mul(weight, norm_scale).map_err(|failure| {
            invalid_calculation(
                ErrorKind::Options,
                "weighted beta prior normalized weight",
                None,
                failure.value,
                "weighted beta prior variance quantile produced non-finite normalized weight",
            )
        })?;
        if normalized_weight <= 0.0 {
            return Err(invalid_calculation(
                ErrorKind::Options,
                "weighted beta prior normalized weight",
                None,
                normalized_weight,
                "weighted beta prior variance quantile needs positive normalized weight",
            ));
        }
        if let Some((last_value, last_weight)) = pairs[..unique].last_mut()
            && *last_value == value
        {
            let Some(next_weight) = checked::add(*last_weight, normalized_weight) else {
                return Err(DeseqError::InvalidOptions {
                    reason:
                        "weighted beta prior variance quantile produced non-finite merged weight"
                            .to_string(),
                });
            };
            *last_weight = next_weight;
            continue;
        }
        pairs[unique] = (value, normalized_weight);
        unique += 1;
    }
    pairs.truncate(unique);
    if pairs.len() == 1 {
        return Ok(pairs[0].0);
    }

    let Some(n) = checked::sum(pairs.iter().map(|(_, weight)| *weight)) else {
        return Err(DeseqError::InvalidOptions {
            reason:
                "weighted beta prior variance quantile produced non-finite total normalized weight"
                    .to_string(),
        });
    };
    let order = 1.0 + (n - 1.0) * probability;
    if !order.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite order".to_string(),
        });
    }
    let low = order.floor().max(1.0);
    let high = (low + 1.0).min(n);
    let fraction = order.fract();

    let mut cumulative = 0.0;
    for (_, weight) in &mut pairs {
        let Some(next_cumulative) = checked::add(cumulative, *weight) else {
            return Err(DeseqError::InvalidOptions {
                reason:
                    "weighted beta prior variance quantile produced non-finite cumulative weight"
                        .to_string(),
            });
        };
        cumulative = next_cumulative;
        *weight = cumulative;
    }

    let order_statistic = |position| {
        let index = pairs.partition_point(|&(_, cumulative)| cumulative < position);
        pairs[index.min(pairs.len() - 1)].0
    };
    stable_interpolate(order_statistic(low), order_statistic(high), fraction).map_err(|failure| {
        invalid_calculation(
            ErrorKind::Options,
            "weighted beta prior quantile",
            None,
            failure.value,
            "weighted beta prior variance quantile produced non-finite quantile",
        )
    })
}

pub(super) fn stable_interpolate(left: f64, right: f64, fraction: f64) -> checked::Checked {
    if left == right {
        return if left.is_finite() {
            Ok(left)
        } else {
            checked::reject(left)
        };
    }
    let delta = checked::try_add(right, -left)?;
    let offset = checked::try_mul(fraction, delta)?;
    checked::try_add(left, offset)
}
