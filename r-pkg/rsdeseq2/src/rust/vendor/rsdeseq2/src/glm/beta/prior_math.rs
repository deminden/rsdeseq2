// Prior input checks, stable arithmetic and quantiles.

fn validate_beta_prior_inputs(
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
        return Err(DeseqError::InvalidOptions {
            reason: "wide beta prior variance must be finite and positive".to_string(),
        });
    }
    if !options.finite_beta_cutoff.is_finite() || options.finite_beta_cutoff <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "finite beta cutoff must be finite and positive".to_string(),
        });
    }
    Ok(())
}

fn beta_prior_weights(
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

fn checked_mean_dispersion_working_weight(
    mean: f64,
    dispersion: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let Some(inv_mean) = checked_div2(1.0, mean) else {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: f64::NAN,
        });
    };
    let Some(denominator) = checked_sum2(inv_mean, dispersion) else {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: f64::NAN,
        });
    };
    checked_div2(1.0, denominator).ok_or_else(|| DeseqError::NonFiniteValue {
        context: context.to_string(),
        index: Some(index),
        value: f64::NAN,
    })
}

fn checked_mean(values: &[f64]) -> Option<f64> {
    checked_scaled_mean(values)
}

fn checked_scaled_sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    let values = values.into_iter().collect::<Vec<_>>();
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let normalized_sum = checked_sum(
        values
            .iter()
            .copied()
            .map(|value| checked_div2(value, scale))
            .collect::<Option<Vec<_>>>()?,
    )?;
    checked_product2(normalized_sum, scale)
}

fn checked_scaled_mean(values: &[f64]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let normalized_sum = checked_sum(
        values
            .iter()
            .copied()
            .map(|value| checked_div2(value, scale))
            .collect::<Option<Vec<_>>>()?,
    )?;
    checked_product2(checked_div2(normalized_sum, values.len() as f64)?, scale)
}

fn checked_sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    values.into_iter().try_fold(0.0, checked_sum2)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn checked_square(value: f64) -> Option<f64> {
    let square = value * value;
    (value.is_finite() && square.is_finite()).then_some(square)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let product = left * right;
    (left.is_finite() && right.is_finite() && product.is_finite()).then_some(product)
}

fn checked_div2(left: f64, right: f64) -> Option<f64> {
    let quotient = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite())
        .then_some(quotient)
}

fn checked_intercept_beta_se(sigma: f64, gene: usize) -> Result<f64, DeseqError> {
    let se = sigma.sqrt();
    let Some(value) = checked_product2(std::f64::consts::LOG2_E, se) else {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("gene {gene} has non-finite intercept beta standard error"),
        });
    };
    Ok(value)
}

fn checked_intercept_beta_covariance(sigma: f64, gene: usize) -> Result<f64, DeseqError> {
    let log2_e = std::f64::consts::LOG2_E;
    let Some(log2_e_squared) = checked_product2(log2_e, log2_e) else {
        return Err(DeseqError::InvalidDispersion {
            reason: "non-finite log2 covariance scaling factor".to_string(),
        });
    };
    let Some(value) = checked_product2(log2_e_squared, sigma) else {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("gene {gene} has non-finite intercept beta covariance"),
        });
    };
    Ok(value)
}

fn validate_upper_quantile(upper_quantile: f64) -> Result<(), DeseqError> {
    if !upper_quantile.is_finite() || upper_quantile <= 0.0 || upper_quantile >= 1.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "upper quantile must be finite and between 0 and 1".to_string(),
        });
    }
    Ok(())
}

fn finite_abs_values(values: &[f64]) -> Result<Vec<f64>, DeseqError> {
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

fn quantile_type7(mut values: Vec<f64>, probability: f64) -> Result<f64, DeseqError> {
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
        let quantile = stable_interpolate(values[lower - 1], values[lower], fraction);
        if quantile.is_finite() {
            Ok(quantile)
        } else {
            Err(DeseqError::InvalidOptions {
                reason: "beta prior variance quantile produced non-finite quantile".to_string(),
            })
        }
    }
}

fn weighted_abs_quantile(
    betas: &[f64],
    weights: &[f64],
    probability: f64,
) -> Result<f64, DeseqError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted quantile probability must be finite and between 0 and 1".to_string(),
        });
    }
    let mut pairs = Vec::with_capacity(betas.len());
    for (idx, (beta, weight)) in betas
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
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

    let Some(weight_sum) = checked_sum(pairs.iter().map(|(_, weight)| *weight)) else {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs finite total weight".to_string(),
        });
    };
    if weight_sum <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite total weight"
                .to_string(),
        });
    }

    let Some(norm_scale) = checked_div2(pairs.len() as f64, weight_sum) else {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite normalization scale"
                .to_string(),
        });
    };
    let mut unique = Vec::<(f64, f64)>::with_capacity(pairs.len());
    for (value, weight) in pairs {
        let Some(normalized_weight) = checked_product2(weight, norm_scale) else {
            return Err(DeseqError::InvalidOptions {
                reason:
                    "weighted beta prior variance quantile produced non-finite normalized weight"
                        .to_string(),
            });
        };
        if normalized_weight <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: "weighted beta prior variance quantile needs positive normalized weight"
                    .to_string(),
            });
        }
        if let Some((last_value, last_weight)) = unique.last_mut()
            && *last_value == value
        {
            let Some(next_weight) = checked_sum2(*last_weight, normalized_weight) else {
                return Err(DeseqError::InvalidOptions {
                    reason:
                        "weighted beta prior variance quantile produced non-finite merged weight"
                            .to_string(),
                });
            };
            *last_weight = next_weight;
            continue;
        }
        unique.push((value, normalized_weight));
    }

    if unique.len() == 1 {
        return Ok(unique[0].0);
    }

    let Some(n) = checked_sum(unique.iter().map(|(_, weight)| *weight)) else {
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

    let mut cumulative_weights = Vec::with_capacity(unique.len());
    let mut cumulative = 0.0;
    for (_, weight) in &unique {
        let Some(next_cumulative) = checked_sum2(cumulative, *weight) else {
            return Err(DeseqError::InvalidOptions {
                reason:
                    "weighted beta prior variance quantile produced non-finite cumulative weight"
                        .to_string(),
            });
        };
        cumulative = next_cumulative;
        cumulative_weights.push(cumulative);
    }

    let low_quantile = weighted_order_statistic(&unique, &cumulative_weights, low);
    let high_quantile = weighted_order_statistic(&unique, &cumulative_weights, high);
    let quantile = stable_interpolate(low_quantile, high_quantile, fraction);
    if !quantile.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite quantile"
                .to_string(),
        });
    }
    Ok(quantile)
}

fn stable_interpolate(left: f64, right: f64, fraction: f64) -> f64 {
    if left == right {
        return left;
    }
    let Some(delta) = checked_sum2(right, -left) else {
        return f64::NAN;
    };
    let Some(offset) = checked_product2(fraction, delta) else {
        return f64::NAN;
    };
    checked_sum2(left, offset).unwrap_or(f64::NAN)
}

fn weighted_order_statistic(
    values_and_weights: &[(f64, f64)],
    cumulative_weights: &[f64],
    position: f64,
) -> f64 {
    for (idx, cumulative) in cumulative_weights.iter().copied().enumerate() {
        if position <= cumulative {
            return values_and_weights[idx].0;
        }
    }
    values_and_weights[values_and_weights.len() - 1].0
}
