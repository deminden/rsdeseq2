// Working weights and responses.

fn working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(mu.len());
    working_weights_into(mu, dispersion, weights, &mut output)?;
    Ok(output)
}

fn working_weights_into(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
    output: &mut Vec<f64>,
) -> Result<(), DeseqError> {
    validate_positive_finite(dispersion, "dispersion", 0)?;
    output.clear();
    let values = mu
        .iter()
        .copied()
        .enumerate()
        .map(|(sample, value)| -> Result<f64, DeseqError> {
            validate_positive_finite(value, "mu", sample)?;
            let working_weight = stable_working_weight(value, dispersion).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "IRLS working weight".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                }
            })?;
            Ok(match weights {
                Some(weights) => {
                    let weight = weights[sample];
                    validate_nonnegative_finite(weight, "weight", sample)?;
                    checked_product2(weight, working_weight).ok_or_else(|| {
                        DeseqError::NonFiniteValue {
                            context: "IRLS weighted working weight".to_string(),
                            index: Some(sample),
                            value: f64::NAN,
                        }
                    })?
                }
                None => working_weight,
            })
        });
    for value in values {
        output.push(value?);
    }
    Ok(())
}

fn stable_working_weight(mu: f64, dispersion: f64) -> Option<f64> {
    validate_positive_finite(mu, "mu", 0).ok()?;
    validate_positive_finite(dispersion, "dispersion", 0).ok()?;
    let disp_mu = dispersion * mu;
    if disp_mu.is_infinite() && disp_mu.is_sign_positive() {
        return Some(dispersion.recip());
    }
    let denominator = checked_sum2(1.0, disp_mu)?;
    checked_div2(mu, denominator)
}

#[cfg(test)]
fn working_response(mu: &[f64], nf: &[f64], y: &[f64]) -> Result<Vec<f64>, DeseqError> {
    working_response_with_log_factors(mu, nf, y, None)
}

#[cfg(test)]
fn working_response_with_log_factors(
    mu: &[f64],
    nf: &[f64],
    y: &[f64],
    log_nf: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(mu.len());
    working_response_into(mu, nf, y, log_nf, &mut output)?;
    Ok(output)
}

fn working_response_into(
    mu: &[f64],
    nf: &[f64],
    y: &[f64],
    log_nf: Option<&[f64]>,
    output: &mut Vec<f64>,
) -> Result<(), DeseqError> {
    output.clear();
    let values = mu
        .iter()
        .copied()
        .zip(nf.iter().copied())
        .zip(y.iter().copied())
        .enumerate()
        .map(
            |(sample, ((mu, factor), count))| -> Result<f64, DeseqError> {
                validate_positive_finite(mu, "mu", sample)?;
                validate_positive_finite(factor, "normalization factor", sample)?;
                let normalized_mean = match log_nf {
                    Some(values) => checked_sum2(mu.ln(), -values[sample]),
                    None => stable_log_normalized_mean(mu, factor),
                };
                let log_normalized_mu =
                    normalized_mean.ok_or_else(|| DeseqError::NonFiniteValue {
                        context: "IRLS working response normalized mean".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    })?;
                let residual =
                    checked_sum2(count, -mu).ok_or_else(|| DeseqError::NonFiniteValue {
                        context: "IRLS working response residual".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    })?;
                let scaled_residual =
                    checked_div2(residual, mu).ok_or_else(|| DeseqError::NonFiniteValue {
                        context: "IRLS working response scaled residual".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    })?;
                checked_sum2(log_normalized_mu, scaled_residual).ok_or_else(|| {
                    DeseqError::NonFiniteValue {
                        context: "IRLS working response".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    }
                })
            },
        );
    for value in values {
        output.push(value?);
    }
    Ok(())
}

fn stable_log_normalized_mean(mu: f64, factor: f64) -> Option<f64> {
    validate_positive_finite(mu, "mu", 0).ok()?;
    validate_positive_finite(factor, "normalization factor", 0).ok()?;
    checked_sum2(mu.ln(), -factor.ln())
}
