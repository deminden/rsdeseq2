use super::*;
use crate::numeric::math::ProfileMath;
// Working weights and responses.

pub(super) fn working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(mu.len());
    working_weights_into(mu, dispersion, weights, &mut output)?;
    Ok(output)
}

pub(super) fn working_weights_into(
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
            let working_weight = stable_working_weight(value, dispersion)
                .nonfinite_at("IRLS working weight", Some(sample))?;
            Ok(match weights {
                Some(weights) => {
                    let weight = weights[sample];
                    validate_nonnegative_finite(weight, "weight", sample)?;
                    checked::try_mul(weight, working_weight)
                        .nonfinite_at("IRLS weighted working weight", Some(sample))?
                }
                None => working_weight,
            })
        });
    for value in values {
        output.push(value?);
    }
    Ok(())
}

pub(super) fn stable_working_weight(mu: f64, dispersion: f64) -> checked::Checked {
    for value in [mu, dispersion] {
        if !value.is_finite() || value <= 0.0 {
            return checked::reject(value);
        }
    }
    let disp_mu = dispersion * mu;
    if disp_mu.is_infinite() && disp_mu.is_sign_positive() {
        return Ok(dispersion.recip());
    }
    let denominator = checked::try_add(1.0, disp_mu)?;
    checked::try_div(mu, denominator)
}

#[cfg(test)]
pub(super) fn working_response(mu: &[f64], nf: &[f64], y: &[f64]) -> Result<Vec<f64>, DeseqError> {
    working_response_with_log_factors(mu, nf, y, None)
}

#[cfg(test)]
pub(super) fn working_response_with_log_factors(
    mu: &[f64],
    nf: &[f64],
    y: &[f64],
    log_nf: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(mu.len());
    working_response_into(mu, nf, y, log_nf, &mut output)?;
    Ok(output)
}

pub(super) fn working_response_into(
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
                    Some(values) => checked::try_add(mu.ln_p(), -values[sample]),
                    None => stable_log_normalized_mean(mu, factor),
                };
                let log_normalized_mu = normalized_mean
                    .nonfinite_at("IRLS working response normalized mean", Some(sample))?;
                let residual = checked::try_add(count, -mu)
                    .nonfinite_at("IRLS working response residual", Some(sample))?;
                let scaled_residual = checked::try_div(residual, mu)
                    .nonfinite_at("IRLS working response scaled residual", Some(sample))?;
                checked::try_add(log_normalized_mu, scaled_residual)
                    .nonfinite_at("IRLS working response", Some(sample))
            },
        );
    for value in values {
        output.push(value?);
    }
    Ok(())
}

pub(super) fn stable_log_normalized_mean(mu: f64, factor: f64) -> checked::Checked {
    for value in [mu, factor] {
        if !value.is_finite() || value <= 0.0 {
            return checked::reject(value);
        }
    }
    checked::try_add(mu.ln_p(), -factor.ln_p())
}
