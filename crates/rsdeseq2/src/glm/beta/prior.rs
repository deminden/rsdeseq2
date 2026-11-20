use super::*;
use crate::errors::{CalculationLocation, ErrorKind, invalid_calculation};
// Prior-variance estimation and upper-quantile matching.

/// Estimate DESeq2-style beta prior variances from unshrunken MLE betas.
///
/// This mirrors the computational core of DESeq2 `estimateBetaPriorVar` for
/// already-built primitive matrices: each coefficient uses MLE betas whose
/// absolute value is below the finite-beta cutoff, then matches an upper
/// absolute-beta quantile to a zero-centered Normal. The weighted method uses
/// DESeq2's `1 / (1 / baseMean + dispFit)` row weights. Intercept columns are
/// assigned the configured wide prior variance.
pub fn estimate_beta_prior_variance(
    beta_matrix: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    coefficient_names: Option<&[String]>,
    options: BetaPriorVarianceOptions,
) -> Result<Vec<f64>, DeseqError> {
    validate_beta_prior_inputs(beta_matrix, base_mean, disp_fit, coefficient_names, options)?;
    let weights = match options.method {
        BetaPriorVarianceMethod::Weighted => Some(beta_prior_weights(
            base_mean,
            disp_fit,
            beta_matrix.n_rows(),
        )?),
        BetaPriorVarianceMethod::Quantile => None,
    };

    let mut prior_variance = Vec::with_capacity(beta_matrix.n_cols());
    for coefficient in 0..beta_matrix.n_cols() {
        let value = if beta_matrix.n_rows() == 1 {
            let beta = beta_matrix.row(0)?[coefficient];
            let is_intercept = coefficient_names.is_some_and(|names| {
                matches!(names[coefficient].as_str(), "Intercept" | "(Intercept)")
            });
            if beta.is_finite() && !is_intercept {
                checked::try_mul(beta, beta).map_err(|failure| {
                    invalid_calculation(
                        ErrorKind::Options,
                        "single-gene beta prior variance",
                        Some(CalculationLocation::Coefficient(coefficient)),
                        failure.value,
                        "single-gene beta prior produced non-finite variance",
                    )
                })?
            } else {
                options.wide_prior_variance
            }
        } else {
            let mut betas = Vec::new();
            let mut selected_weights = Vec::new();
            for row in 0..beta_matrix.n_rows() {
                let beta = beta_matrix.row(row)?[coefficient];
                if beta.is_finite() && beta.abs() < options.finite_beta_cutoff {
                    betas.push(beta);
                    if let Some(weights) = &weights {
                        selected_weights.push(weights[row]);
                    }
                }
            }
            if betas.is_empty() {
                options.wide_prior_variance
            } else {
                let variance = match options.method {
                    BetaPriorVarianceMethod::Quantile => {
                        match_upper_quantile_for_variance(&betas, options.upper_quantile)
                    }
                    BetaPriorVarianceMethod::Weighted => {
                        match_weighted_upper_quantile_for_variance(
                            &betas,
                            &selected_weights,
                            options.upper_quantile,
                        )
                    }
                };
                variance.map_err(|mut error| {
                    if let DeseqError::InvalidCalculation { location, .. } = &mut error {
                        *location = Some(CalculationLocation::Coefficient(coefficient));
                    }
                    error
                })?
            }
        };
        prior_variance.push(value);
    }

    if let Some(names) = coefficient_names {
        for (idx, name) in names.iter().enumerate() {
            if name == "Intercept" || name == "(Intercept)" {
                prior_variance[idx] = options.wide_prior_variance;
            }
        }
    }
    Ok(prior_variance)
}

/// Match an upper absolute-beta quantile to a zero-centered Normal variance.
pub fn match_upper_quantile_for_variance(
    betas: &[f64],
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    validate_upper_quantile(upper_quantile)?;
    let abs_betas = finite_abs_values(betas)?;
    let quantile = quantile_type7_stable(abs_betas, 1.0 - upper_quantile)?;
    variance_from_quantile(quantile, upper_quantile, "beta prior variance")
}

/// Weighted version of [`match_upper_quantile_for_variance`].
pub fn match_weighted_upper_quantile_for_variance(
    betas: &[f64],
    weights: &[f64],
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    validate_upper_quantile(upper_quantile)?;
    if betas.len() != weights.len() {
        return Err(invalid_dimensions(
            "beta prior variance weights",
            betas.len(),
            weights.len(),
        ));
    }
    let weighted_quantile = weighted_abs_quantile(betas, weights, 1.0 - upper_quantile)?;
    variance_from_quantile(
        weighted_quantile,
        upper_quantile,
        "weighted beta prior variance",
    )
}

fn variance_from_quantile(
    quantile: f64,
    upper_quantile: f64,
    operation: &'static str,
) -> Result<f64, DeseqError> {
    let normal_quantile = standard_normal_inverse_cdf(1.0 - upper_quantile / 2.0);
    let scale = quantile / normal_quantile;
    checked::try_mul(scale, scale).map_err(|failure| {
        invalid_calculation(
            ErrorKind::Options,
            operation,
            None,
            failure.value,
            format_args!("{operation} quantile produced non-finite variance"),
        )
    })
}
