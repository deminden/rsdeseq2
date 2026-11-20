use super::*;
use crate::errors::{ErrorKind, Requirement, invalid_parameter, validate_matrix_shape};
use crate::numeric::math::ProfileMath;
pub(super) fn best_log_alpha(
    objective: DispersionObjectiveInput<'_>,
    grid: &[f64],
) -> Result<(f64, f64), DeseqError> {
    let mut best_log = grid[0];
    let mut best_score = dispersion_log_posterior_objective(objective, best_log)?;
    for log_alpha in grid.iter().copied().skip(1) {
        let score = dispersion_log_posterior_objective(objective, log_alpha)?;
        if score > best_score {
            best_log = log_alpha;
            best_score = score;
        }
    }
    Ok((best_log, best_score))
}

pub(super) fn linspace(start: f64, end: f64, len: usize) -> Result<Vec<f64>, DeseqError> {
    if len == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion grid points".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !start.is_finite() || !end.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion grid endpoint".to_string(),
            index: None,
            value: if start.is_finite() { end } else { start },
        });
    }
    if len == 1 {
        return Ok(vec![start]);
    }
    let span = checked::try_sub(end, start).nonfinite_at("dispersion grid span", Some(0))?;
    let step = span / (len as f64 - 1.0);
    if !step.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion grid step".to_string(),
            index: None,
            value: step,
        });
    }
    (0..len)
        .map(|idx| {
            let offset = checked::try_mul(step, idx as f64)
                .nonfinite_at("dispersion grid offset", Some(idx))?;
            checked::try_add(start, offset).nonfinite_at("dispersion grid value", Some(idx))
        })
        .collect()
}

pub(super) fn validate_gene_est_inputs(
    input: &GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<(), DeseqError> {
    if input.design.n_samples() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion design samples",
            input.counts.n_samples(),
            input.design.n_samples(),
        ));
    }
    validate_matrix_shape(
        "dispersion normalized counts",
        (input.counts.n_genes(), input.counts.n_samples()),
        (
            input.normalized_counts.n_rows(),
            input.normalized_counts.n_cols(),
        ),
    )?;
    if input.base_mean.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseMean",
            input.counts.n_genes(),
            input.base_mean.len(),
        ));
    }
    if input.base_var.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseVar",
            input.counts.n_genes(),
            input.base_var.len(),
        ));
    }
    if input.all_zero.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion allZero",
            input.counts.n_genes(),
            input.all_zero.len(),
        ));
    }
    if input.size_factors.len() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion size factors",
            input.counts.n_samples(),
            input.size_factors.len(),
        ));
    }
    validate_size_factor_values(input.size_factors)?;
    if let Some(normalization_factors) = input.normalization_factors {
        validate_normalization_factor_matrix(
            normalization_factors,
            input.counts.n_genes(),
            input.counts.n_samples(),
            "dispersion normalization factors",
        )?;
    }
    if let Some(observation_weights) = input.observation_weights {
        validate_matrix_shape(
            "dispersion observation weights",
            (input.counts.n_genes(), input.counts.n_samples()),
            (observation_weights.n_rows(), observation_weights.n_cols()),
        )?;
        for (idx, weight) in observation_weights.as_slice().iter().copied().enumerate() {
            if !weight.is_finite() || weight < 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "dispersion observation weight".to_string(),
                    index: Some(idx),
                    value: weight,
                });
            }
        }
    }
    validate_gene_est_options(options)?;
    validate_dispersion_bounds(
        options.min_disp,
        max_dispersion(options, input.counts.n_samples()),
    )?;
    Ok(())
}

pub(super) fn validate_gene_est_options(
    options: GeneWiseDispersionOptions,
) -> Result<(), DeseqError> {
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "min_disp",
            Requirement::FiniteAndPositive,
        ));
    }
    if (options.min_disp / 10.0).ln_p() <= -30.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "log(min_disp / 10) must be above -30 for numerical stability".to_string(),
        });
    }
    if !options.min_mu.is_finite() || options.min_mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion min_mu".to_string(),
            index: None,
            value: options.min_mu,
        });
    }
    if options.grid_points < 3 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion grid points".to_string(),
            expected: 3,
            actual: options.grid_points,
        });
    }
    if !options.kappa_0.is_finite() || options.kappa_0 <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "kappa_0",
            Requirement::FiniteAndPositive,
        ));
    }
    if !options.disp_tol.is_finite() || options.disp_tol <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "disp_tol",
            Requirement::FiniteAndPositive,
        ));
    }
    if options.maxit == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion maxit".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if options.niter == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion niter".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    validate_weight_threshold(options.weight_threshold, "dispersion weight_threshold")?;
    Ok(())
}

pub(super) fn validate_weight_threshold(value: f64, context: &str) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            context,
            Requirement::FiniteAndNonNegative,
        ));
    }
    Ok(())
}

pub(super) fn validate_dispersion_prior(prior: Option<DispersionPrior>) -> Result<(), DeseqError> {
    if let Some(prior) = prior {
        if !prior.log_mean.is_finite() {
            return Err(invalid_parameter(
                ErrorKind::Dispersion,
                "dispersion prior log_mean",
                Requirement::Finite,
            ));
        }
        if !prior.variance.is_finite() || prior.variance <= 0.0 {
            return Err(invalid_parameter(
                ErrorKind::Dispersion,
                "dispersion prior variance",
                Requirement::FiniteAndPositive,
            ));
        }
    }
    Ok(())
}

#[inline(always)]
pub(super) fn validate_positive_mu(mu: f64, sample: usize) -> Result<(), DeseqError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion mean".to_string(),
            index: Some(sample),
            value: mu,
        });
    }
    Ok(())
}

#[inline(always)]
pub(super) fn checked_matrix_add_assign(
    sum: &mut f64,
    term: f64,
    index: usize,
    context: &str,
) -> Result<(), DeseqError> {
    let value = *sum + term;
    if value.is_finite() {
        *sum = value;
        Ok(())
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

pub(super) fn checked_sum_indexed(
    values: impl IntoIterator<Item = f64>,
    context: &str,
) -> Result<f64, DeseqError> {
    let mut sum = 0.0;
    for (idx, value) in values.into_iter().enumerate() {
        checked_matrix_add_assign(&mut sum, value, idx, context)?;
    }
    Ok(sum)
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MuAlphaTerms {
    pub(super) log1p: f64,
    pub(super) ratio: f64,
    pub(super) alpha_over_one_plus: f64,
    pub(super) mu_squared_alpha_over_one_plus_squared: f64,
    pub(super) inv_one_plus_squared: f64,
}

pub(super) fn mu_alpha_terms(
    mu: f64,
    alpha: f64,
    index: usize,
    context: &str,
) -> Result<MuAlphaTerms, DeseqError> {
    mu_alpha_terms_impl::<true>(mu, alpha, index, context)
}

#[inline(always)]
pub(super) fn mu_alpha_terms_impl<const NEED_LOG: bool>(
    mu: f64,
    alpha: f64,
    index: usize,
    context: &str,
) -> Result<MuAlphaTerms, DeseqError> {
    let mu_alpha = mu * alpha;
    let terms = if mu_alpha.is_finite() {
        let one_plus = 1.0 + mu_alpha;
        let inv_one_plus = one_plus.recip();
        let ratio = mu_alpha * inv_one_plus;
        let alpha_over_one_plus = alpha * inv_one_plus;
        let inv_one_plus_squared = checked_mu_alpha_product(
            inv_one_plus,
            inv_one_plus,
            index,
            context,
            "inverse denominator square",
        )?;
        let ratio_squared =
            checked_mu_alpha_product(ratio, ratio, index, context, "mean-dispersion ratio square")?;
        let mu_squared_alpha_over_one_plus_squared = checked_mu_alpha_product(
            ratio_squared,
            alpha.recip(),
            index,
            context,
            "mean curvature term",
        )?;
        MuAlphaTerms {
            log1p: if NEED_LOG {
                if crate::numeric_profile() == crate::NumericProfile::Portable {
                    mu_alpha.ln_1p_p()
                } else {
                    (1.0 + mu_alpha).ln_p()
                }
            } else {
                0.0
            },
            ratio,
            alpha_over_one_plus,
            mu_squared_alpha_over_one_plus_squared,
            inv_one_plus_squared,
        }
    } else {
        let log1p = if NEED_LOG {
            mu.ln_p() + alpha.ln_p()
        } else {
            0.0
        };
        let alpha_over_one_plus = mu.recip();
        let mu_squared_alpha_over_one_plus_squared = alpha.recip();
        MuAlphaTerms {
            log1p,
            ratio: 1.0,
            alpha_over_one_plus,
            mu_squared_alpha_over_one_plus_squared,
            inv_one_plus_squared: 0.0,
        }
    };
    for value in [
        terms.log1p,
        terms.ratio,
        terms.alpha_over_one_plus,
        terms.mu_squared_alpha_over_one_plus_squared,
        terms.inv_one_plus_squared,
    ] {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(index),
                value,
            });
        }
    }
    Ok(terms)
}

#[inline(always)]
pub(super) fn checked_mu_alpha_product(
    left: f64,
    right: f64,
    index: usize,
    context: &str,
    suffix: &str,
) -> Result<f64, DeseqError> {
    let value = left * right;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: format!("{context} {suffix}"),
            index: Some(index),
            value,
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct CoxReidWeightTerms {
    pub(super) weight: f64,
    pub(super) d_weight: f64,
    pub(super) d2_weight: f64,
}

#[inline(always)]
pub(super) fn cox_reid_weight_terms(
    mu: f64,
    alpha: f64,
    index: usize,
) -> Result<CoxReidWeightTerms, DeseqError> {
    // Positive finite inputs make the unused logarithm finite. Keep all
    // remaining curvature/weight checks and the original invalid-input path.
    let mu_alpha = if mu.is_finite() && mu > 0.0 && alpha.is_finite() && alpha > 0.0 {
        mu_alpha_terms_impl::<false>(mu, alpha, index, "Cox-Reid working weight")?
    } else {
        mu_alpha_terms(mu, alpha, index, "Cox-Reid working weight")?
    };
    let weight = checked::try_mul(mu_alpha.ratio, alpha.recip())
        .nonfinite_at("Cox-Reid working weight", Some(index))?;
    let weight_square = checked::try_mul(weight, weight)
        .nonfinite_at("Cox-Reid working weight square", Some(index))?;
    let d_weight = -weight_square;
    let d2_weight = checked::try_mul(
        2.0,
        checked::try_mul(weight_square, weight)
            .nonfinite_at("Cox-Reid working weight cube", Some(index))?,
    )
    .nonfinite_at("Cox-Reid working second derivative weight", Some(index))?;
    Ok(CoxReidWeightTerms {
        weight,
        d_weight,
        d2_weight,
    })
}

pub(super) fn checked_log_alpha_first_derivative(
    inv_alpha: f64,
    derivative_alpha: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    checked::try_mul(inv_alpha, derivative_alpha).nonfinite_at(context, Some(0))
}

pub(super) fn checked_log_alpha_second_derivative(
    second_alpha_sum: f64,
    inv_alpha: f64,
    first_alpha_sum: f64,
    first_log_alpha: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    let alpha_first_scale = checked::try_mul(2.0, inv_alpha).nonfinite_at(context, Some(0))?;
    let alpha_first_term =
        -checked::try_mul(alpha_first_scale, first_alpha_sum).nonfinite_at(context, Some(0))?;
    checked::try_scaled_sum(&[second_alpha_sum, alpha_first_term, first_log_alpha])
        .nonfinite_at(context, None)
}

pub(super) fn checked_cox_reid_log_alpha_second_derivative(
    second_alpha: f64,
    alpha: f64,
    first_log_alpha: f64,
) -> Result<f64, DeseqError> {
    let alpha_squared = checked::try_mul(alpha, alpha)
        .nonfinite_at("Cox-Reid log-alpha second derivative alpha square", Some(0))?;
    let alpha_term = checked::try_mul(second_alpha, alpha_squared)
        .nonfinite_at("Cox-Reid log-alpha second derivative alpha term", Some(0))?;
    checked::try_scaled_sum(&[alpha_term, first_log_alpha])
        .nonfinite_at("Cox-Reid log-alpha second derivative", None)
}

pub(super) fn validate_observation_weight_slice(
    weights: Option<&[f64]>,
    expected_len: usize,
    context: &str,
) -> Result<(), DeseqError> {
    let Some(weights) = weights else {
        return Ok(());
    };
    if weights.len() != expected_len {
        return Err(invalid_dimensions(context, expected_len, weights.len()));
    }
    for (idx, weight) in weights.iter().copied().enumerate() {
        if !weight.is_finite() || weight < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value: weight,
            });
        }
    }
    Ok(())
}

pub(super) fn validate_dispersion_bounds(min_disp: f64, max_disp: f64) -> Result<(), DeseqError> {
    if !max_disp.is_finite() || max_disp <= min_disp {
        return Err(DeseqError::InvalidDispersion {
            reason: "max dispersion must be finite and greater than min dispersion".to_string(),
        });
    }
    Ok(())
}

pub(super) fn max_dispersion(options: GeneWiseDispersionOptions, n_samples: usize) -> f64 {
    options
        .max_disp
        .unwrap_or_else(|| 10.0_f64.max(n_samples as f64))
}
