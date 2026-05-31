/// DESeq2's alpha-dependent NB log-likelihood kernel.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer.
pub fn dispersion_nb_log_likelihood_kernel(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_weighted(counts, mu, log_alpha, None)
}

/// DESeq2's alpha-dependent NB log-likelihood kernel with optional observation weights.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer. When supplied, observation weights multiply
/// the per-sample terms.
pub fn dispersion_nb_log_likelihood_kernel_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(weights, counts.len(), "dispersion objective weights")?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let mut total = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu * alpha;
        let mu_plus_inv_alpha = checked_add(
            mu,
            inv_alpha,
            sample,
            "dispersion objective mean plus inverse alpha",
        )?;
        let term = ln_gamma(y + inv_alpha)
            - ln_gamma(inv_alpha)
            - y * mu_plus_inv_alpha.ln()
            - inv_alpha * mu_alpha.ln_1p();
        if !term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective likelihood term".to_string(),
                index: Some(sample),
                value: term,
            });
        }
        checked_matrix_add_assign(
            &mut total,
            checked_mul(
                observation_weight,
                term,
                sample,
                "dispersion objective weighted likelihood term",
            )?,
            sample,
            "dispersion objective likelihood sum",
        )?;
    }
    Ok(total)
}

/// Derivative of DESeq2's alpha-dependent NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, None)
}

/// Derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective derivative weights",
    )?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let mut derivative_alpha = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu_alpha_terms(mu, alpha, sample, "dispersion objective derivative")?;
        let term = digamma(inv_alpha) + mu_alpha.log1p - mu_alpha.ratio - digamma(y + inv_alpha)
            + y * mu_alpha.alpha_over_one_plus;
        if !term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective derivative term".to_string(),
                index: Some(sample),
                value: term,
            });
        }
        checked_matrix_add_assign(
            &mut derivative_alpha,
            checked_mul(
                observation_weight,
                term,
                sample,
                "dispersion objective weighted derivative term",
            )?,
            sample,
            "dispersion objective derivative sum",
        )?;
    }
    checked_log_alpha_first_derivative(
        inv_alpha,
        derivative_alpha,
        "dispersion objective log-alpha derivative",
    )
}

/// Second derivative of DESeq2's NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_second_derivative_weighted(counts, mu, log_alpha, None)
}

/// Second derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective second derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective second derivative weights",
    )?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let inv_alpha_squared = inv_alpha * inv_alpha;
    let mut first_alpha_sum = 0.0;
    let mut second_alpha_sum = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu_alpha_terms(mu, alpha, sample, "dispersion objective second derivative")?;
        let first_term =
            digamma(inv_alpha) + mu_alpha.log1p - mu_alpha.ratio - digamma(y + inv_alpha)
                + y * mu_alpha.alpha_over_one_plus;
        let second_term = -inv_alpha_squared * trigamma(inv_alpha)?
            + mu_alpha.mu_squared_alpha_over_one_plus_squared
            + inv_alpha_squared * trigamma(y + inv_alpha)?
            + y * mu_alpha.inv_one_plus_squared;
        if !first_term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective second derivative first term".to_string(),
                index: Some(sample),
                value: first_term,
            });
        }
        if !second_term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective second derivative term".to_string(),
                index: Some(sample),
                value: second_term,
            });
        }
        checked_matrix_add_assign(
            &mut first_alpha_sum,
            checked_mul(
                observation_weight,
                first_term,
                sample,
                "dispersion objective weighted first derivative term",
            )?,
            sample,
            "dispersion objective first derivative sum",
        )?;
        checked_matrix_add_assign(
            &mut second_alpha_sum,
            checked_mul(
                observation_weight,
                second_term,
                sample,
                "dispersion objective weighted second derivative term",
            )?,
            sample,
            "dispersion objective second derivative sum",
        )?;
    }
    let first_log_alpha =
        dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, weights)?;
    checked_log_alpha_second_derivative(
        second_alpha_sum,
        inv_alpha,
        first_alpha_sum,
        first_log_alpha,
        "dispersion objective log-alpha second derivative",
    )
}
