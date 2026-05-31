/// DESeq2's log-dispersion prior kernel, omitting additive constants.
pub fn dispersion_prior_log_density(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    let residual = checked_sub(
        log_alpha,
        prior.log_mean,
        0,
        "dispersion prior log residual",
    )?;
    let residual_square = checked_mul(
        residual,
        residual,
        0,
        "dispersion prior log residual square",
    )?;
    Ok(-0.5 * residual_square / prior.variance)
}

/// Derivative of the log-dispersion prior kernel with respect to log alpha.
pub fn dispersion_prior_derivative(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    Ok(-(log_alpha - prior.log_mean) / prior.variance)
}

/// Second derivative of the log-dispersion prior kernel with respect to log alpha.
pub fn dispersion_prior_second_derivative(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    Ok(-prior.variance.recip())
}

/// Dispersion log posterior without prior and with optional Cox-Reid correction.
pub fn dispersion_log_posterior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_with_prior(counts, mu, design, log_alpha, use_cox_reid, None)
}

/// Dispersion log posterior with optional Cox-Reid correction and log-alpha prior.
pub fn dispersion_log_posterior_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Dispersion log posterior with optional Cox-Reid correction, log-alpha prior, and weights.
pub fn dispersion_log_posterior_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let posterior = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion objective requires a design matrix".to_string(),
            });
        };
        checked_scaled_sum(
            &[
                likelihood,
                cox_reid_adjustment_weighted_with_threshold(
                    design,
                    input.mu,
                    log_alpha,
                    input.weights,
                    input.weight_threshold,
                )?,
            ],
            "dispersion log posterior Cox-Reid sum",
        )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        checked_scaled_sum(
            &[posterior, dispersion_prior_log_density(log_alpha, prior)?],
            "dispersion log posterior prior sum",
        )
    } else {
        Ok(posterior)
    }
}

/// Derivative of the dispersion log posterior with respect to log alpha.
pub fn dispersion_log_posterior_derivative(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_with_prior(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        None,
    )
}

/// Derivative of the dispersion log posterior with an optional log-alpha prior.
pub fn dispersion_log_posterior_derivative_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Derivative of the dispersion log posterior with optional prior and weights.
pub fn dispersion_log_posterior_derivative_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_derivative_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_derivative_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let derivative = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion derivative requires a design matrix".to_string(),
            });
        };
        checked_scaled_sum(
            &[
                likelihood,
                cox_reid_adjustment_derivative_weighted_with_threshold(
                    design,
                    input.mu,
                    log_alpha,
                    input.weights,
                    input.weight_threshold,
                )?,
            ],
            "dispersion log posterior derivative Cox-Reid sum",
        )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        checked_scaled_sum(
            &[derivative, dispersion_prior_derivative(log_alpha, prior)?],
            "dispersion log posterior derivative prior sum",
        )
    } else {
        Ok(derivative)
    }
}

/// Second derivative of the dispersion log posterior with respect to log alpha.
pub fn dispersion_log_posterior_second_derivative(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_with_prior(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        None,
    )
}

/// Second derivative of the dispersion log posterior with an optional log-alpha prior.
pub fn dispersion_log_posterior_second_derivative_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Second derivative of the dispersion log posterior with optional prior and weights.
pub fn dispersion_log_posterior_second_derivative_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_second_derivative_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let second_derivative = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion second derivative requires a design matrix"
                    .to_string(),
            });
        };
        checked_scaled_sum(
            &[
                likelihood,
                cox_reid_adjustment_second_derivative_weighted_with_threshold(
                    design,
                    input.mu,
                    log_alpha,
                    input.weights,
                    input.weight_threshold,
                )?,
            ],
            "dispersion log posterior second derivative Cox-Reid sum",
        )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        checked_scaled_sum(
            &[
                second_derivative,
                dispersion_prior_second_derivative(log_alpha, prior)?,
            ],
            "dispersion log posterior second derivative prior sum",
        )
    } else {
        Ok(second_derivative)
    }
}
