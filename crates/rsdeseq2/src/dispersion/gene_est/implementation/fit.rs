/// Combine rough and moments estimates using DESeq2's bounded start shape.
pub fn initial_dispersion_estimates(
    rough_disp: &[f64],
    moments_disp: &[f64],
    min_disp: f64,
    max_disp: f64,
) -> Result<Vec<f64>, DeseqError> {
    if rough_disp.len() != moments_disp.len() {
        return Err(invalid_dimensions(
            "initial dispersion starts",
            rough_disp.len(),
            moments_disp.len(),
        ));
    }
    validate_dispersion_bounds(min_disp, max_disp)?;
    Ok(rough_disp
        .iter()
        .copied()
        .zip(moments_disp.iter().copied())
        .map(|(rough, moments)| {
            if !rough.is_finite() || !moments.is_finite() {
                f64::NAN
            } else {
                rough.min(moments).clamp(min_disp, max_disp)
            }
        })
        .collect())
}

fn fit_dispersion_for_gene(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize, bool), DeseqError> {
    let diagnostics = fit_dispersion_for_gene_detailed(
        counts,
        mu,
        design,
        initial_dispersion,
        options,
        n_samples,
    )?;
    Ok((
        diagnostics.estimate,
        diagnostics.iterations,
        diagnostics.converged,
    ))
}

fn fit_dispersion_for_gene_detailed(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<GeneDispersionFitDiagnostics, DeseqError> {
    fit_dispersion_for_gene_detailed_with_weights(
        counts,
        mu,
        design,
        initial_dispersion,
        options,
        n_samples,
        None,
    )
}

fn fit_dispersion_for_gene_detailed_with_weights(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    weights: Option<&[f64]>,
) -> Result<GeneDispersionFitDiagnostics, DeseqError> {
    match options.fit_method {
        GeneWiseDispersionFitMethod::Grid => {
            let (dispersion, evaluations) = fit_dispersion_grid_inner(DispersionOptimizerInput {
                counts,
                mu,
                design: Some(design),
                initial_dispersion,
                options,
                n_samples,
                prior: None,
                weights,
            })?;
            let last_lp = dispersion_log_posterior_objective(
                DispersionObjectiveInput {
                    counts,
                    mu,
                    design: Some(design),
                    use_cox_reid: options.use_cox_reid,
                    prior: None,
                    weights,
                    weight_threshold: options.weight_threshold,
                },
                dispersion.ln(),
            )?;
            Ok(GeneDispersionFitDiagnostics {
                estimate: dispersion,
                iterations: evaluations,
                converged: true,
                initial_lp: last_lp,
                last_lp,
            })
        }
        GeneWiseDispersionFitMethod::LineSearch => {
            let line_search = fit_dispersion_line_search_inner(DispersionOptimizerInput {
                counts,
                mu,
                design: Some(design),
                initial_dispersion,
                options,
                n_samples,
                prior: None,
                weights,
            })?;
            let mut dispersion = line_search.dispersion;
            if !line_search.converged && dispersion > options.min_disp * 10.0 {
                dispersion = fit_dispersion_grid_inner(DispersionOptimizerInput {
                    counts,
                    mu,
                    design: Some(design),
                    initial_dispersion: dispersion,
                    options,
                    n_samples,
                    prior: None,
                    weights,
                })?
                .0;
            }
            Ok(GeneDispersionFitDiagnostics {
                estimate: dispersion,
                iterations: line_search.iter,
                converged: line_search.converged,
                initial_lp: line_search.initial_lp,
                last_lp: line_search.last_lp,
            })
        }
    }
}

/// Fit one dispersion by DESeq2's Armijo line-search shape.
pub fn fit_dispersion_line_search(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit one dispersion by DESeq2's Armijo line-search shape with a log-alpha prior.
pub fn fit_dispersion_line_search_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit one dispersion by DESeq2's weighted Armijo line-search shape with a log-alpha prior.
pub fn fit_dispersion_line_search_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts: input.counts,
        mu: input.mu,
        design: Some(input.design),
        initial_dispersion: input.initial_dispersion,
        options: input.options,
        n_samples: input.n_samples,
        prior: Some(input.prior),
        weights: Some(input.weights),
    })
}

/// Fit one dispersion by line search without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit one dispersion by line search with a log-alpha prior and without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit one dispersion by weighted line search with a log-alpha prior and without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: Some(weights),
    })
}

#[derive(Clone, Copy)]
struct DispersionOptimizerInput<'a> {
    counts: &'a [u32],
    mu: &'a [f64],
    design: Option<&'a DesignMatrix>,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: Option<DispersionPrior>,
    weights: Option<&'a [f64]>,
}

#[derive(Clone, Copy)]
struct DispersionObjectiveInput<'a> {
    counts: &'a [u32],
    mu: &'a [f64],
    design: Option<&'a DesignMatrix>,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&'a [f64]>,
    weight_threshold: f64,
}

fn fit_dispersion_line_search_inner(
    input: DispersionOptimizerInput<'_>,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    let counts = input.counts;
    let mu = input.mu;
    let design = input.design;
    let initial_dispersion = input.initial_dispersion;
    let mut options = input.options;
    let n_samples = input.n_samples;
    let prior = input.prior;
    let weights = input.weights;
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion line-search mu",
            counts.len(),
            mu.len(),
        ));
    }
    if design.is_none() {
        options.use_cox_reid = false;
    }
    validate_gene_est_options(options)?;
    validate_dispersion_prior(prior)?;
    let max_disp = max_dispersion(options, n_samples);
    validate_dispersion_bounds(options.min_disp, max_disp)?;
    if !initial_dispersion.is_finite() || initial_dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "initial dispersion must be finite and positive".to_string(),
        });
    }

    let min_log_alpha = (options.min_disp / 10.0).ln();
    let mut log_alpha = initial_dispersion.max(options.min_disp).ln();
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
    };
    let mut lp = dispersion_log_posterior_objective(objective, log_alpha)?;
    let mut dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
    let initial_lp = lp;
    let initial_dlp = dlp;
    let mut kappa = options.kappa_0;
    let mut iter = 0_usize;
    let mut iter_accept = 0_usize;
    let mut last_change = -1.0;
    let epsilon = 1.0e-4;

    for _ in 0..options.maxit {
        iter += 1;
        if !dlp.is_finite() || !kappa.is_finite() || kappa <= 0.0 {
            break;
        }

        let Some((proposed_log_alpha, effective_kappa)) =
            bounded_log_alpha_proposal(log_alpha, dlp, kappa, -30.0, 10.0)
        else {
            break;
        };

        let theta_kappa = -dispersion_log_posterior_objective(objective, proposed_log_alpha)?;
        let theta_hat_kappa = checked_line_search_armijo_bound(lp, effective_kappa, epsilon, dlp)?;
        if theta_kappa <= theta_hat_kappa {
            iter_accept += 1;
            log_alpha = proposed_log_alpha;
            let lp_new = dispersion_log_posterior_objective(objective, log_alpha)?;
            last_change = checked_sub(
                lp_new,
                lp,
                0,
                "dispersion line-search accepted objective change",
            )?;
            lp = lp_new;
            if last_change < options.disp_tol {
                break;
            }
            if log_alpha < min_log_alpha {
                break;
            }
            dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
            kappa = (effective_kappa * 1.1).min(options.kappa_0);
            if iter_accept.is_multiple_of(5) {
                kappa /= 2.0;
            }
        } else {
            kappa = effective_kappa / 2.0;
        }
    }

    let dispersion = log_alpha.exp().clamp(options.min_disp, max_disp);
    let last_dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
    let last_d2lp = dispersion_log_posterior_second_derivative_objective(objective, log_alpha)?;
    Ok(DispersionLineSearchOutput {
        dispersion,
        log_alpha,
        iter,
        iter_accept,
        initial_lp,
        initial_dlp,
        last_lp: lp,
        last_dlp,
        last_d2lp,
        last_change,
        converged: iter < options.maxit && iter != 1,
    })
}

fn bounded_log_alpha_proposal(
    log_alpha: f64,
    direction: f64,
    step: f64,
    lower: f64,
    upper: f64,
) -> Option<(f64, f64)> {
    if !log_alpha.is_finite()
        || !direction.is_finite()
        || !step.is_finite()
        || !lower.is_finite()
        || !upper.is_finite()
        || direction == 0.0
        || step <= 0.0
        || lower >= upper
    {
        return None;
    }
    let unclamped = checked_mul(step, direction, 0, "dispersion line-search proposal step")
        .ok()
        .and_then(|movement| {
            checked_add(log_alpha, movement, 0, "dispersion line-search proposal").ok()
        })?;
    let clamped = unclamped.clamp(lower, upper);
    let effective_step = checked_sub(
        clamped,
        log_alpha,
        0,
        "dispersion line-search effective proposal movement",
    )
    .ok()
    .and_then(|movement| {
        checked_div(
            movement,
            direction,
            0,
            "dispersion line-search effective proposal step",
        )
        .ok()
    })?;
    if effective_step > 0.0 {
        Some((clamped, effective_step))
    } else {
        None
    }
}

fn checked_line_search_armijo_bound(
    lp: f64,
    effective_kappa: f64,
    epsilon: f64,
    dlp: f64,
) -> Result<f64, DeseqError> {
    let dlp_square = checked_mul(dlp, dlp, 0, "dispersion line-search Armijo slope square")?;
    let scaled_slope = checked_mul(
        effective_kappa,
        epsilon,
        0,
        "dispersion line-search Armijo scale",
    )
    .and_then(|scale| checked_mul(scale, dlp_square, 0, "dispersion line-search Armijo scale"))?;
    checked_sub(-lp, scaled_slope, 0, "dispersion line-search Armijo bound")
}

/// Fit a dispersion for one gene by DESeq2-style two-pass log-alpha grid search.
pub fn fit_dispersion_grid(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit a dispersion by two-pass log-alpha grid search with a log-alpha prior.
pub fn fit_dispersion_grid_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit a dispersion by weighted two-pass log-alpha grid search with a log-alpha prior.
pub fn fit_dispersion_grid_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts: input.counts,
        mu: input.mu,
        design: Some(input.design),
        initial_dispersion: input.initial_dispersion,
        options: input.options,
        n_samples: input.n_samples,
        prior: Some(input.prior),
        weights: Some(input.weights),
    })
}

/// Fit a dispersion for one gene without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit a dispersion by two-pass log-alpha grid search with a prior and without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit a dispersion by weighted two-pass log-alpha grid search with a prior and without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: Some(weights),
    })
}

fn fit_dispersion_grid_inner(
    input: DispersionOptimizerInput<'_>,
) -> Result<(f64, usize), DeseqError> {
    let counts = input.counts;
    let mu = input.mu;
    let design = input.design;
    let mut options = input.options;
    let n_samples = input.n_samples;
    let prior = input.prior;
    let weights = input.weights;
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion grid mu",
            counts.len(),
            mu.len(),
        ));
    }
    if design.is_none() {
        options.use_cox_reid = false;
    }
    validate_gene_est_options(options)?;
    validate_dispersion_prior(prior)?;
    if options.use_cox_reid && design.is_none() {
        return Err(DeseqError::UnsupportedFeature {
            feature: "Cox-Reid dispersion fitting requires a design matrix".to_string(),
        });
    }
    let max_disp = max_dispersion(options, n_samples);
    validate_dispersion_bounds(options.min_disp, max_disp)?;
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
    };
    let min_log = options.min_disp.ln();
    let max_log = max_disp.ln();
    let coarse = linspace(min_log, max_log, options.grid_points)?;
    let (best_log, _) = best_log_alpha(objective, &coarse)?;
    let delta = checked_sub(coarse[1], coarse[0], 1, "dispersion grid step")?;
    let fine_lower = checked_sub(best_log, delta, 0, "dispersion fine grid lower bound")?;
    let fine_upper = checked_add(best_log, delta, 0, "dispersion fine grid upper bound")?;
    let fine = linspace(fine_lower, fine_upper, options.grid_points)?;
    let (best_fine_log, _) = best_log_alpha(objective, &fine)?;
    Ok((
        best_fine_log.exp().clamp(options.min_disp, max_disp),
        options.grid_points * 2,
    ))
}
