use super::*;
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

pub(super) fn fit_dispersion_for_gene(
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

pub(super) fn fit_dispersion_for_gene_detailed(
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

pub(super) fn fit_dispersion_for_gene_detailed_with_weights(
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
            let (dispersion, evaluations) = grid_search_dispersion(DispersionFitInput {
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
                    prepared_cox: None,
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
            let line_search = line_search_dispersion(DispersionFitInput {
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
                dispersion = grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).design(design))"
)]
pub fn fit_dispersion_line_search(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).design(design).prior(prior))"
)]
pub fn fit_dispersion_line_search_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(..).design(design).prior(prior).weights(weights))"
)]
pub fn fit_dispersion_line_search_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples))"
)]
pub fn fit_dispersion_line_search_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).prior(prior))"
)]
pub fn fit_dispersion_line_search_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use line_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).prior(prior).weights(weights))"
)]
pub fn fit_dispersion_line_search_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<DispersionLineSearchOutput, DeseqError> {
    line_search_dispersion(DispersionFitInput {
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

/// One gene's dispersion fit for [`line_search_dispersion`] and
/// [`grid_search_dispersion`].
///
/// Without a design the Cox-Reid adjustment is off whatever
/// `options.use_cox_reid` says. Build with [`DispersionFitInput::new`] and the
/// `design`, `prior` and `weights` setters.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct DispersionFitInput<'a> {
    /// One gene's counts.
    pub counts: &'a [u32],
    /// Fitted means for the same samples.
    pub mu: &'a [f64],
    /// Design for the Cox-Reid adjustment.
    pub design: Option<&'a DesignMatrix>,
    /// Starting dispersion.
    pub initial_dispersion: f64,
    /// Optimizer options.
    pub options: GeneWiseDispersionOptions,
    /// Number of samples, for the default upper dispersion bound.
    pub n_samples: usize,
    /// Optional normal prior on log dispersion (MAP fit).
    pub prior: Option<DispersionPrior>,
    /// Optional observation weights.
    pub weights: Option<&'a [f64]>,
}

impl<'a> DispersionFitInput<'a> {
    /// A gene-wise fit without a design, prior or weights.
    pub fn new(
        counts: &'a [u32],
        mu: &'a [f64],
        initial_dispersion: f64,
        options: GeneWiseDispersionOptions,
        n_samples: usize,
    ) -> Self {
        Self {
            counts,
            mu,
            design: None,
            initial_dispersion,
            options,
            n_samples,
            prior: None,
            weights: None,
        }
    }

    /// Use this design for the Cox-Reid adjustment.
    pub fn design(mut self, design: &'a DesignMatrix) -> Self {
        self.design = Some(design);
        self
    }

    /// Add a normal prior on log dispersion.
    pub fn prior(mut self, prior: DispersionPrior) -> Self {
        self.prior = Some(prior);
        self
    }

    /// Weight the observations.
    pub fn weights(mut self, weights: &'a [f64]) -> Self {
        self.weights = Some(weights);
        self
    }
}

#[derive(Clone, Copy)]
pub(super) struct DispersionObjectiveInput<'a> {
    pub(super) counts: &'a [u32],
    pub(super) mu: &'a [f64],
    pub(super) design: Option<&'a DesignMatrix>,
    pub(super) use_cox_reid: bool,
    pub(super) prior: Option<DispersionPrior>,
    pub(super) weights: Option<&'a [f64]>,
    pub(super) weight_threshold: f64,
    pub(super) prepared_cox: Option<&'a PreparedCoxReidDesign>,
}

/// Fit one dispersion by DESeq2's Armijo line search on log dispersion
/// (`fitDispWrapper`).
pub fn line_search_dispersion(
    input: DispersionFitInput<'_>,
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
    let prepared_cox = prepare_unweighted_cox_reid_design(design, weights, options.use_cox_reid);
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
        prepared_cox: prepared_cox.as_deref(),
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

        // Keep the proposal and Armijo arithmetic in DESeq2's source order.
        // One-ULP differences here can flip near-flat high-dispersion rows.
        let mut proposed_log_alpha = log_alpha + kappa * dlp;
        if !proposed_log_alpha.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion line-search proposal".to_string(),
                index: Some(0),
                value: proposed_log_alpha,
            });
        }
        if proposed_log_alpha < -30.0 {
            kappa = (-30.0 - log_alpha) / dlp;
            proposed_log_alpha = log_alpha + kappa * dlp;
        }
        if proposed_log_alpha > 10.0 {
            kappa = (10.0 - log_alpha) / dlp;
            proposed_log_alpha = log_alpha + kappa * dlp;
        }
        if !kappa.is_finite() || kappa <= 0.0 || !proposed_log_alpha.is_finite() {
            break;
        }

        let theta_kappa = -dispersion_log_posterior_objective(objective, proposed_log_alpha)?;
        let theta_hat_kappa = -lp - kappa * epsilon * dlp * dlp;
        if !theta_hat_kappa.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion line-search Armijo bound".to_string(),
                index: Some(0),
                value: theta_hat_kappa,
            });
        }
        if theta_kappa <= theta_hat_kappa {
            iter_accept += 1;
            log_alpha = proposed_log_alpha;
            // The Armijo test already evaluated this exact point. Reuse its
            // value without changing the accepted-step arithmetic.
            let lp_new = -theta_kappa;
            last_change = checked::try_sub(lp_new, lp)
                .nonfinite_at("dispersion line-search accepted objective change", Some(0))?;
            lp = lp_new;
            if last_change < options.disp_tol {
                break;
            }
            if log_alpha < min_log_alpha {
                break;
            }
            dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
            kappa = (kappa * 1.1).min(options.kappa_0);
            if iter_accept.is_multiple_of(5) {
                kappa /= 2.0;
            }
        } else {
            kappa /= 2.0;
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

#[cfg(test)]
pub(super) fn bounded_log_alpha_proposal(
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
    let unclamped = checked::try_mul(step, direction)
        .nonfinite_at("dispersion line-search proposal step", Some(0))
        .ok()
        .and_then(|movement| {
            checked::try_add(log_alpha, movement)
                .nonfinite_at("dispersion line-search proposal", Some(0))
                .ok()
        })?;
    let clamped = unclamped.clamp(lower, upper);
    if clamped == unclamped {
        return Some((unclamped, step));
    }
    let effective_step = checked::try_sub(clamped, log_alpha)
        .nonfinite_at(
            "dispersion line-search effective proposal movement",
            Some(0),
        )
        .ok()
        .and_then(|movement| {
            checked::try_div(movement, direction)
                .nonfinite_at("dispersion line-search effective proposal step", Some(0))
                .ok()
        })?;
    if effective_step > 0.0 {
        Some((clamped, effective_step))
    } else {
        None
    }
}

#[cfg(test)]
pub(super) fn checked_line_search_armijo_bound(
    lp: f64,
    effective_kappa: f64,
    epsilon: f64,
    dlp: f64,
) -> Result<f64, DeseqError> {
    let scaled_direction = checked::try_mul(effective_kappa, epsilon)
        .nonfinite_at("dispersion line-search Armijo scale", Some(0))?;
    let scaled_slope = checked::try_mul(scaled_direction, dlp)
        .nonfinite_at("dispersion line-search Armijo scale", Some(0))
        .and_then(|scale| {
            checked::try_mul(scale, dlp)
                .nonfinite_at("dispersion line-search Armijo scale", Some(0))
        })?;
    checked::try_sub(-lp, scaled_slope).nonfinite_at("dispersion line-search Armijo bound", Some(0))
}

/// Fit a dispersion for one gene by DESeq2-style two-pass log-alpha grid search.
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).design(design))"
)]
pub fn fit_dispersion_grid(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).design(design).prior(prior))"
)]
pub fn fit_dispersion_grid_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(..).design(design).prior(prior).weights(weights))"
)]
pub fn fit_dispersion_grid_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples))"
)]
pub fn fit_dispersion_grid_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).prior(prior))"
)]
pub fn fit_dispersion_grid_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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
#[deprecated(
    since = "0.6.0",
    note = "use grid_search_dispersion(DispersionFitInput::new(counts, mu, initial_dispersion, options, n_samples).prior(prior).weights(weights))"
)]
pub fn fit_dispersion_grid_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<(f64, usize), DeseqError> {
    grid_search_dispersion(DispersionFitInput {
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

/// Fit one dispersion by DESeq2's two-pass grid search on log dispersion
/// (`fitDispGridWrapper`); returns the dispersion and the number of
/// objective evaluations.
pub fn grid_search_dispersion(input: DispersionFitInput<'_>) -> Result<(f64, usize), DeseqError> {
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
    let prepared_cox = prepare_unweighted_cox_reid_design(design, weights, options.use_cox_reid);
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
        prepared_cox: prepared_cox.as_deref(),
    };
    let min_log = options.min_disp.ln();
    let max_log = max_disp.ln();
    let coarse = linspace(min_log, max_log, options.grid_points)?;
    let (best_log, _) = best_log_alpha(objective, &coarse)?;
    let delta =
        checked::try_sub(coarse[1], coarse[0]).nonfinite_at("dispersion grid step", Some(1))?;
    let fine_lower = checked::try_sub(best_log, delta)
        .nonfinite_at("dispersion fine grid lower bound", Some(0))?;
    let fine_upper = checked::try_add(best_log, delta)
        .nonfinite_at("dispersion fine grid upper bound", Some(0))?;
    let fine = linspace(fine_lower, fine_upper, options.grid_points)?;
    let (best_fine_log, _) = best_log_alpha(objective, &fine)?;
    Ok((
        best_fine_log.exp().clamp(options.min_disp, max_disp),
        options.grid_points * 2,
    ))
}
