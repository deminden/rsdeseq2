use super::*;
// Bounded-optimizer fallback, starting values and acceptance checks.

#[allow(clippy::too_many_arguments)]
pub(super) fn refit_optim_fallback_rows<const SELECTED: bool>(
    rows: &[usize],
    beta_values: &mut [f64],
    beta_var_values: &mut [f64],
    beta_optim_start_values: &mut [f64],
    beta_covariance_values: &mut [f64],
    mu_values: &mut [f64],
    beta_converged: &mut [bool],
    optim_log_like: &mut [Option<f64>],
    beta_optim_iter: &mut [f64],
    beta_optim_start_objective: &mut [f64],
    beta_optim_objective: &mut [f64],
    beta_optim_gradient_norm: &mut [f64],
    input: OptimFallbackInput<'_, SELECTED>,
) -> Result<(), DeseqError> {
    let p = input.x.ncols();
    let n = input.x.nrows();
    let mut covariance_scratch = CovarianceScratch::default();
    for gene in rows.iter().copied() {
        let source_gene = input.rows.source(gene);
        let counts_row = input.rows.input.counts.row(source_gene)?;
        let nf = input.rows.input.offsets.row(source_gene)?;
        let weight_row = input.rows.weights(source_gene)?;
        let dispersion = input.rows.input.dispersions[source_gene];
        let beta_start = optim_start_beta_log2(
            &beta_values[gene * p..(gene + 1) * p],
            input.x,
            counts_row,
            nf,
            input.options.max_beta_abs,
        )?;
        let beta_input = BetaOptimInput {
            x: input.x,
            counts: counts_row,
            nf,
            dispersion,
            weights: weight_row,
            ridge_lambda: input.ridge_lambda,
        };
        let start_objective = if input.options.r_optim_compat {
            beta_log2_r_objective(&beta_input, &beta_start)?
        } else {
            beta_log2_objective(&beta_input, &beta_start)?
        };
        let output = optimize_beta_log2(beta_input, &beta_start, input.options)?;
        let (_, final_gradient, _) =
            beta_log2_objective_gradient_hessian(&beta_input, &output.parameters)?;
        let final_gradient_norm = projected_gradient_norm(
            &output.parameters,
            &final_gradient,
            -input.options.max_beta_abs,
            input.options.max_beta_abs,
        )
        .unwrap_or(f64::NAN);

        let mu_unfloored = fitted_mu_log2_unfloored(input.x, &output.parameters, nf)?;
        let mu_for_inference = mu_unfloored
            .iter()
            .copied()
            .map(|value| value.max(input.options.min_mu))
            .collect::<Vec<_>>();
        let w = working_weights(&mu_for_inference, dispersion, weight_row)?;
        // Fallback preserves the IRLS leverage for Cook's distance.
        let Some(beta_covariance) = covariance_and_hat_diagonal(
            &mut covariance_scratch,
            input.x,
            &w,
            input.ridge_lambda,
            None,
            input.max_design,
            !beta_covariance_values.is_empty(),
        )?
        else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "optim fallback covariance is singular".to_string(),
            });
        };

        let beta_row = gene * p..(gene + 1) * p;
        beta_values[beta_row.clone()].copy_from_slice(&output.parameters);
        if !beta_optim_start_values.is_empty() {
            beta_optim_start_values[beta_row].copy_from_slice(&beta_start);
        }
        if let FitCovariance::Full(covariance) = beta_covariance {
            let start_index = gene * p * p;
            let destination = if beta_covariance_values.is_empty() {
                None
            } else {
                for col in 0..p {
                    beta_var_values[gene * p + col] = checked_log2_standard_error(
                        covariance[(col, col)],
                        gene * p + col,
                        "optim fallback beta standard error",
                    )?;
                }
                Some(&mut beta_covariance_values[start_index..start_index + p * p])
            };
            retain_covariance_on_log2_scale(
                covariance,
                destination,
                start_index,
                "optim fallback beta covariance log2 scale",
            )?;
        }
        for (sample, value) in mu_unfloored.iter().copied().enumerate() {
            mu_values[gene * n + sample] = value;
        }
        // DESeq2 updates beta, SE, mu, and log-likelihood after L-BFGS-B
        // fallback, but preserves the original IRLS hat diagonals returned by
        // fitBetaWrapper. Cook's replacement relies on that pre-optim leverage.
        beta_converged[gene] = output.converged;
        beta_optim_iter[gene] = output.iterations as f64;
        beta_optim_start_objective[gene] = start_objective;
        beta_optim_objective[gene] = output.value;
        beta_optim_gradient_norm[gene] = final_gradient_norm;
        optim_log_like[gene] = Some(nbinom_log_likelihood_weighted(
            counts_row,
            &mu_for_inference,
            dispersion,
            weight_row,
        )?);
    }
    Ok(())
}

pub(super) fn optim_start_beta_log2(
    current_beta_log2: &[f64],
    x: &DMatrix<f64>,
    counts: &[u32],
    nf: &[f64],
    bound: f64,
) -> Result<Vec<f64>, DeseqError> {
    if current_beta_log2
        .iter()
        .copied()
        .all(|value| value.is_finite() && value.abs() <= bound)
    {
        return Ok(current_beta_log2
            .iter()
            .copied()
            .map(|value| value.clamp(-bound, bound))
            .collect());
    }

    let y = counts.iter().copied().map(f64::from).collect::<Vec<_>>();
    Ok(initial_beta(x, &y, nf)?
        .iter()
        .copied()
        .map(|value| value.clamp(-bound, bound))
        .collect())
}

#[derive(Clone, Copy)]
pub(super) struct BetaOptimInput<'a> {
    pub(super) x: &'a DMatrix<f64>,
    pub(super) counts: &'a [u32],
    pub(super) nf: &'a [f64],
    pub(super) dispersion: f64,
    pub(super) weights: Option<&'a [f64]>,
    pub(super) ridge_lambda: &'a [f64],
}

pub(super) fn optimize_beta_log2(
    input: BetaOptimInput<'_>,
    start: &[f64],
    options: &IrlsOptions,
) -> Result<BoundedOptimizationOutput, DeseqError> {
    minimize_beta_log2_lbfgsb(
        input,
        start,
        -options.max_beta_abs,
        options.max_beta_abs,
        options.optim_maxit,
        options.optim_tol,
        options.r_optim_compat,
    )
}

pub(super) fn minimize_beta_log2_lbfgsb(
    input: BetaOptimInput<'_>,
    start: &[f64],
    lower: f64,
    upper: f64,
    maxit: usize,
    gradient_tol: f64,
    r_optim_compat: bool,
) -> Result<BoundedOptimizationOutput, DeseqError> {
    let parameters = start
        .iter()
        .copied()
        .map(|value| value.clamp(lower, upper))
        .collect::<Vec<_>>();
    let start_parameters = parameters.clone();
    let deferred_error = RefCell::new(None);
    let mut objective = |beta: &[f64]| match if r_optim_compat {
        beta_log2_r_objective(&input, beta)
    } else {
        beta_log2_objective(&input, beta)
    } {
        Ok(value) => value,
        Err(error) => {
            deferred_error.borrow_mut().get_or_insert(error);
            1.0e300
        }
    };
    let start_value = if r_optim_compat {
        beta_log2_r_objective(&input, &parameters)?
    } else {
        beta_log2_objective(&input, &parameters)?
    };
    let bounds = Bounds::new(vec![lower; parameters.len()], vec![upper; parameters.len()])
        .map_err(|error| DeseqError::UnsupportedFeature {
            feature: format!("optim fallback L-BFGS-B bounds failed: {error}"),
        })?;
    let control = optim_fallback_control(parameters.len(), maxit, gradient_tol);
    let solution = if r_optim_compat {
        let numeric_fit = if crate::numeric_profile() == crate::NumericProfile::Portable {
            let mut gradient = |beta: &[f64]| match beta_log2_r_centered_gradient(
                &input,
                beta,
                lower,
                upper,
                &control.ndeps,
            ) {
                Ok(gradient) => gradient,
                Err(error) => {
                    deferred_error.borrow_mut().get_or_insert(error);
                    vec![0.0; beta.len()]
                }
            };
            optim_lbfgsb_with_gradient(
                parameters.clone(),
                bounds.clone(),
                &mut objective,
                &mut gradient,
                control.clone(),
            )
        } else {
            optim_lbfgsb(
                parameters.clone(),
                bounds.clone(),
                &mut objective,
                control.clone(),
            )
        };
        let numeric = match numeric_fit {
            Ok(solution) => solution,
            Err(error) => {
                return Err(DeseqError::UnsupportedFeature {
                    feature: format!("optim fallback L-BFGS-B failed: {error}"),
                });
            }
        };
        let mut chosen = numeric;
        let mut analytic_objective =
            |beta: &[f64]| beta_log2_objective(&input, beta).unwrap_or(1.0e300);
        let mut analytic_gradient = |beta: &[f64]| {
            beta_log2_value_gradient(&input, beta)
                .map(|(_, gradient)| gradient)
                .unwrap_or_else(|_| vec![0.0; beta.len()])
        };
        if let Ok(mut analytic) = optim_lbfgsb_with_gradient(
            parameters,
            bounds,
            &mut analytic_objective,
            &mut analytic_gradient,
            control.clone(),
        ) {
            let numeric_reduced = beta_log2_objective(&input, &chosen.par)?;
            let analytic_reduced = beta_log2_objective(&input, &analytic.par)?;
            if analytic.is_success()
                && !r_optim_candidate_within_analytic_tolerance(
                    numeric_reduced,
                    analytic_reduced,
                    control.factr,
                )
                && r_optim_candidates_materially_disagree(
                    &chosen.par,
                    &analytic.par,
                    &control.ndeps,
                    &control.parscale,
                )
            {
                analytic.value = beta_log2_r_objective(&input, &analytic.par)?;
                chosen = analytic;
            }
        }
        Ok(chosen)
    } else {
        let mut gradient = |beta: &[f64]| match beta_log2_value_gradient(&input, beta) {
            Ok((_, gradient)) => gradient,
            Err(error) => {
                deferred_error.borrow_mut().get_or_insert(error);
                vec![0.0; beta.len()]
            }
        };
        optim_lbfgsb_with_gradient(parameters, bounds, &mut objective, &mut gradient, control)
    }
    .map_err(|error| DeseqError::UnsupportedFeature {
        feature: format!("optim fallback L-BFGS-B failed: {error}"),
    })?;
    if let Some(error) = deferred_error.into_inner() {
        return Err(error);
    }
    let (_, final_gradient) = beta_log2_value_gradient(&input, &solution.par)?;
    let final_gradient_norm = projected_gradient_norm(&solution.par, &final_gradient, lower, upper);
    let converged_by_gradient = final_gradient_norm
        .is_some_and(|norm| norm <= gradient_tol.max(1.0e-4) && solution.value.is_finite());
    let converged = solution.is_success() || converged_by_gradient;
    let iterations = solution.gradient_count();
    if !converged && final_gradient_norm.is_some_and(|norm| norm > 0.1) {
        return Ok(BoundedOptimizationOutput {
            parameters: start_parameters,
            value: start_value,
            converged: false,
            iterations,
        });
    }
    Ok(BoundedOptimizationOutput {
        parameters: solution.par,
        value: solution.value,
        converged,
        iterations,
    })
}

pub(super) fn optim_fallback_control(
    dimension: usize,
    maxit: usize,
    gradient_tol: f64,
) -> OptimControl {
    let mut control = OptimControl::default_for_dimension(dimension);
    control.maxit = maxit;
    control.factr = 1.0e7;
    control.pgtol = gradient_tol;
    control.lmm = 5;
    control
}

pub(super) fn r_optim_candidate_within_analytic_tolerance(
    compatible_objective: f64,
    analytic_objective: f64,
    factr: f64,
) -> bool {
    if !compatible_objective.is_finite()
        || !analytic_objective.is_finite()
        || !factr.is_finite()
        || factr < 0.0
    {
        return false;
    }
    // `factr * epsilon` is the L-BFGS-B relative-reduction tolerance. Allow a
    // tenfold band for the deliberately finite-difference-compatible endpoint;
    // beyond it the callback path and analytic solution disagree materially.
    let scale = compatible_objective
        .abs()
        .max(analytic_objective.abs())
        .max(1.0);
    let tolerance = 10.0 * factr * f64::EPSILON * scale;
    compatible_objective <= analytic_objective + tolerance
}

pub(super) fn r_optim_candidates_materially_disagree(
    compatible: &[f64],
    analytic: &[f64],
    ndeps: &[f64],
    parscale: &[f64],
) -> bool {
    if compatible.len() != analytic.len()
        || compatible.len() != ndeps.len()
        || compatible.len() != parscale.len()
    {
        return true;
    }
    compatible
        .iter()
        .zip(analytic)
        .zip(ndeps.iter().zip(parscale))
        .any(|((&r_value, &analytic_value), (&step, &scale))| {
            let resolution = step * scale;
            !r_value.is_finite()
                || !analytic_value.is_finite()
                || !resolution.is_finite()
                || resolution <= 0.0
                || (r_value - analytic_value).abs() > 10.0 * resolution
        })
}
