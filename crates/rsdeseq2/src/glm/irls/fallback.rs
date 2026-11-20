use super::*;
use crate::glm::nb::nbinom_log_likelihood_weighted;
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
    let lower = -input.options.max_beta_abs;
    let upper = input.options.max_beta_abs;
    let control = optim_fallback_control(p, input.options.optim_maxit, input.options.optim_tol);
    let mut objective_scratch = BetaObjectiveScratch::new(p)?;
    let mut stencil = try_filled_vec(p, 0.0, "IRLS refit stencil coefficients")?;
    let mut inference_mu = try_filled_vec(n, 0.0, "IRLS refit inference means")?;
    let mut working_weights = try_alloc_matrix(n, 1, "IRLS refit working weights")?;
    let mut scaled_factors = try_alloc_matrix(
        if input.rows.offset_scales.is_some() {
            n
        } else {
            0
        },
        1,
        "IRLS refit scaled offsets",
    )?;
    let mut covariance_scratch = CovarianceScratch::default();
    for gene in rows.iter().copied() {
        let source_gene = input.rows.source(gene);
        let counts_row = input.rows.input.counts.row(source_gene)?;
        let nf = input.rows.factors(gene, &mut scaled_factors)?;
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
        let mut objective = objective_scratch.for_row(beta_input);
        let start_objective = if input.options.r_optim_compat {
            beta_log2_r_objective(&beta_input, &beta_start)?
        } else {
            objective.value(&beta_start)?
        };
        let output = minimize_beta_log2(
            &mut objective,
            &beta_start,
            (lower, upper),
            &control,
            input.options.r_optim_compat,
            &mut stencil,
        )?;
        let final_gradient = objective.gradient(&output.parameters)?;
        let final_gradient_norm =
            projected_gradient_norm(&output.parameters, final_gradient, lower, upper)
                .unwrap_or(f64::NAN);

        let mu_row = &mut mu_values[gene * n..(gene + 1) * n];
        fitted_mu_log2_unfloored_into(input.x, &output.parameters, nf, mu_row)?;
        for (value, mu) in inference_mu.iter_mut().zip(mu_row.iter().copied()) {
            *value = mu.max(input.options.min_mu);
        }
        working_weights_into(&inference_mu, dispersion, weight_row, &mut working_weights)?;
        // Fallback preserves the IRLS leverage for Cook's distance.
        let Some(beta_covariance) = covariance_and_hat_diagonal(
            &mut covariance_scratch,
            input.x,
            &working_weights,
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
            &inference_mu,
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

#[cfg(test)]
pub(super) fn minimize_beta_log2_lbfgsb(
    input: BetaOptimInput<'_>,
    start: &[f64],
    lower: f64,
    upper: f64,
    maxit: usize,
    gradient_tol: f64,
    r_optim_compat: bool,
) -> Result<BoundedOptimizationOutput, DeseqError> {
    let mut scratch = BetaObjectiveScratch::new(input.x.ncols())?;
    let mut stencil = vec![0.0; start.len()];
    minimize_beta_log2(
        &mut scratch.for_row(input),
        start,
        (lower, upper),
        &optim_fallback_control(start.len(), maxit, gradient_tol),
        r_optim_compat,
        &mut stencil,
    )
}

fn minimize_beta_log2(
    evaluation: &mut BetaObjective<'_, '_>,
    start: &[f64],
    (lower, upper): (f64, f64),
    control: &OptimControl,
    r_optim_compat: bool,
    stencil: &mut [f64],
) -> Result<BoundedOptimizationOutput, DeseqError> {
    let input = evaluation.input;
    let mut evaluator = RefCell::new(evaluation);
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
        evaluator.borrow_mut().value(beta)
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
        evaluator.borrow_mut().value(&parameters)?
    };
    let bounds = Bounds::new(vec![lower; parameters.len()], vec![upper; parameters.len()])
        .map_err(|error| DeseqError::UnsupportedFeature {
            feature: format!("optim fallback L-BFGS-B bounds failed: {error}"),
        })?;
    let solution = if r_optim_compat {
        let numeric_fit = if crate::numeric_profile() == crate::NumericProfile::Portable {
            let mut gradient = |beta: &[f64]| match beta_log2_r_centered_gradient_into(
                &input,
                beta,
                lower,
                upper,
                &control.ndeps,
                stencil,
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
            |beta: &[f64]| evaluator.borrow_mut().value(beta).unwrap_or(1.0e300);
        let mut analytic_gradient = |beta: &[f64]| {
            evaluator
                .borrow_mut()
                .gradient(beta)
                .map(<[f64]>::to_vec)
                .unwrap_or_else(|_| vec![0.0; beta.len()])
        };
        if let Ok(mut analytic) = optim_lbfgsb_with_gradient(
            parameters,
            bounds,
            &mut analytic_objective,
            &mut analytic_gradient,
            control.clone(),
        ) {
            let numeric_reduced = evaluator.borrow_mut().value(&chosen.par)?;
            let analytic_reduced = evaluator.borrow_mut().value(&analytic.par)?;
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
        let mut gradient = |beta: &[f64]| match evaluator.borrow_mut().gradient(beta) {
            Ok(gradient) => gradient.to_vec(),
            Err(error) => {
                deferred_error.borrow_mut().get_or_insert(error);
                vec![0.0; beta.len()]
            }
        };
        optim_lbfgsb_with_gradient(
            parameters,
            bounds,
            &mut objective,
            &mut gradient,
            control.clone(),
        )
    }
    .map_err(|error| DeseqError::UnsupportedFeature {
        feature: format!("optim fallback L-BFGS-B failed: {error}"),
    })?;
    if let Some(error) = deferred_error.into_inner() {
        return Err(error);
    }
    let final_gradient = evaluator.get_mut().gradient(&solution.par)?;
    let final_gradient_norm = projected_gradient_norm(&solution.par, final_gradient, lower, upper);
    let converged_by_gradient = final_gradient_norm
        .is_some_and(|norm| norm <= control.pgtol.max(1.0e-4) && solution.value.is_finite());
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
