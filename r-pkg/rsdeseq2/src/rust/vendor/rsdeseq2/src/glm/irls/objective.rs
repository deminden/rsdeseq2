// Optimizer objectives, gradients and Hessians.

fn beta_log2_value_gradient(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<(f64, Vec<f64>), DeseqError> {
    let (value, gradient, _) = beta_log2_objective_gradient_hessian(input, beta)?;
    Ok((value, gradient))
}

fn beta_log2_objective_gradient_hessian(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<(f64, Vec<f64>, DMatrix<f64>), DeseqError> {
    let p = input.x.ncols();
    if beta.len() != p {
        return Err(invalid_dimensions("optim beta coefficients", p, beta.len()));
    }
    if input.ridge_lambda.len() != p {
        return Err(invalid_dimensions(
            "optim ridge coefficients",
            p,
            input.ridge_lambda.len(),
        ));
    }
    if input.counts.len() != input.x.nrows() || input.nf.len() != input.x.nrows() {
        return Err(invalid_dimensions(
            "optim samples",
            input.x.nrows(),
            input.counts.len().min(input.nf.len()),
        ));
    }
    if let Some(weights) = input.weights
        && weights.len() != input.x.nrows()
    {
        return Err(invalid_dimensions(
            "optim weights",
            input.x.nrows(),
            weights.len(),
        ));
    }
    validate_positive_finite(input.dispersion, "dispersion", 0)?;

    let mut objective = 0.0;
    let mut gradient = vec![0.0; p];
    let mut hessian = DMatrix::zeros(p, p);
    let ln2 = std::f64::consts::LN_2;
    let ln2_squared = ln2 * ln2;
    for sample in 0..input.x.nrows() {
        validate_positive_finite(input.nf[sample], "normalization factor", sample)?;
        let weight = input.weights.map_or(1.0, |weights| weights[sample]);
        validate_nonnegative_finite(weight, "weight", sample)?;
        let mut eta = 0.0;
        for (col, beta_value) in beta.iter().copied().enumerate().take(p) {
            let next_eta = input.x[(sample, col)].mul_add(beta_value, eta);
            if !next_eta.is_finite() {
                return Ok(beta_optim_penalty(p));
            }
            eta = next_eta;
        }
        let mu = input.nf[sample] * 2.0_f64.powf(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Ok(beta_optim_penalty(p));
        }
        let objective_term =
            match beta_log2_nb_objective_term(input.counts[sample], mu, input.dispersion) {
                Some(value) => value,
                None => return Ok(beta_optim_penalty(p)),
            };
        let Some(weighted_objective_term) = checked_product2(weight, objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_objective) = checked_sum2(objective, weighted_objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        objective = next_objective;
        let Some(disp_mu) = checked_product2(input.dispersion, mu) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(one_plus_disp_mu) = checked_sum2(1.0, disp_mu) else {
            return Ok(beta_optim_penalty(p));
        };
        if !one_plus_disp_mu.is_finite() || one_plus_disp_mu <= 0.0 {
            return Ok(beta_optim_penalty(p));
        }
        let inv_one_plus_disp_mu = one_plus_disp_mu.recip();
        let count_residual = f64::from(input.counts[sample]) - mu;
        if !count_residual.is_finite() {
            return Ok(beta_optim_penalty(p));
        }
        let Some(sample_score) =
            checked_product4(weight, ln2, count_residual, inv_one_plus_disp_mu)
        else {
            return Ok(beta_optim_penalty(p));
        };
        for (col, gradient_value) in gradient.iter_mut().enumerate().take(p) {
            let Some(term) = checked_product2(-input.x[(sample, col)], sample_score) else {
                return Ok(beta_optim_penalty(p));
            };
            let Some(next_gradient) = checked_sum2(*gradient_value, term) else {
                return Ok(beta_optim_penalty(p));
            };
            *gradient_value = next_gradient;
        }
        let Some(disp_count) = checked_product2(input.dispersion, f64::from(input.counts[sample]))
        else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(one_plus_disp_count) = checked_sum2(1.0, disp_count) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(sample_hessian_weight) = checked_product6(
            weight,
            ln2_squared,
            mu,
            one_plus_disp_count,
            inv_one_plus_disp_mu,
            inv_one_plus_disp_mu,
        ) else {
            return Ok(beta_optim_penalty(p));
        };
        for row in 0..p {
            for col in 0..p {
                let Some(term) = checked_product3(
                    input.x[(sample, row)],
                    sample_hessian_weight,
                    input.x[(sample, col)],
                ) else {
                    return Ok(beta_optim_penalty(p));
                };
                let Some(next_hessian) = checked_sum2(hessian[(row, col)], term) else {
                    return Ok(beta_optim_penalty(p));
                };
                hessian[(row, col)] = next_hessian;
            }
        }
    }

    if !objective.is_finite() {
        return Ok(beta_optim_penalty(p));
    }
    for col in 0..p {
        validate_nonnegative_finite(input.ridge_lambda[col], "ridge lambda", col)?;
        let Some(ridge_log2) = checked_product2(input.ridge_lambda[col], ln2_squared) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(objective_term) = checked_product4(0.5, ridge_log2, beta[col], beta[col]) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(gradient_term) = checked_product2(ridge_log2, beta[col]) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_objective) = checked_sum2(objective, objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_gradient) = checked_sum2(gradient[col], gradient_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_hessian) = checked_sum2(hessian[(col, col)], ridge_log2) else {
            return Ok(beta_optim_penalty(p));
        };
        objective = next_objective;
        gradient[col] = next_gradient;
        hessian[(col, col)] = next_hessian;
    }
    Ok((objective, gradient, hessian))
}

fn beta_log2_nb_objective_term(count: u32, mu: f64, dispersion: f64) -> Option<f64> {
    validate_positive_finite(mu, "optim fallback mean", 0).ok()?;
    validate_positive_finite(dispersion, "dispersion", 0).ok()?;
    let count = f64::from(count);
    let size = dispersion.recip();
    if !size.is_finite() || size <= 0.0 {
        return None;
    }
    let size_plus_mu = checked_sum2(size, mu)?;
    let count_plus_size = checked_sum2(count, size)?;
    let first = checked_product2(count_plus_size, size_plus_mu.ln())?;
    let second = checked_product2(size, size.ln())?;
    let third = checked_product2(count, mu.ln())?;
    checked_sum2(checked_sum2(first, -second)?, -third)
}

fn beta_optim_penalty(p: usize) -> (f64, Vec<f64>, DMatrix<f64>) {
    (1.0e300, vec![0.0; p], DMatrix::identity(p, p))
}

fn beta_log2_objective(input: &BetaOptimInput<'_>, beta: &[f64]) -> Result<f64, DeseqError> {
    beta_log2_objective_gradient_hessian(input, beta).map(|(objective, _, _)| objective)
}

fn beta_log2_r_objective(input: &BetaOptimInput<'_>, beta: &[f64]) -> Result<f64, DeseqError> {
    let mut ll = 0.0;
    for sample in 0..input.x.nrows() {
        let mut eta = 0.0;
        for (col, b) in beta.iter().copied().enumerate() {
            eta = input.x[(sample, col)].mul_add(b, eta)
        }
        let mu = input.nf[sample] * 2.0_f64.powf(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Ok(1e300);
        }
        ll += input.weights.map_or(1.0, |w| w[sample])
            * nbinom_log(input.counts[sample], mu, input.dispersion)
    }
    let l2 = std::f64::consts::LN_2.powi(2);
    let mut lp = 0.0;
    for (col, b) in beta.iter().copied().enumerate() {
        let precision = input.ridge_lambda[col] * l2;
        if precision > 0.0 {
            lp += normal_log(b, precision);
        }
    }
    Ok(-(ll + lp))
}

#[cfg(test)]
fn beta_log2_numeric_gradient(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
    lower: f64,
    upper: f64,
) -> Result<Vec<f64>, DeseqError> {
    const R_OPTIM_NDEPS: f64 = 1.0e-3;
    let mut gradient = Vec::with_capacity(beta.len());
    for col in 0..beta.len() {
        let forward = (beta[col] + R_OPTIM_NDEPS).min(upper);
        let forward_eps = forward - beta[col];
        let backward = (beta[col] - R_OPTIM_NDEPS).max(lower);
        let backward_eps = beta[col] - backward;
        let denominator = forward_eps + backward_eps;
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(DeseqError::UnsupportedFeature {
                feature: "optim fallback finite-difference gradient at degenerate bounds"
                    .to_string(),
            });
        }
        let mut forward_beta = beta.to_vec();
        forward_beta[col] = forward;
        let mut backward_beta = beta.to_vec();
        backward_beta[col] = backward;
        let forward_value = beta_log2_objective(input, &forward_beta)?;
        let backward_value = beta_log2_objective(input, &backward_beta)?;
        let Some(difference) = checked_sum2(forward_value, -backward_value) else {
            return Ok(vec![0.0; beta.len()]);
        };
        let Some(value) = checked_div2(difference, denominator) else {
            return Ok(vec![0.0; beta.len()]);
        };
        gradient.push(value);
    }
    Ok(gradient)
}
