use super::*;
use crate::numeric::math::ProfileMath;
// Optimizer objectives, gradients and Hessians.

pub(super) struct BetaObjectiveScratch {
    last_beta: Vec<f64>,
    gradient: Vec<f64>,
    hessian: DenseMatrix,
}

impl BetaObjectiveScratch {
    pub(super) fn new(coefficients: usize) -> Result<Self, DeseqError> {
        Ok(Self {
            last_beta: try_filled_vec(coefficients, 0.0, "IRLS refit objective coefficients")?,
            gradient: try_filled_vec(coefficients, 0.0, "IRLS refit gradient")?,
            hessian: try_zero_matrix(coefficients, coefficients, "IRLS refit Hessian")?,
        })
    }

    pub(super) fn for_row<'row>(&mut self, input: BetaOptimInput<'row>) -> BetaObjective<'row, '_> {
        BetaObjective {
            input,
            scratch: self,
            cached_value: None,
        }
    }
}

/// The objective and gradient callbacks share an evaluation for a single gene.
pub(super) struct BetaObjective<'row, 'scratch> {
    pub(super) input: BetaOptimInput<'row>,
    scratch: &'scratch mut BetaObjectiveScratch,
    cached_value: Option<f64>,
}

impl BetaObjective<'_, '_> {
    pub(super) fn value(&mut self, beta: &[f64]) -> Result<f64, DeseqError> {
        if let Some(value) = self.cached_value
            && beta.len() == self.scratch.last_beta.len()
            && beta
                .iter()
                .zip(&self.scratch.last_beta)
                .all(|(left, right)| left.to_bits() == right.to_bits())
        {
            return Ok(value);
        }
        self.cached_value = None;
        let value = evaluate_beta_log2(
            &self.input,
            beta,
            &mut self.scratch.gradient,
            &mut self.scratch.hessian,
        )?;
        self.scratch.last_beta.copy_from_slice(beta);
        self.cached_value = Some(value);
        Ok(value)
    }

    pub(super) fn gradient(&mut self, beta: &[f64]) -> Result<&[f64], DeseqError> {
        self.value(beta)?;
        Ok(&self.scratch.gradient)
    }
}

#[cfg(test)]
pub(super) fn beta_log2_objective_gradient_hessian(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<(f64, Vec<f64>, DenseMatrix), DeseqError> {
    let mut scratch = BetaObjectiveScratch::new(input.x.ncols())?;
    let value = scratch.for_row(*input).value(beta)?;
    Ok((value, scratch.gradient, scratch.hessian))
}

fn evaluate_beta_log2(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
    gradient: &mut [f64],
    hessian: &mut DenseMatrix,
) -> Result<f64, DeseqError> {
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
    gradient.fill(0.0);
    hessian.fill(0.0);
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
                return Ok(beta_optim_penalty(gradient, hessian));
            }
            eta = next_eta;
        }
        let mu = input.nf[sample] * 2.0_f64.powf_p(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Ok(beta_optim_penalty(gradient, hessian));
        }
        let objective_term =
            match beta_log2_nb_objective_term(input.counts[sample], mu, input.dispersion) {
                Some(value) => value,
                None => return Ok(beta_optim_penalty(gradient, hessian)),
            };
        let Some(weighted_objective_term) = checked::mul(weight, objective_term) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(next_objective) = checked::add(objective, weighted_objective_term) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        objective = next_objective;
        let Some(disp_mu) = checked::mul(input.dispersion, mu) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(one_plus_disp_mu) = checked::add(1.0, disp_mu) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        if !one_plus_disp_mu.is_finite() || one_plus_disp_mu <= 0.0 {
            return Ok(beta_optim_penalty(gradient, hessian));
        }
        let inv_one_plus_disp_mu = one_plus_disp_mu.recip();
        let count_residual = f64::from(input.counts[sample]) - mu;
        if !count_residual.is_finite() {
            return Ok(beta_optim_penalty(gradient, hessian));
        }
        let Some(sample_score) =
            checked_product4(weight, ln2, count_residual, inv_one_plus_disp_mu)
        else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        for (col, gradient_value) in gradient.iter_mut().enumerate().take(p) {
            let Some(term) = checked::mul(-input.x[(sample, col)], sample_score) else {
                return Ok(beta_optim_penalty(gradient, hessian));
            };
            let Some(next_gradient) = checked::add(*gradient_value, term) else {
                return Ok(beta_optim_penalty(gradient, hessian));
            };
            *gradient_value = next_gradient;
        }
        let Some(disp_count) = checked::mul(input.dispersion, f64::from(input.counts[sample]))
        else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(one_plus_disp_count) = checked::add(1.0, disp_count) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(sample_hessian_weight) = checked_product6(
            weight,
            ln2_squared,
            mu,
            one_plus_disp_count,
            inv_one_plus_disp_mu,
            inv_one_plus_disp_mu,
        ) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        for row in 0..p {
            for col in 0..p {
                let Some(term) = checked_product3(
                    input.x[(sample, row)],
                    sample_hessian_weight,
                    input.x[(sample, col)],
                ) else {
                    return Ok(beta_optim_penalty(gradient, hessian));
                };
                let Some(next_hessian) = checked::add(hessian[(row, col)], term) else {
                    return Ok(beta_optim_penalty(gradient, hessian));
                };
                hessian[(row, col)] = next_hessian;
            }
        }
    }

    if !objective.is_finite() {
        return Ok(beta_optim_penalty(gradient, hessian));
    }
    for col in 0..p {
        validate_nonnegative_finite(input.ridge_lambda[col], "ridge lambda", col)?;
        let Some(ridge_log2) = checked::mul(input.ridge_lambda[col], ln2_squared) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(objective_term) = checked_product4(0.5, ridge_log2, beta[col], beta[col]) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(gradient_term) = checked::mul(ridge_log2, beta[col]) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(next_objective) = checked::add(objective, objective_term) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(next_gradient) = checked::add(gradient[col], gradient_term) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        let Some(next_hessian) = checked::add(hessian[(col, col)], ridge_log2) else {
            return Ok(beta_optim_penalty(gradient, hessian));
        };
        objective = next_objective;
        gradient[col] = next_gradient;
        hessian[(col, col)] = next_hessian;
    }
    Ok(objective)
}

pub(super) fn beta_log2_nb_objective_term(count: u32, mu: f64, dispersion: f64) -> Option<f64> {
    validate_positive_finite(mu, "optim fallback mean", 0).ok()?;
    validate_positive_finite(dispersion, "dispersion", 0).ok()?;
    let count = f64::from(count);
    let size = dispersion.recip();
    if !size.is_finite() || size <= 0.0 {
        return None;
    }
    let size_plus_mu = checked::add(size, mu)?;
    let count_plus_size = checked::add(count, size)?;
    let first = checked::mul(count_plus_size, size_plus_mu.ln_p())?;
    let second = checked::mul(size, size.ln_p())?;
    let third = checked::mul(count, mu.ln_p())?;
    checked::add(checked::add(first, -second)?, -third)
}

fn beta_optim_penalty(gradient: &mut [f64], hessian: &mut DenseMatrix) -> f64 {
    gradient.fill(0.0);
    hessian.fill(0.0);
    for diagonal in 0..hessian.nrows() {
        hessian[(diagonal, diagonal)] = 1.0;
    }
    1.0e300
}

#[cfg(test)]
pub(super) fn beta_log2_objective(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<f64, DeseqError> {
    beta_log2_objective_gradient_hessian(input, beta).map(|(objective, _, _)| objective)
}

pub(super) fn beta_log2_r_objective(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<f64, DeseqError> {
    let mut ll = 0.0;
    for sample in 0..input.x.nrows() {
        let mut eta = 0.0;
        for (col, b) in beta.iter().copied().enumerate() {
            eta = input.x[(sample, col)].mul_add(b, eta)
        }
        let mu = input.nf[sample] * 2.0_f64.powf_p(eta);
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
pub(super) fn beta_log2_r_centered_gradient(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
    lower: f64,
    upper: f64,
    ndeps: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    let mut candidate = vec![0.0; beta.len()];
    beta_log2_r_centered_gradient_into(input, beta, lower, upper, ndeps, &mut candidate)
}

pub(super) fn beta_log2_r_centered_gradient_into(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
    lower: f64,
    upper: f64,
    ndeps: &[f64],
    candidate: &mut [f64],
) -> Result<Vec<f64>, DeseqError> {
    let ln2 = std::f64::consts::LN_2;
    let size = input.dispersion.recip();
    candidate.copy_from_slice(beta);
    let mut gradient = try_filled_vec(beta.len(), 0.0, "IRLS refit centered gradient")?;
    for (col, value) in gradient.iter_mut().enumerate() {
        let forward_step = ndeps[col].min(upper - beta[col]);
        let backward_step = ndeps[col].min(beta[col] - lower);
        let denominator = forward_step + backward_step;
        if lower == upper {
            continue;
        }
        let backward = beta[col] - backward_step;
        let forward = beta[col] + forward_step;
        if denominator <= 0.0 || (forward == beta[col] && backward == beta[col]) {
            return Err(DeseqError::UnsupportedFeature {
                feature: format!("optim fallback stencil cannot perturb coordinate {col}"),
            });
        }
        candidate[col] = backward;
        let coordinate_change = (forward - backward) * ln2;
        let precision = input.ridge_lambda[col] * ln2.powi(2);
        let ridge_change = (0.5 * precision) * (forward - backward) * (forward + backward);
        let mut roundoff = 0.0;
        let centered = (0..input.x.nrows())
            .try_fold(0.0, |difference, sample| {
                let term = beta_log2_sample_loss_change(
                    input,
                    candidate,
                    col,
                    sample,
                    size,
                    coordinate_change,
                )?;
                let weight = input.weights.map_or(1.0, |weights| weights[sample]);
                // The signed sample changes can nearly cancel at an optimum.
                let increment = weight * term - roundoff;
                let next = difference + increment;
                roundoff = (next - difference) - increment;
                Some(next)
            })
            .map(|difference| difference + ridge_change)
            .filter(|difference| difference.is_finite());
        let difference = match centered {
            Some(difference) => difference,
            None => {
                // Keep the R objective's behaviour when a stencil leaves the
                // ordinary floating-point range.
                let before = beta_log2_r_objective(input, candidate)?;
                candidate[col] = forward;
                beta_log2_r_objective(input, candidate)? - before
            }
        };
        candidate[col] = beta[col];
        *value = difference / denominator;
    }
    Ok(gradient)
}

fn beta_log2_sample_loss_change(
    input: &BetaOptimInput<'_>,
    backward: &[f64],
    col: usize,
    sample: usize,
    size: f64,
    coordinate_change: f64,
) -> Option<f64> {
    let eta = backward
        .iter()
        .copied()
        .enumerate()
        .fold(0.0, |eta, (j, coefficient)| {
            input.x[(sample, j)].mul_add(coefficient, eta)
        });
    let delta = input.x[(sample, col)] * coordinate_change;
    // Leave room for exp2 and the offset product at the other end of the
    // stencil. The direct R objective handles wider steps and extreme means.
    if eta.abs() > 1020.0 || delta.abs() > 1.0 {
        return None;
    }
    let mu = input.nf[sample] * 2.0_f64.powf_p(eta);
    let total = mu + size;
    if mu < 3.0 * f64::MIN_POSITIVE || mu > f64::MAX / 3.0 || !total.is_finite() {
        return None;
    }
    let share = mu.min(size) / total;
    if !share.is_normal() {
        return None;
    }
    let count = f64::from(input.counts[sample]);
    // The means satisfy mu_forward = mu * exp(delta). Evaluate the change in
    // log(size + mu) after factoring out the larger contribution. Terms that
    // depend only on the count or dispersion cancel exactly.
    let term = if mu <= size {
        let shifted = share * delta.exp_m1_p();
        (count + size).mul_add(shifted.ln_1p_p(), -count * delta)
    } else {
        let shifted = share * (-delta).exp_m1_p();
        (count + size).mul_add(shifted.ln_1p_p(), size * delta)
    };
    Some(term)
}

#[cfg(test)]
pub(super) fn beta_log2_numeric_gradient(
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
        let Some(difference) = checked::add(forward_value, -backward_value) else {
            return Ok(vec![0.0; beta.len()]);
        };
        let Some(value) = checked::div(difference, denominator) else {
            return Ok(vec![0.0; beta.len()]);
        };
        gradient.push(value);
    }
    Ok(gradient)
}
