use super::*;
use crate::numeric::math::ProfileMath;
// Per-gene IRLS iterations and scratch reuse.

#[derive(Default)]
pub(super) struct IrlsScratch {
    counts: Vec<f64>,
    log_offsets: Vec<f64>,
    fitted_means: Vec<f64>,
    weights: Vec<f64>,
    response: Vec<f64>,
    covariance: Vec<f64>,
}

pub(super) fn fit_fixed_dispersion_irls_row<const SELECTED: bool>(
    scratch: &mut IrlsScratch,
    input: IrlsRowInput<'_, SELECTED>,
    output: IrlsRowOutput<'_>,
) -> Result<IrlsRowStatus, DeseqError> {
    let gene = input.gene;
    let p = input.x.ncols();
    let source_gene = input.rows.source(gene);
    let count_row = input.rows.input.counts.row(source_gene)?;
    if count_row.iter().all(|&count| count == 0) {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"),
        });
    }
    let IrlsScratch {
        counts: y,
        log_offsets,
        fitted_means: mu,
        weights: w,
        response: z,
        covariance: transient_covariance,
    } = scratch;
    y.clear();
    y.extend(count_row.iter().copied().map(f64::from));
    let nf = input.rows.input.offsets.row(source_gene)?;
    let log_nf = match input.shared_log_offsets {
        Some(shared) => shared,
        None => {
            log_offsets.clear();
            log_offsets.extend(nf.iter().copied().map(f64::ln_p));
            log_offsets.as_slice()
        }
    };
    let weight_row = input.rows.weights(source_gene)?;
    let dispersion = input.rows.input.dispersions[source_gene];
    let mut beta =
        initial_beta_for_rule(input.x, y, nf, Some(input.initial_lu), input.initial_beta)?;
    fitted_mu_into(
        input.x,
        &beta,
        nf,
        Some(input.options.min_mu),
        Some(log_nf),
        mu,
    )?;
    let likelihood = FixedDispersionLikelihood::new(count_row, dispersion);
    let mut dev_old = 0.0;
    let mut iter = 0_usize;
    let mut converged = false;

    for t in 0..input.options.maxit {
        iter += 1;
        working_weights_into(mu, dispersion, weight_row, w)?;
        working_response_into(mu, nf, y, Some(log_nf), z)?;
        let Some(next_beta) =
            solve_weighted_least_squares(input.x, w, z, input.ridge_lambda, input.options.solver)
        else {
            iter = input.options.maxit;
            break;
        };
        beta = next_beta;
        if beta
            .iter()
            .any(|value| !value.is_finite() || value.abs() > input.options.max_beta_abs)
        {
            iter = input.options.maxit;
            break;
        }
        fitted_mu_into(
            input.x,
            &beta,
            nf,
            Some(input.options.min_mu),
            Some(log_nf),
            mu,
        )?;
        let dev = -2.0 * likelihood.log_likelihood(mu, weight_row)?;
        let Some(conv_test) = irls_deviance_convergence_stat(dev, dev_old) else {
            iter = input.options.maxit;
            break;
        };
        if t > 0 && conv_test < input.options.beta_tol {
            converged = true;
            break;
        }
        dev_old = dev;
    }

    working_weights_into(mu, dispersion, weight_row, w)?;
    let covariance = covariance_and_hat_diagonal(input.x, w, input.ridge_lambda, output.hat);

    for (col, (out, value)) in output.beta.iter_mut().zip(beta.iter().copied()).enumerate() {
        *out = checked_log2_scale(value, col, "IRLS beta log2 scale")?;
    }
    let scaled_covariance = output.covariance.unwrap_or_else(|| {
        transient_covariance.resize(p * p, f64::NAN);
        transient_covariance
    });
    if let Some(covariance) = covariance {
        for (diagonal, out) in output.beta_se.iter_mut().enumerate() {
            let value = covariance[diagonal * p + diagonal];
            *out = checked_log2_standard_error(value, diagonal, "IRLS beta standard error")?;
        }
        for (idx, (out, value)) in scaled_covariance.iter_mut().zip(covariance).enumerate() {
            *out = checked_log2_covariance(value, idx, "IRLS beta covariance log2 scale")?;
        }
    } else {
        iter = input.options.maxit;
        scaled_covariance.fill(f64::NAN);
        // Computing the hat diagonal can fail partway through a row.
        output.hat.fill(f64::NAN);
    }
    fitted_mu_into(input.x, &beta, nf, None, Some(log_nf), mu)?;
    output.mu.copy_from_slice(mu);
    let converged = converged && iter < input.options.maxit;
    Ok(IrlsRowStatus {
        iter,
        converged,
        refit: requires_optim_fallback(
            converged,
            output.beta,
            scaled_covariance,
            input.options.use_optim,
            input.options.force_optim,
        ),
    })
}
