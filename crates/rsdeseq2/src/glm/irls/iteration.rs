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
    coefficients: DVector<f64>,
    least_squares: LeastSquaresScratch,
    covariance: CovarianceScratch,
}

pub(super) fn fit_fixed_dispersion_irls_row<const SELECTED: bool>(
    scratch: &mut IrlsScratch,
    input: IrlsRowInput<'_, SELECTED>,
    mut output: IrlsRowOutput<'_>,
) -> Result<IrlsRowStatus, DeseqError> {
    let gene = input.gene;
    let source_gene = input.rows.source(gene);
    let count_row = input.rows.input.counts.row(source_gene)?;
    if count_row.iter().all(|&count| count == 0) {
        return Err(DeseqError::AllZeroGene { gene, source_gene });
    }
    let IrlsScratch {
        counts: y,
        log_offsets,
        fitted_means: mu,
        weights: w,
        response: z,
        coefficients: beta,
        least_squares,
        covariance: covariance_scratch,
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
    if beta.len() != input.x.ncols() {
        *beta = DVector::zeros(input.x.ncols());
    }
    input.initial_beta.initialize(y, nf, z, beta)?;
    fitted_mu_into(
        input.x,
        beta,
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
        if solve_weighted_least_squares(
            least_squares,
            input.x,
            w,
            z,
            input.ridge_lambda,
            input.options.solver,
            beta,
        )
        .is_none()
        {
            iter = input.options.maxit;
            break;
        }
        if beta
            .iter()
            .any(|value| !value.is_finite() || value.abs() > input.options.max_beta_abs)
        {
            iter = input.options.maxit;
            break;
        }
        fitted_mu_into(
            input.x,
            beta,
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
    let covariance = covariance_and_hat_diagonal(
        covariance_scratch,
        input.x,
        w,
        input.ridge_lambda,
        output.hat.as_deref_mut(),
        input.max_design,
        output.covariance.is_some(),
    );

    for (col, (out, value)) in output.beta.iter_mut().zip(beta.iter().copied()).enumerate() {
        *out = checked_log2_scale(value, col, "IRLS beta log2 scale")?;
    }
    let variance_positive = match covariance {
        Some(FitCovariance::PositiveDiagonal) => true,
        Some(FitCovariance::Full(covariance)) => {
            if let Some(beta_se) = output.beta_se {
                for (diagonal, out) in beta_se.iter_mut().enumerate() {
                    let value = covariance[(diagonal, diagonal)];
                    *out =
                        checked_log2_standard_error(value, diagonal, "IRLS beta standard error")?;
                }
            }
            retain_covariance_on_log2_scale(
                covariance,
                output.covariance.as_deref_mut(),
                0,
                "IRLS beta covariance log2 scale",
            )?
        }
        None => {
            iter = input.options.maxit;
            if let Some(covariance) = output.covariance {
                covariance.fill(f64::NAN);
            }
            // Computing the hat diagonal can fail partway through a row.
            if let Some(hat) = output.hat {
                hat.fill(f64::NAN);
            }
            false
        }
    };
    fitted_mu_into(input.x, beta, nf, None, Some(log_nf), mu)?;
    output.mu.copy_from_slice(mu);
    let converged = converged && iter < input.options.maxit;
    Ok(IrlsRowStatus {
        iter,
        converged,
        refit: requires_optim_fallback(
            converged,
            output.beta,
            variance_positive,
            input.options.use_optim,
            input.options.force_optim,
        ),
    })
}
