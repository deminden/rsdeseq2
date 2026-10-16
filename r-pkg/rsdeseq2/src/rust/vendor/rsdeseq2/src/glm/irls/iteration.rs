// Per-gene IRLS iterations and scratch reuse.

fn fit_fixed_dispersion_irls_row(input: IrlsRowInput<'_>) -> Result<IrlsRowOutput, DeseqError> {
    let gene = input.gene;
    let p = input.x.ncols();
    if input.counts.is_all_zero_gene(gene)? {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"),
        });
    }
    let y = input
        .counts
        .row(gene)?
        .iter()
        .copied()
        .map(f64::from)
        .collect::<Vec<_>>();
    let nf = input.normalization_factors.row(gene)?;
    let log_nf = nf.iter().copied().map(f64::ln).collect::<Vec<_>>();
    let weight_row = input.weights.map(|matrix| matrix.row(gene)).transpose()?;
    let dispersion = input.dispersions[gene];
    let mut beta = initial_beta_impl(input.x, &y, nf, Some(input.initial_lu))?;
    let mut mu = fitted_mu_impl_with_logs(
        input.x,
        &beta,
        nf,
        Some(input.options.min_mu),
        Some(&log_nf),
    )?;
    let likelihood = FixedDispersionLikelihood::new(input.counts.row(gene)?, dispersion);
    let mut dev_old = 0.0;
    let mut dev = 0.0;
    let mut iter = 0_usize;
    let mut converged = false;

    let mut w = Vec::with_capacity(mu.len());
    let mut z = Vec::with_capacity(mu.len());
    for t in 0..input.options.maxit {
        iter += 1;
        working_weights_into(&mu, dispersion, weight_row, &mut w)?;
        working_response_into(&mu, nf, &y, Some(&log_nf), &mut z)?;
        let Some(next_beta) =
            solve_weighted_least_squares(input.x, &w, &z, input.ridge_lambda, input.options.solver)
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
            Some(&log_nf),
            &mut mu,
        )?;
        dev = -2.0 * likelihood.log_likelihood(&mu, weight_row)?;
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

    let mut row = IrlsRowOutput {
        beta: Vec::with_capacity(p),
        beta_se: Vec::with_capacity(p),
        beta_covariance: Vec::with_capacity(p * p),
        mu: Vec::with_capacity(input.counts.n_samples()),
        hat: Vec::with_capacity(input.counts.n_samples()),
        iter,
        converged: false,
    };

    let w = working_weights(&mu, dispersion, weight_row)?;
    let Some((beta_covariance, hat_diag)) =
        covariance_and_hat_diagonal(input.x, &w, input.ridge_lambda)
    else {
        row.iter = input.options.maxit;
        row.beta_se.extend(std::iter::repeat_n(f64::NAN, p));
        row.beta_covariance
            .extend(std::iter::repeat_n(f64::NAN, p * p));
        row.hat
            .extend(std::iter::repeat_n(f64::NAN, input.counts.n_samples()));
        for (col, value) in beta.iter().copied().enumerate() {
            row.beta
                .push(checked_log2_scale(value, col, "IRLS beta log2 scale")?);
        }
        row.mu.extend(
            fitted_mu_impl_with_logs(input.x, &beta, nf, None, Some(&log_nf))?
                .iter()
                .copied(),
        );
        return Ok(row);
    };

    for (col, value) in beta.iter().copied().enumerate() {
        row.beta
            .push(checked_log2_scale(value, col, "IRLS beta log2 scale")?);
    }
    for diagonal in 0..p {
        let value = beta_covariance[diagonal * p + diagonal];
        row.beta_se.push(checked_log2_standard_error(
            value,
            diagonal,
            "IRLS beta standard error",
        )?);
    }
    for (idx, value) in beta_covariance.into_iter().enumerate() {
        row.beta_covariance.push(checked_log2_covariance(
            value,
            idx,
            "IRLS beta covariance log2 scale",
        )?);
    }
    row.mu.extend(
        fitted_mu_impl_with_logs(input.x, &beta, nf, None, Some(&log_nf))?
            .iter()
            .copied(),
    );
    row.hat.extend(hat_diag);
    row.converged = converged && row.iter < input.options.maxit;
    let _ = dev;
    Ok(row)
}
