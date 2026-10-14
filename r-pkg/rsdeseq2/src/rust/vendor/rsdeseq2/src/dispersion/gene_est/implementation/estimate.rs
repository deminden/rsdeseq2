/// Estimate gene-wise dispersions using DESeq2's linear-mu branch shape.
///
/// This implements a clean Rust subset of `estimateDispersionsGeneEst`:
/// base normalized counts are projected through the supplied design matrix,
/// raw means are reconstructed from size factors or gene/sample normalization
/// factors, rough/moments starts are bounded, and each gene's dispersion is
/// optimized on a two-pass log-alpha grid with optional Cox-Reid correction and
/// without priors. General iterative GLM mean refitting remains future work.
pub fn estimate_gene_wise_dispersions_linear_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let max_disp = max_dispersion(options, input.counts.n_samples());
    let normalized_mu = linear_model_mu(input.normalized_counts, input.design)?;
    let rough_disp = rough_dispersion_estimates(input.normalized_counts, input.design)?;
    let moments_disp = match input.normalization_factors {
        Some(normalization_factors) => moments_dispersion_estimates_with_normalization_factors(
            input.base_mean,
            input.base_var,
            normalization_factors,
            Some(input.all_zero),
        )?,
        None => moments_dispersion_estimates(input.base_mean, input.base_var, input.size_factors)?,
    };
    let initial_disp =
        initial_dispersion_estimates(&rough_disp, &moments_disp, options.min_disp, max_disp)?;

    let mut mu_values = vec![f64::NAN; input.counts.n_genes() * input.counts.n_samples()];
    let mut disp_gene_est = vec![f64::NAN; input.counts.n_genes()];
    let mut disp_iter = vec![0; input.counts.n_genes()];
    let mut converged = vec![false; input.counts.n_genes()];

    for gene in 0..input.counts.n_genes() {
        if input.all_zero[gene] {
            continue;
        }
        let mu_start = gene * input.counts.n_samples();
        let normalization_factor_row = input
            .normalization_factors
            .map(|normalization_factors| normalization_factors.row(gene))
            .transpose()?;
        for sample in 0..input.counts.n_samples() {
            let factor = normalization_factor_row
                .map(|row| row[sample])
                .unwrap_or(input.size_factors[sample]);
            let value = normalized_mu.row(gene)?[sample] * factor;
            mu_values[mu_start + sample] = value.max(options.min_mu);
        }
        let row_mu = &mu_values[mu_start..mu_start + input.counts.n_samples()];
        let (estimate, iterations, is_converged) = fit_dispersion_for_gene(
            input.counts.row(gene)?,
            row_mu,
            input.design,
            initial_disp[gene],
            options,
            input.counts.n_samples(),
        )?;
        disp_gene_est[gene] = estimate.clamp(options.min_disp, max_disp);
        disp_iter[gene] = iterations;
        converged[gene] = is_converged;
    }

    Ok(GeneWiseDispersionOutput {
        disp_gene_est,
        disp_iter,
        rough_disp,
        moments_disp,
        initial_disp,
        mu: RowMajorMatrix::from_row_major(
            input.counts.n_genes(),
            input.counts.n_samples(),
            mu_values,
        )?,
        converged,
    })
}

/// Estimate gene-wise dispersions using one or more GLM mean-refit iterations.
///
/// This follows the non-`linearMu` branch shape of DESeq2's
/// `estimateDispersionsGeneEst`: rough/moments estimates initialize
/// `alpha_hat`, non-all-zero rows alternate between fixed-dispersion NB GLM
/// mean fitting and fixed-mean dispersion optimization, and rows stop
/// refitting when the log-dispersion move is at most `0.05`. When
/// row-normalized observation weights are supplied, they are used in the
/// fixed-dispersion IRLS mean fit and the fixed-mean likelihood objective;
/// Cox-Reid terms use DESeq2's thresholded weighted sample subset.
/// glmGamPoi fitting remains a future high-level branch.
pub fn estimate_gene_wise_dispersions_glm_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
    irls_options: IrlsOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let max_disp = max_dispersion(options, input.counts.n_samples());
    let rough_disp = rough_dispersion_estimates(input.normalized_counts, input.design)?;
    let moments_disp = match input.normalization_factors {
        Some(normalization_factors) => moments_dispersion_estimates_with_normalization_factors(
            input.base_mean,
            input.base_var,
            normalization_factors,
            Some(input.all_zero),
        )?,
        None => moments_dispersion_estimates(input.base_mean, input.base_var, input.size_factors)?,
    };
    let initial_disp =
        initial_dispersion_estimates(&rough_disp, &moments_disp, options.min_disp, max_disp)?;

    let mut alpha_hat = initial_disp.clone();
    let mut alpha_hat_new = initial_disp.clone();
    let alpha_init = initial_disp.clone();
    let fitting_gene_order = input
        .all_zero
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(gene, all_zero)| (!all_zero).then_some(gene))
        .collect::<Vec<_>>();
    let mut fitidx = input
        .all_zero
        .iter()
        .map(|all_zero| !all_zero)
        .collect::<Vec<_>>();
    let mut mu_values = vec![f64::NAN; input.counts.n_genes() * input.counts.n_samples()];
    let mut disp_iter = vec![0; input.counts.n_genes()];
    let mut initial_lp = vec![f64::NAN; input.counts.n_genes()];
    let mut last_lp = vec![f64::NAN; input.counts.n_genes()];

    let mut mean_options = irls_options;
    mean_options.min_mu = options.min_mu;
    mean_options.r_optim_compat = false;

    for _ in 0..options.niter {
        let fit_genes = fitidx
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(gene, should_fit)| should_fit.then_some(gene))
            .collect::<Vec<_>>();
        if fit_genes.is_empty() {
            break;
        }

        let compact_counts = compact_counts_rows(input.counts, &fit_genes)?;
        let compact_disp = compact_gene_values(&alpha_hat, &fit_genes)?;
        let compact_weights = input
            .observation_weights
            .map(|weights| compact_matrix_rows(weights, &fit_genes))
            .transpose()?;
        let fit = match input.normalization_factors {
            Some(normalization_factors) => {
                let compact_factors = compact_matrix_rows(normalization_factors, &fit_genes)?;
                fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
                    &compact_counts,
                    input.design,
                    &compact_factors,
                    &compact_disp,
                    compact_weights.as_ref(),
                    mean_options.clone(),
                )?
            }
            None => fit_fixed_dispersion_irls_with_weights(
                &compact_counts,
                input.design,
                input.size_factors,
                &compact_disp,
                compact_weights.as_ref(),
                mean_options.clone(),
            )?,
        };

        let updates = fit_genes
            .par_iter()
            .copied()
            .enumerate()
            .map(
                |(compact_row, gene)| -> Result<
                    (usize, Vec<f64>, GeneDispersionFitDiagnostics),
                    DeseqError,
                > {
                    let fit_mu_raw = fit.mu.row(compact_row)?;
                    let fit_mu = fit_mu_raw
                        .iter()
                        .copied()
                        .map(|value| value.max(options.min_mu))
                        .collect::<Vec<_>>();
                    // DESeq2 passes the full non-all-zero weight matrix into fitDisp
                    // even when counts/mu are subset by fitidx; the C++ then indexes
                    // weights by compact row position.
                    let weight_row = input
                        .observation_weights
                        .map(|weights| weights.row(fitting_gene_order[compact_row]))
                        .transpose()?;
                    let diagnostics = fit_dispersion_for_gene_detailed_with_weights(
                        input.counts.row(gene)?,
                        &fit_mu,
                        input.design,
                        alpha_hat[gene],
                        options,
                        input.counts.n_samples(),
                        weight_row,
                    )?;
                    Ok((gene, fit_mu, diagnostics))
                },
            )
            .collect::<Result<Vec<_>, DeseqError>>()?;

        for (gene, fit_mu, diagnostics) in updates {
            let start = gene * input.counts.n_samples();
            mu_values[start..start + input.counts.n_samples()].copy_from_slice(&fit_mu);
            alpha_hat_new[gene] = diagnostics.estimate.min(max_disp);
            disp_iter[gene] = diagnostics.iterations;
            initial_lp[gene] = diagnostics.initial_lp;
            last_lp[gene] = diagnostics.last_lp;
        }

        fitidx = input
            .all_zero
            .iter()
            .copied()
            .enumerate()
            .map(|(gene, all_zero)| {
                if all_zero {
                    return false;
                }
                let move_size = (alpha_hat_new[gene] / alpha_hat[gene]).ln().abs();
                move_size.is_finite() && move_size > 0.05
            })
            .collect();
        alpha_hat.clone_from(&alpha_hat_new);
        if !fitidx.iter().any(|should_fit| *should_fit) {
            break;
        }
    }

    let mut disp_gene_est = alpha_hat;
    if options.niter == 1 {
        for gene in 0..input.counts.n_genes() {
            if input.all_zero[gene] || !initial_lp[gene].is_finite() || !last_lp[gene].is_finite() {
                continue;
            }
            if last_lp[gene] < initial_lp[gene] + initial_lp[gene].abs() / 1.0e6 {
                disp_gene_est[gene] = alpha_init[gene];
            }
        }
    }

    let final_rows = (0..input.counts.n_genes())
        .into_par_iter()
        .map(|gene| -> Result<(usize, f64, bool), DeseqError> {
            if input.all_zero[gene] {
                return Ok((gene, f64::NAN, false));
            }
            let converged = disp_iter[gene] < options.maxit && disp_iter[gene] != 1;
            let mut estimate = disp_gene_est[gene];
            if !converged && estimate > options.min_disp * 10.0 {
                let mu = &mu_values[gene * input.counts.n_samples()
                    ..gene * input.counts.n_samples() + input.counts.n_samples()];
                let weight_row = input
                    .observation_weights
                    .map(|weights| weights.row(gene))
                    .transpose()?;
                estimate = fit_dispersion_grid_inner(DispersionOptimizerInput {
                    counts: input.counts.row(gene)?,
                    mu,
                    design: Some(input.design),
                    initial_dispersion: estimate,
                    options,
                    n_samples: input.counts.n_samples(),
                    prior: None,
                    weights: weight_row,
                })?
                .0;
            }
            Ok((gene, estimate.clamp(options.min_disp, max_disp), converged))
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;

    let mut converged = vec![false; input.counts.n_genes()];
    for (gene, estimate, row_converged) in final_rows {
        disp_gene_est[gene] = estimate;
        converged[gene] = row_converged;
    }

    Ok(GeneWiseDispersionOutput {
        disp_gene_est,
        disp_iter,
        rough_disp,
        moments_disp,
        initial_disp,
        mu: RowMajorMatrix::from_row_major(
            input.counts.n_genes(),
            input.counts.n_samples(),
            mu_values,
        )?,
        converged,
    })
}

/// Project normalized counts onto the supplied design matrix.
///
/// This is the Rust analogue of DESeq2's `linearModelMu` helper for row-wise
/// fitted values, using `Y X (X'X)^-1 X'`.
pub fn linear_model_mu(
    normalized_counts: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if normalized_counts.n_cols() != design.n_samples() {
        return Err(invalid_dimensions(
            "linear mu samples",
            design.n_samples(),
            normalized_counts.n_cols(),
        ));
    }
    let y = DMatrix::from_row_slice(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        normalized_counts.as_slice(),
    );
    let x = DMatrix::from_row_slice(
        design.n_samples(),
        design.n_coefficients(),
        design.matrix().as_slice(),
    );
    let xtx = x.transpose() * &x;
    let Some(xtx_inverse) = xtx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "linear mu design rank".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    let hat = &x * xtx_inverse * x.transpose();
    let fitted = y * hat;
    let mut values = Vec::with_capacity(normalized_counts.n_rows() * normalized_counts.n_cols());
    for row in 0..normalized_counts.n_rows() {
        for col in 0..normalized_counts.n_cols() {
            values.push(fitted[(row, col)]);
        }
    }
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// DESeq2-style rough dispersion estimates from normalized counts.
pub fn rough_dispersion_estimates(
    normalized_counts: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
) -> Result<Vec<f64>, DeseqError> {
    if design.n_samples() <= design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion residual degrees of freedom".to_string(),
            expected: design.n_coefficients() + 1,
            actual: design.n_samples(),
        });
    }
    let mu = linear_model_mu(normalized_counts, design)?;
    let residual_df = (design.n_samples() - design.n_coefficients()) as f64;
    let mut estimates = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let y = normalized_counts.row(gene)?;
        let mu = mu.row(gene)?;
        let mut sum = 0.0;
        for (sample, (count, fitted)) in y.iter().copied().zip(mu.iter().copied()).enumerate() {
            let fitted = fitted.max(1.0);
            let residual = checked_sub(count, fitted, sample, "rough dispersion residual")?;
            let inv_fitted = fitted.recip();
            let relative_residual = checked_mul(
                residual,
                inv_fitted,
                sample,
                "rough dispersion relative residual",
            )?;
            let relative_square = checked_mul(
                relative_residual,
                relative_residual,
                sample,
                "rough dispersion relative residual square",
            )?;
            let term = checked_sub(relative_square, inv_fitted, sample, "rough dispersion term")?;
            checked_matrix_add_assign(&mut sum, term, sample, "rough dispersion row sum")?;
        }
        let average = checked_div(sum, residual_df, gene, "rough dispersion row mean")?;
        estimates.push(average.max(0.0));
    }
    Ok(estimates)
}

/// DESeq2-style moments dispersion estimates.
pub fn moments_dispersion_estimates(
    base_mean: &[f64],
    base_var: &[f64],
    size_factors: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    validate_size_factors(size_factors)?;
    let inverse_sum = checked_sum_indexed(
        size_factors.iter().copied().map(f64::recip),
        "moments dispersion inverse size-factor sum",
    )?;
    let xim = checked_div(
        inverse_sum,
        size_factors.len() as f64,
        0,
        "moments dispersion inverse size-factor mean",
    )?;
    moments_dispersion_estimates_with_xim(base_mean, base_var, xim)
}

/// DESeq2-style moments dispersion estimates with gene/sample normalization factors.
///
/// This follows `momentsDispEstimate`: when normalization factors are present,
/// `xim = mean(1 / colMeans(normalizationFactors))`. If `all_zero` is supplied,
/// all-zero rows are excluded from the column means, matching the fact that
/// DESeq2 calls the helper on `objectNZ`.
pub fn moments_dispersion_estimates_with_normalization_factors(
    base_mean: &[f64],
    base_var: &[f64],
    normalization_factors: &RowMajorMatrix<f64>,
    all_zero: Option<&[bool]>,
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    if normalization_factors.n_rows() != base_mean.len() {
        return Err(invalid_dimensions(
            "moments dispersion normalization-factor rows",
            base_mean.len(),
            normalization_factors.n_rows(),
        ));
    }
    if let Some(all_zero) = all_zero
        && all_zero.len() != base_mean.len() {
            return Err(invalid_dimensions(
                "moments dispersion allZero",
                base_mean.len(),
                all_zero.len(),
            ));
        }
    let xim = normalization_factor_moments_xim(normalization_factors, all_zero)?;
    moments_dispersion_estimates_with_xim(base_mean, base_var, xim)
}

fn moments_dispersion_estimates_with_xim(
    base_mean: &[f64],
    base_var: &[f64],
    xim: f64,
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    if !xim.is_finite() || xim <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "moments dispersion normalization factor summary must be finite and positive"
                .to_string(),
        });
    }
    let mut estimates = Vec::with_capacity(base_mean.len());
    for (gene, (mean, variance)) in base_mean
        .iter()
        .copied()
        .zip(base_var.iter().copied())
        .enumerate()
    {
        if mean > 0.0 {
            let inv_mean = mean.recip();
            let xim_mean = checked_mul(xim, mean, gene, "moments dispersion xim mean")?;
            let centered = variance - xim_mean;
            if !centered.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "moments dispersion centered variance".to_string(),
                    index: Some(gene),
                    value: centered,
                });
            }
            let inv_square = checked_mul(
                inv_mean,
                inv_mean,
                gene,
                "moments dispersion inverse mean square",
            )?;
            estimates.push(checked_mul(
                centered,
                inv_square,
                gene,
                "moments dispersion estimate",
            )?);
        } else {
            estimates.push(f64::NAN);
        }
    }
    Ok(estimates)
}

fn normalization_factor_moments_xim(
    normalization_factors: &RowMajorMatrix<f64>,
    all_zero: Option<&[bool]>,
) -> Result<f64, DeseqError> {
    let mut col_sums = vec![0.0; normalization_factors.n_cols()];
    let mut n_rows_used = 0_usize;
    for row in 0..normalization_factors.n_rows() {
        if all_zero.is_some_and(|flags| flags[row]) {
            continue;
        }
        for (sample, value) in normalization_factors.row(row)?.iter().copied().enumerate() {
            validate_normalization_factor(value, sample)?;
            checked_matrix_add_assign(
                &mut col_sums[sample],
                value,
                sample,
                "moments dispersion normalization-factor column sum",
            )?;
        }
        n_rows_used += 1;
    }
    if n_rows_used == 0 {
        return Err(DeseqError::InvalidCounts {
            reason: "no non-all-zero rows available for normalization-factor moments estimate"
                .to_string(),
        });
    }
    let mut inverse_col_mean_sum = 0.0;
    for (sample, sum) in col_sums.iter().copied().enumerate() {
        let col_mean = checked_div(
            sum,
            n_rows_used as f64,
            sample,
            "moments dispersion normalization-factor column mean",
        )?;
        if !col_mean.is_finite() || col_mean <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!(
                    "normalization-factor column mean at sample {sample} must be finite and positive"
                ),
            });
        }
        checked_matrix_add_assign(
            &mut inverse_col_mean_sum,
            col_mean.recip(),
            sample,
            "moments dispersion inverse normalization-factor mean sum",
        )?;
    }
    checked_div(
        inverse_col_mean_sum,
        normalization_factors.n_cols() as f64,
        0,
        "moments dispersion inverse normalization-factor mean",
    )
}
