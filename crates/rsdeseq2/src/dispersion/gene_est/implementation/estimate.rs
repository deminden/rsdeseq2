use super::*;
use crate::errors::{ErrorKind, Requirement, invalid_parameter};
use crate::matrix::{try_alloc_matrix, try_dmatrix_from_row_major, try_zero_dmatrix};
use crate::numeric::math::ProfileMath;
use nalgebra::DMatrixView;
/// Estimate gene-wise dispersions using DESeq2's linear-mu branch shape.
///
/// This implements a clean Rust subset of `estimateDispersionsGeneEst`:
/// base normalized counts are projected through the supplied design matrix,
/// raw means are reconstructed from size factors or gene/sample normalization
/// factors, rough/moments starts are bounded, and each gene's dispersion is
/// optimized with the configured search and optional Cox-Reid correction.
pub fn estimate_gene_wise_dispersions_linear_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let max_disp = max_dispersion(options, input.counts.n_samples());
    let normalized_mu = linear_model_mu(input.normalized_counts, input.design)?;
    let residual_df = rough_dispersion_residual_df(input.design)?;
    let rough_disp =
        rough_dispersion_from_means(input.normalized_counts, &normalized_mu, residual_df)?;
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

    let mut mu_values = normalized_mu.into_values();
    let mut disp_gene_est = vec![f64::NAN; input.counts.n_genes()];
    let mut disp_iter = vec![0; input.counts.n_genes()];
    let mut converged = vec![false; input.counts.n_genes()];

    for (gene, row_mu) in mu_values
        .chunks_exact_mut(input.counts.n_samples())
        .enumerate()
    {
        if input.all_zero[gene] {
            row_mu.fill(f64::NAN);
            continue;
        }
        let factors = match input.normalization_factors {
            Some(factors) => factors.row(gene)?,
            None => input.size_factors,
        };
        for (value, factor) in row_mu.iter_mut().zip(factors) {
            *value = (*value * factor).max(options.min_mu);
        }
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
pub fn estimate_gene_wise_dispersions_glm_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
    irls_options: IrlsOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    estimate_gene_wise_dispersions_with_mean_fit(input, options, irls_options, false)
}

/// Estimate dispersions with DESeq2's automatic fitted-mean selection.
///
/// Unweighted designs with as many distinct rows as coefficients use linear
/// fitted means. Other designs and all weighted fits use negative-binomial
/// GLM means. Both branches share the same dispersion optimization, likelihood
/// improvement check, and grid fallback.
pub fn estimate_gene_wise_dispersions_auto(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
    irls_options: IrlsOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    let linear_mu =
        input.observation_weights.is_none() && design_uses_linear_dispersion_means(input.design)?;
    estimate_gene_wise_dispersions_with_mean_fit(input, options, irls_options, linear_mu)
}

/// Whether the design is saturated across its distinct sample groups.
pub fn design_uses_linear_dispersion_means(design: &DesignMatrix) -> Result<bool, DeseqError> {
    design.validate_full_rank("automatic dispersion mean selection")?;
    let mut groups: Vec<&[f64]> = Vec::new();
    for sample in 0..design.n_samples() {
        let row = design.matrix().row(sample)?;
        if !groups.contains(&row) {
            groups.push(row);
            if groups.len() > design.n_coefficients() {
                return Ok(false);
            }
        }
    }
    Ok(groups.len() == design.n_coefficients())
}

pub(super) fn estimate_gene_wise_dispersions_with_mean_fit(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
    irls_options: IrlsOptions,
    linear_mu: bool,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let (rough_disp, linear_means) = if linear_mu {
        let means = linear_model_mu(input.normalized_counts, input.design)?;
        let residual_df = rough_dispersion_residual_df(input.design)?;
        let rough_disp = rough_dispersion_from_means(input.normalized_counts, &means, residual_df)?;
        let mut values = means.into_values();
        let samples = input.counts.n_samples();
        for (gene, row) in values.chunks_exact_mut(samples).enumerate() {
            let factors = match input.normalization_factors {
                Some(factors) => factors.row(gene)?,
                None => input.size_factors,
            };
            for (value, factor) in row.iter_mut().zip(factors) {
                *value *= factor;
            }
        }
        let means = RowMajorMatrix::from_row_major(input.counts.n_genes(), samples, values)?;
        (rough_disp, Some(means))
    } else {
        (
            rough_dispersion_estimates(input.normalized_counts, input.design)?,
            None,
        )
    };
    let max_disp = max_dispersion(options, input.counts.n_samples());
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
    let mut mu_values = {
        let mut values = try_alloc_matrix(
            input.counts.n_genes(),
            input.counts.n_samples(),
            "gene-wise dispersion means",
        )?;
        values.resize(input.counts.n_genes() * input.counts.n_samples(), f64::NAN);
        values
    };
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

        let irls_means;
        let fit_mu = if let Some(means) = &linear_means {
            means
        } else {
            let offsets = match input.normalization_factors {
                Some(factors) => Offsets::PerCell(factors),
                None => Offsets::PerSample(input.size_factors),
            };
            irls_means = fit_selected_irls_means(
                GlmInput {
                    counts: input.counts,
                    offsets,
                    dispersions: &alpha_hat,
                    weights: input.observation_weights,
                },
                input.design,
                mean_options.clone(),
                &fit_genes,
            )?;
            &irls_means
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
                    let mean_row = if linear_means.is_some() { gene } else { compact_row };
                    let fit_mu_raw = fit_mu.row(mean_row)?;
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
                let move_size = (alpha_hat_new[gene] / alpha_hat[gene]).ln_p().abs();
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
                estimate = grid_search_dispersion(DispersionFitInput {
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
    let samples = design.n_samples();
    let coefficients = design.n_coefficients();
    let x = try_dmatrix_from_row_major(design.matrix(), "linear mean design")?;
    // The row-major design slice has the column-major layout of its transpose.
    let xt = DMatrixView::from_slice(design.matrix().as_slice(), coefficients, samples);
    let mut xtx = try_zero_dmatrix(coefficients, coefficients, "linear mean crossproduct")?;
    xt.mul_to(&x, &mut xtx);
    if !xtx.try_inverse_mut() {
        return Err(DeseqError::InvalidDimensions {
            context: "linear mu design rank".to_string(),
            expected: coefficients,
            actual: 0,
        });
    }
    let mut left = try_zero_dmatrix(samples, coefficients, "linear mean projection left factor")?;
    x.mul_to(&xtx, &mut left);
    let mut hat = try_zero_dmatrix(samples, samples, "linear mean projection")?;
    left.mul_to(&xt, &mut hat);
    let y = try_dmatrix_from_row_major(normalized_counts, "linear mean normalized counts")?;
    let mut fitted = try_zero_dmatrix(normalized_counts.n_rows(), samples, "linear fitted means")?;
    y.mul_to(&hat, &mut fitted);
    drop(y);
    drop(hat);
    let mut values = try_alloc_matrix(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        "linear fitted means",
    )?;
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
    let residual_df = rough_dispersion_residual_df(design)?;
    let mu = linear_model_mu(normalized_counts, design)?;
    rough_dispersion_from_means(normalized_counts, &mu, residual_df)
}

fn rough_dispersion_residual_df(design: &DesignMatrix) -> Result<f64, DeseqError> {
    if design.n_samples() <= design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion residual degrees of freedom".to_string(),
            expected: design.n_coefficients() + 1,
            actual: design.n_samples(),
        });
    }
    Ok((design.n_samples() - design.n_coefficients()) as f64)
}

fn rough_dispersion_from_means(
    normalized_counts: &RowMajorMatrix<f64>,
    mu: &RowMajorMatrix<f64>,
    residual_df: f64,
) -> Result<Vec<f64>, DeseqError> {
    let mut estimates = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let y = normalized_counts.row(gene)?;
        let mu = mu.row(gene)?;
        let mut sum = 0.0;
        for (sample, (count, fitted)) in y.iter().copied().zip(mu.iter().copied()).enumerate() {
            let fitted = fitted.max(1.0);
            let residual = checked::try_sub(count, fitted)
                .nonfinite_at("rough dispersion residual", Some(sample))?;
            let inv_fitted = fitted.recip();
            let relative_residual = checked::try_mul(residual, inv_fitted)
                .nonfinite_at("rough dispersion relative residual", Some(sample))?;
            let relative_square = checked::try_mul(relative_residual, relative_residual)
                .nonfinite_at("rough dispersion relative residual square", Some(sample))?;
            let term = checked::try_sub(relative_square, inv_fitted)
                .nonfinite_at("rough dispersion term", Some(sample))?;
            checked_matrix_add_assign(&mut sum, term, sample, "rough dispersion row sum")?;
        }
        let average = checked::try_div(sum, residual_df)
            .nonfinite_at("rough dispersion row mean", Some(gene))?;
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
    validate_size_factor_values(size_factors)?;
    let inverse_sum = checked_sum_indexed(
        size_factors.iter().copied().map(f64::recip),
        "moments dispersion inverse size-factor sum",
    )?;
    let xim = checked::try_div(inverse_sum, size_factors.len() as f64)
        .nonfinite_at("moments dispersion inverse size-factor mean", Some(0))?;
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
        && all_zero.len() != base_mean.len()
    {
        return Err(invalid_dimensions(
            "moments dispersion allZero",
            base_mean.len(),
            all_zero.len(),
        ));
    }
    let xim = normalization_factor_moments_xim(normalization_factors, all_zero)?;
    moments_dispersion_estimates_with_xim(base_mean, base_var, xim)
}

pub(super) fn moments_dispersion_estimates_with_xim(
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
        return Err(invalid_parameter(
            ErrorKind::Normalization,
            "moments dispersion normalization factor summary",
            Requirement::FiniteAndPositive,
        ));
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
            let xim_mean = checked::try_mul(xim, mean)
                .nonfinite_at("moments dispersion xim mean", Some(gene))?;
            let centered = variance - xim_mean;
            if !centered.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "moments dispersion centered variance".to_string(),
                    index: Some(gene),
                    value: centered,
                });
            }
            let inv_square = checked::try_mul(inv_mean, inv_mean)
                .nonfinite_at("moments dispersion inverse mean square", Some(gene))?;
            estimates.push(
                checked::try_mul(centered, inv_square)
                    .nonfinite_at("moments dispersion estimate", Some(gene))?,
            );
        } else {
            estimates.push(f64::NAN);
        }
    }
    Ok(estimates)
}

pub(super) fn normalization_factor_moments_xim(
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
            validate_normalization_factor_value(value, sample)?;
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
        let col_mean = checked::try_div(sum, n_rows_used as f64).nonfinite_at(
            "moments dispersion normalization-factor column mean",
            Some(sample),
        )?;
        if !col_mean.is_finite() || col_mean <= 0.0 {
            return Err(invalid_parameter(
                ErrorKind::Normalization,
                format!("normalization-factor column mean at sample {sample}"),
                Requirement::FiniteAndPositive,
            ));
        }
        checked_matrix_add_assign(
            &mut inverse_col_mean_sum,
            col_mean.recip(),
            sample,
            "moments dispersion inverse normalization-factor mean sum",
        )?;
    }
    checked::try_div(inverse_col_mean_sum, normalization_factors.n_cols() as f64).nonfinite_at(
        "moments dispersion inverse normalization-factor mean",
        Some(0),
    )
}
