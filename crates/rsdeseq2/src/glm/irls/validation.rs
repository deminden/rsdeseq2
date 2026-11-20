use super::*;
use crate::errors::validate_matrix_shape;
// Input validation, normalization factors and output scaling.

/// Validate the sample factors for a shared IRLS offset.
pub(in crate::glm) fn validate_size_factor_offsets(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<(), DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    if counts.n_genes() > 0 {
        for (idx, factor) in size_factors.iter().copied().enumerate() {
            validate_positive_finite(factor, "size factor", idx)?;
        }
    }
    Ok(())
}

pub(super) fn validate_nf_irls_inputs<const SELECTED: bool>(
    fit_rows: IrlsRows<'_, SELECTED>,
    design: &DesignMatrix,
    options: &IrlsOptions,
) -> Result<(), DeseqError> {
    let GlmInput {
        counts,
        offsets,
        dispersions,
        weights,
    } = fit_rows.input;
    if design.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "design rows",
            counts.n_samples(),
            design.n_samples(),
        ));
    }
    let shape = match offsets {
        // Size factors stand for the genes x samples matrix they used to fill.
        Offsets::PerSample(factors) => (counts.n_genes(), factors.len()),
        Offsets::PerCell(matrix) => (matrix.n_rows(), matrix.n_cols()),
    };
    validate_matrix_shape(
        "normalization factors",
        (counts.n_genes(), counts.n_samples()),
        shape,
    )?;
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    if SELECTED {
        for (gene, &source_gene) in fit_rows.selection.iter().enumerate() {
            let dispersion = *dispersions.get(source_gene).ok_or_else(|| {
                invalid_dimensions(
                    "IRLS gene index",
                    dispersions.len(),
                    source_gene.saturating_add(1),
                )
            })?;
            validate_positive_finite(dispersion, "dispersion", gene)?;
        }
    } else {
        for (gene, &dispersion) in dispersions.iter().enumerate() {
            validate_positive_finite(dispersion, "dispersion", gene)?;
        }
    }
    if let Some(weights) = weights {
        validate_matrix_shape(
            "weights",
            (counts.n_genes(), counts.n_samples()),
            (weights.n_rows(), weights.n_cols()),
        )?;
        if SELECTED {
            for (gene, &source_gene) in fit_rows.selection.iter().enumerate() {
                let row = weights.row(source_gene)?;
                for (sample, &weight) in row.iter().enumerate() {
                    validate_nonnegative_finite(
                        weight,
                        "weight",
                        gene * weights.n_cols() + sample,
                    )?;
                }
            }
        } else {
            for (index, &weight) in weights.as_slice().iter().enumerate() {
                validate_nonnegative_finite(weight, "weight", index)?;
            }
        }
    }
    if !options.beta_tol.is_finite()
        || options.beta_tol <= 0.0
        || options.maxit == 0
        || options.optim_maxit == 0
        || !options.optim_tol.is_finite()
        || options.optim_tol < 0.0
        || !options.min_mu.is_finite()
        || options.min_mu <= 0.0
        || !options.max_beta_abs.is_finite()
        || options.max_beta_abs <= 0.0
        || !options.ridge_lambda.is_finite()
        || options.ridge_lambda < 0.0
    {
        return Err(DeseqError::UnsupportedFeature {
            feature: "invalid IRLS options".to_string(),
        });
    }
    match offsets {
        // Already checked finite and positive by validate_size_factor_offsets.
        Offsets::PerSample(_) => {}
        Offsets::PerCell(matrix) if SELECTED => {
            for (gene, &source_gene) in fit_rows.selection.iter().enumerate() {
                let row = matrix.row(source_gene)?;
                for (sample, &value) in row.iter().enumerate() {
                    if !value.is_finite() {
                        return Err(DeseqError::NonFiniteValue {
                            context: "normalization factors".into(),
                            index: Some(gene * matrix.n_cols() + sample),
                            value,
                        });
                    }
                }
            }
        }
        Offsets::PerCell(matrix) => matrix.validate_finite("normalization factors")?,
    }
    Ok(())
}

pub(super) fn checked_log2_scale(
    value: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let scaled = std::f64::consts::LOG2_E * value;
    if scaled.is_finite() {
        Ok(scaled)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: scaled,
        })
    }
}

pub(super) fn checked_log2_standard_error(
    covariance: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let se = covariance.max(0.0).sqrt();
    checked_log2_scale(se, index, context)
}

pub(super) fn checked_log2_covariance(
    value: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let log2_e = std::f64::consts::LOG2_E;
    let log2_e_squared = checked_log2_scale(log2_e, index, context)?;
    let scaled = log2_e_squared * value;
    if scaled.is_finite() {
        Ok(scaled)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: scaled,
        })
    }
}

pub(super) fn retain_covariance_on_log2_scale(
    covariance: &DMatrix<f64>,
    destination: Option<&mut [f64]>,
    start_index: usize,
    context: &str,
) -> Result<bool, DeseqError> {
    let p = covariance.ncols();
    let mut positive = true;
    if let Some(destination) = destination {
        for (row, values) in destination.chunks_exact_mut(p).enumerate() {
            for (col, value) in values.iter_mut().enumerate() {
                *value = checked_log2_covariance(
                    covariance[(row, col)],
                    start_index + row * p + col,
                    context,
                )?;
            }
            positive &= values[row] > 0.0;
        }
        return Ok(positive);
    }

    // The caller checks that the matrix is finite. LOG2_E squared is less
    // than four, so only the diagonal needs scaling in the regular range.
    if covariance.iter().any(|value| value.abs() > f64::MAX / 4.0) {
        for row in 0..p {
            for col in 0..p {
                checked_log2_covariance(
                    covariance[(row, col)],
                    start_index + row * p + col,
                    context,
                )?;
            }
        }
    }
    for index in 0..p {
        let diagonal = checked_log2_covariance(
            covariance[(index, index)],
            start_index + index * (p + 1),
            context,
        )?;
        positive &= diagonal > 0.0;
    }
    Ok(positive)
}
