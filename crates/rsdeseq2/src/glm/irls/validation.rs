use super::*;
// Input validation, normalization factors and output scaling.

/// Size factors repeated per gene, rejecting non-positive values as
/// `NonFiniteValue` (shared by the IRLS and beta-prior fits).
pub(in crate::glm) fn normalization_factors_from_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _ in 0..counts.n_genes() {
        for (idx, factor) in size_factors.iter().copied().enumerate() {
            validate_positive_finite(factor, "size factor", idx)?;
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

/// The checks `normalization_factors_from_size_factors` applied before the
/// IRLS fit, without building the matrix.
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

pub(super) fn validate_nf_irls_inputs(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: Offsets<'_>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: &IrlsOptions,
) -> Result<(), DeseqError> {
    if design.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "design rows",
            counts.n_samples(),
            design.n_samples(),
        ));
    }
    let (rows, cols, len) = match normalization_factors {
        // Size factors stand for the genes x samples matrix they used to fill.
        Offsets::PerSample(factors) => (
            counts.n_genes(),
            factors.len(),
            counts.n_genes() * factors.len(),
        ),
        Offsets::PerCell(matrix) => (matrix.n_rows(), matrix.n_cols(), matrix.len()),
    };
    if rows != counts.n_genes() || cols != counts.n_samples() {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factors".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: len,
        });
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    for (idx, dispersion) in dispersions.iter().copied().enumerate() {
        validate_positive_finite(dispersion, "dispersion", idx)?;
    }
    if let Some(weights) = weights {
        if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
            return Err(DeseqError::InvalidDimensions {
                context: "weights".to_string(),
                expected: counts.n_genes() * counts.n_samples(),
                actual: weights.len(),
            });
        }
        for (idx, weight) in weights.as_slice().iter().copied().enumerate() {
            validate_nonnegative_finite(weight, "weight", idx)?;
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
    match normalization_factors {
        // Already checked finite and positive by validate_size_factor_offsets.
        Offsets::PerSample(_) => {}
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
