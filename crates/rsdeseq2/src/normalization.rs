use crate::core::CountMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::math::median::median;
use crate::matrix::RowMajorMatrix;
use crate::options::SizeFactorMethod;

/// Estimate size factors with the selected method.
pub fn estimate_size_factors(
    counts: &CountMatrix,
    method: SizeFactorMethod,
) -> Result<Vec<f64>, DeseqError> {
    estimate_size_factors_with_options(counts, method, None, None)
}

/// Estimate size factors with optional geometric means and control genes.
pub fn estimate_size_factors_with_options(
    counts: &CountMatrix,
    method: SizeFactorMethod,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    match method {
        SizeFactorMethod::Ratio => {
            estimate_size_factors_ratio_with_options(counts, geo_means, control_genes)
        }
        SizeFactorMethod::PosCounts => {
            estimate_size_factors_poscounts_with_options(counts, geo_means, control_genes)
        }
    }
}

/// Estimate DESeq2 median-ratio size factors.
pub fn estimate_size_factors_ratio(counts: &CountMatrix) -> Result<Vec<f64>, DeseqError> {
    estimate_size_factors_ratio_with_options(counts, None, None)
}

/// Estimate DESeq2 median-ratio size factors with optional geometric means and controls.
///
/// When `geo_means` is supplied, size factors are stabilized to geometric mean 1,
/// matching DESeq2's frozen-size-factor behavior.
pub fn estimate_size_factors_ratio_with_options(
    counts: &CountMatrix,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    let incoming_geo_means = geo_means.is_some();
    let log_geo_means = match geo_means {
        Some(values) => supplied_log_geo_means(counts, values)?,
        None => ratio_log_geo_means(counts)?,
    };
    estimate_from_log_geo_means(counts, &log_geo_means, control_genes, incoming_geo_means)
}

/// Estimate DESeq2 `poscounts` size factors.
pub fn estimate_size_factors_poscounts(counts: &CountMatrix) -> Result<Vec<f64>, DeseqError> {
    estimate_size_factors_poscounts_with_options(counts, None, None)
}

/// Estimate DESeq2 `poscounts` size factors with optional geometric means and controls.
///
/// DESeq2's `poscounts` method computes each gene's pseudo-geometric mean as
/// the nth root of the product of positive counts, where n is the total sample
/// count. All-zero genes are skipped.
pub fn estimate_size_factors_poscounts_with_options(
    counts: &CountMatrix,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    let incoming_geo_means = geo_means.is_some();
    let log_geo_means = match geo_means {
        Some(values) => supplied_log_geo_means(counts, values)?,
        None => poscounts_log_geo_means(counts)?,
    };
    estimate_from_log_geo_means(counts, &log_geo_means, control_genes, incoming_geo_means)
}

/// Divide counts by sample-specific size factors.
pub fn normalized_counts(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_size_factors(size_factors, counts.n_samples())?;
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        for (sample, count) in counts.row_values(gene).iter().copied().enumerate() {
            values.push(f64::from(count) / size_factors[sample]);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

/// Expand sample-specific size factors into a genes x samples normalization matrix.
///
/// This is equivalent to DESeq2's `getSizeOrNormFactors()` fallback when no
/// gene/sample normalization factors are present.
pub fn normalization_factors_from_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_size_factors(size_factors, counts.n_samples())?;
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _ in 0..counts.n_genes() {
        values.extend_from_slice(size_factors);
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

/// Divide counts by gene/sample normalization factors.
///
/// DESeq2's `counts(dds, normalized=TRUE)` uses `normalizationFactors(dds)`
/// when present and lets those factors preempt size factors.
pub fn normalized_counts_with_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_normalization_factors(counts, normalization_factors)?;
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let count_row = counts.row_values(gene);
        let factor_row = normalization_factors.row(gene)?;
        for (count, factor) in count_row.iter().copied().zip(factor_row.iter().copied()) {
            values.push(f64::from(count) / factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

/// Validate that normalization factors can be used as DESeq2-style count-scale factors.
pub fn validate_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != counts.n_genes()
        || normalization_factors.n_cols() != counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factors".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: normalization_factors.len(),
        });
    }
    for (idx, value) in normalization_factors.as_slice().iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("normalization factor at index {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

/// Calculate per-gene base means from normalized counts.
pub fn base_mean(normalized_counts: &RowMajorMatrix<f64>) -> Result<Vec<f64>, DeseqError> {
    normalized_counts.validate_finite("normalized counts")?;
    let mut means = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        means.push(checked_mean(row).ok_or_else(|| DeseqError::NonFiniteValue {
            context: "baseMean".to_string(),
            index: Some(gene),
            value: f64::NAN,
        })?);
    }
    Ok(means)
}

/// Calculate DESeq2-style weighted base means.
///
/// This mirrors `getBaseMeansAndVariances`: normalized counts are multiplied
/// elementwise by observation weights, then ordinary row means are computed.
pub fn base_mean_with_weights(
    normalized_counts: &RowMajorMatrix<f64>,
    weights: &RowMajorMatrix<f64>,
) -> Result<Vec<f64>, DeseqError> {
    validate_weighted_base_inputs(normalized_counts, weights)?;
    let n_samples = normalized_counts.n_cols();
    let mut means = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let weight_row = weights.row(gene)?;
        let sum = checked_weighted_sum(row, weight_row, "weighted baseMean", gene)?;
        let mean = sum / n_samples as f64;
        if !mean.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "weighted baseMean".to_string(),
                index: Some(gene),
                value: mean,
            });
        }
        means.push(mean);
    }
    Ok(means)
}

/// Calculate per-gene sample variance from normalized counts.
///
/// DESeq2 stores `baseVar` using `matrixStats::rowVars`, which uses sample
/// variance. For one sample, R variance is undefined; this returns `NaN` for
/// each row in that case.
pub fn base_variance(normalized_counts: &RowMajorMatrix<f64>) -> Result<Vec<f64>, DeseqError> {
    normalized_counts.validate_finite("normalized counts")?;
    let n_samples = normalized_counts.n_cols();
    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        if n_samples < 2 {
            variances.push(f64::NAN);
            continue;
        }
        let mean = checked_mean(row).ok_or_else(|| DeseqError::NonFiniteValue {
            context: "baseVar mean".to_string(),
            index: Some(gene),
            value: f64::NAN,
        })?;
        let sum_squares =
            checked_centered_sum_squares(row.iter().copied(), mean).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "baseVar".to_string(),
                    index: Some(gene),
                    value: f64::NAN,
                }
            })?;
        variances.push(sum_squares / (n_samples as f64 - 1.0));
    }
    Ok(variances)
}

/// Calculate DESeq2-style weighted base variances.
///
/// This mirrors `getBaseMeansAndVariances`: normalized counts are multiplied
/// elementwise by observation weights, then ordinary row sample variances are
/// computed with the same shape as `matrixStats::rowVars`.
pub fn base_variance_with_weights(
    normalized_counts: &RowMajorMatrix<f64>,
    weights: &RowMajorMatrix<f64>,
) -> Result<Vec<f64>, DeseqError> {
    validate_weighted_base_inputs(normalized_counts, weights)?;
    let n_samples = normalized_counts.n_cols();
    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let weight_row = weights.row(gene)?;
        if n_samples < 2 {
            variances.push(f64::NAN);
            continue;
        }
        let weighted_values = checked_weighted_values(row, weight_row, "weighted baseVar", gene)?;
        let mean = checked_mean(&weighted_values).ok_or_else(|| DeseqError::NonFiniteValue {
            context: "weighted baseVar mean".to_string(),
            index: Some(gene),
            value: f64::NAN,
        })?;
        let sum_squares = checked_centered_sum_squares(weighted_values.iter().copied(), mean)
            .ok_or_else(|| DeseqError::NonFiniteValue {
                context: "weighted baseVar".to_string(),
                index: Some(gene),
                value: f64::NAN,
            })?;
        variances.push(sum_squares / (n_samples as f64 - 1.0));
    }
    Ok(variances)
}

fn checked_mean(values: &[f64]) -> Option<f64> {
    let sum = checked_sum(values.iter().copied())?;
    let mean = sum / values.len() as f64;
    mean.is_finite().then_some(mean)
}

fn checked_weighted_sum(
    values: &[f64],
    weights: &[f64],
    context: &str,
    gene: usize,
) -> Result<f64, DeseqError> {
    let weighted_values = checked_weighted_values(values, weights, context, gene)?;
    checked_sum(weighted_values).ok_or_else(|| DeseqError::NonFiniteValue {
        context: context.to_string(),
        index: Some(gene),
        value: f64::NAN,
    })
}

fn checked_weighted_values(
    values: &[f64],
    weights: &[f64],
    context: &str,
    gene: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut out = Vec::with_capacity(values.len());
    for (value, weight) in values.iter().copied().zip(weights.iter().copied()) {
        let weighted = value * weight;
        if !weighted.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(gene),
                value: weighted,
            });
        }
        out.push(weighted);
    }
    Ok(out)
}

fn checked_centered_sum_squares(values: impl IntoIterator<Item = f64>, mean: f64) -> Option<f64> {
    let mut sum = 0.0;
    for value in values {
        let centered = value - mean;
        let square = centered * centered;
        let next = checked_sum2(sum, square)?;
        if !centered.is_finite() || !square.is_finite() {
            return None;
        }
        sum = next;
    }
    Some(sum)
}

fn checked_sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    values.into_iter().try_fold(0.0, checked_sum2)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn validate_weighted_base_inputs(
    normalized_counts: &RowMajorMatrix<f64>,
    weights: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    normalized_counts.validate_finite("normalized counts")?;
    if weights.n_rows() != normalized_counts.n_rows()
        || weights.n_cols() != normalized_counts.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "base metadata weights".to_string(),
            expected: normalized_counts.len(),
            actual: weights.len(),
        });
    }
    for (idx, value) in weights.as_slice().iter().copied().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "base metadata weight".to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

fn ratio_log_geo_means(counts: &CountMatrix) -> Result<Vec<f64>, DeseqError> {
    (0..counts.n_genes())
        .map(|gene| {
            let row = counts.row_values(gene);
            if row.contains(&0) {
                Ok(f64::NEG_INFINITY)
            } else {
                let sum = checked_sum(row.iter().map(|count| f64::from(*count).ln())).ok_or_else(
                    || DeseqError::NonFiniteValue {
                        context: "ratio geometric mean log sum".to_string(),
                        index: Some(gene),
                        value: f64::NAN,
                    },
                )?;
                Ok(sum / counts.n_samples() as f64)
            }
        })
        .collect()
}

fn poscounts_log_geo_means(counts: &CountMatrix) -> Result<Vec<f64>, DeseqError> {
    (0..counts.n_genes())
        .map(|gene| {
            let row = counts.row_values(gene);
            if row.iter().all(|count| *count == 0) {
                Ok(f64::NEG_INFINITY)
            } else {
                let sum = checked_sum(
                    row.iter()
                        .filter(|count| **count > 0)
                        .map(|count| f64::from(*count).ln()),
                )
                .ok_or_else(|| DeseqError::NonFiniteValue {
                    context: "poscounts geometric mean log sum".to_string(),
                    index: Some(gene),
                    value: f64::NAN,
                })?;
                Ok(sum / counts.n_samples() as f64)
            }
        })
        .collect()
}

fn supplied_log_geo_means(counts: &CountMatrix, geo_means: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if geo_means.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "geometric means",
            counts.n_genes(),
            geo_means.len(),
        ));
    }
    geo_means
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, value)| {
            if value.is_nan() || value < 0.0 {
                Err(DeseqError::InvalidSizeFactors {
                    reason: format!("geometric mean at row {idx} must be non-negative or infinite"),
                })
            } else if value == 0.0 {
                Ok(f64::NEG_INFINITY)
            } else {
                Ok(value.ln())
            }
        })
        .collect()
}

fn estimate_from_log_geo_means(
    counts: &CountMatrix,
    log_geo_means: &[f64],
    control_genes: Option<&[usize]>,
    stabilize_to_geometric_mean_one: bool,
) -> Result<Vec<f64>, DeseqError> {
    if log_geo_means.iter().all(|value| value.is_infinite()) {
        return Err(DeseqError::NoUsableGenesForSizeFactors);
    }
    let genes = usable_gene_indices(counts, control_genes)?;
    let mut size_factors = Vec::with_capacity(counts.n_samples());
    for sample in 0..counts.n_samples() {
        let mut log_ratios = Vec::new();
        for gene in &genes {
            let log_geo_mean = log_geo_means[*gene];
            let count = counts.row_values(*gene)[sample];
            if log_geo_mean.is_finite() && count > 0 {
                log_ratios.push(f64::from(count).ln() - log_geo_mean);
            }
        }
        if log_ratios.is_empty() {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("sample {sample} has no usable positive count ratios"),
            });
        }
        size_factors.push(median(&log_ratios)?.exp());
    }
    validate_size_factors(&size_factors, counts.n_samples())?;
    if stabilize_to_geometric_mean_one {
        stabilize_size_factors(&mut size_factors)?;
    }
    Ok(size_factors)
}

fn usable_gene_indices(
    counts: &CountMatrix,
    control_genes: Option<&[usize]>,
) -> Result<Vec<usize>, DeseqError> {
    match control_genes {
        Some(indices) => {
            if indices.is_empty() {
                return Err(DeseqError::NoUsableGenesForSizeFactors);
            }
            for index in indices {
                if *index >= counts.n_genes() {
                    return Err(DeseqError::InvalidDimensions {
                        context: "control gene index".to_string(),
                        expected: counts.n_genes().saturating_sub(1),
                        actual: *index,
                    });
                }
            }
            Ok(indices.to_vec())
        }
        None => Ok((0..counts.n_genes()).collect()),
    }
}

fn stabilize_size_factors(size_factors: &mut [f64]) -> Result<(), DeseqError> {
    let log_sum = checked_sum(size_factors.iter().map(|value| value.ln())).ok_or_else(|| {
        DeseqError::InvalidSizeFactors {
            reason: "cannot stabilize size factors to geometric mean one".to_string(),
        }
    })?;
    let log_mean = log_sum / size_factors.len() as f64;
    let scale = log_mean.exp();
    if !scale.is_finite() || scale <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "cannot stabilize size factors to geometric mean one".to_string(),
        });
    }
    for value in size_factors {
        *value /= scale;
    }
    Ok(())
}

fn validate_size_factors(size_factors: &[f64], n_samples: usize) -> Result<(), DeseqError> {
    if size_factors.len() != n_samples {
        return Err(invalid_dimensions(
            "size factors",
            n_samples,
            size_factors.len(),
        ));
    }
    for (idx, value) in size_factors.iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("size factor at sample {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}
