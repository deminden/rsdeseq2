use crate::matrix::try_alloc_matrix;
use crate::numeric::checked::{self, CheckedExt};
use crate::numeric::math::ProfileMath;
mod offsets;
pub use offsets::normalization_factors_with_offsets;

use crate::errors::{DeseqError, invalid_dimensions};
use crate::math::median::median;
use crate::matrix::CountMatrix;
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

/// Estimate size factors the way `estimateSizeFactors(dds, type = ...)` and
/// `DESeq(dds, sfType = ...)` do, rather than `estimateSizeFactorsForMatrix`.
///
/// The two differ only for `poscounts` without supplied geometric means:
/// the `DESeqDataSet` method computes the non-zero geometric means itself and
/// passes them in, so its size factors are rescaled to a geometric mean of 1,
/// while the matrix function's are not. Supplied geometric means are used as
/// given (DESeq2's method would replace them for `poscounts`).
pub fn estimate_size_factors_for_dataset(
    counts: &CountMatrix,
    method: SizeFactorMethod,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    match (method, geo_means) {
        (SizeFactorMethod::PosCounts, None) => {
            // geoMeanNZ, then log(geoMeans) as estimateSizeFactorsForMatrix
            // takes it: exp and log are applied as in R.
            let geo_means = poscounts_log_geo_means(counts)?
                .into_iter()
                .map(f64::exp_p)
                .collect::<Vec<_>>();
            let log_geo_means = supplied_log_geo_means(counts, &geo_means)?;
            estimate_from_log_geo_means(counts, &log_geo_means, control_genes, true)
        }
        _ => estimate_size_factors_with_options(counts, method, geo_means, control_genes),
    }
}

/// Divide counts by sample-specific size factors.
pub fn normalized_counts(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_size_factors(size_factors, counts.n_samples())?;
    let mut values = try_alloc_matrix(counts.n_genes(), counts.n_samples(), "normalized counts")?;
    for gene in 0..counts.n_genes() {
        for (sample, count) in counts.row_values(gene).iter().copied().enumerate() {
            values.push(
                checked::try_div(f64::from(count), size_factors[sample]).map_err(|error| {
                    DeseqError::NonFiniteValue {
                        context: "normalized count".to_string(),
                        index: Some(gene * counts.n_samples() + sample),
                        value: error.value,
                    }
                })?,
            );
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
    let mut values = try_alloc_matrix(counts.n_genes(), counts.n_samples(), "normalized counts")?;
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
    let mut values = try_alloc_matrix(counts.n_genes(), counts.n_samples(), "normalized counts")?;
    for gene in 0..counts.n_genes() {
        let count_row = counts.row_values(gene);
        let factor_row = normalization_factors.row(gene)?;
        for (sample, (count, factor)) in count_row
            .iter()
            .copied()
            .zip(factor_row.iter().copied())
            .enumerate()
        {
            values.push(checked::try_div(f64::from(count), factor).map_err(|error| {
                DeseqError::NonFiniteValue {
                    context: "normalization-factor normalized count".to_string(),
                    index: Some(gene * counts.n_samples() + sample),
                    value: error.value,
                }
            })?);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

/// Validate that normalization factors can be used as DESeq2-style count-scale factors.
pub fn validate_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    validate_normalization_factor_matrix(
        normalization_factors,
        counts.n_genes(),
        counts.n_samples(),
        "normalization factors",
    )
}

/// Check a genes x samples normalization-factor matrix: the shape (reported
/// under `context`) and then every value in row-major order.
pub(crate) fn validate_normalization_factor_matrix(
    normalization_factors: &RowMajorMatrix<f64>,
    n_genes: usize,
    n_samples: usize,
    context: &str,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != n_genes || normalization_factors.n_cols() != n_samples {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: n_genes * n_samples,
            actual: normalization_factors.len(),
        });
    }
    for (idx, value) in normalization_factors.as_slice().iter().copied().enumerate() {
        validate_normalization_factor_value(value, idx)?;
    }
    Ok(())
}

/// Reject a normalization factor that is not finite and positive.
#[inline]
pub(crate) fn validate_normalization_factor_value(
    value: f64,
    index: usize,
) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: format!("normalization factor at index {index} must be finite and positive"),
        });
    }
    Ok(())
}

/// Calculate per-gene base means from normalized counts.
pub fn base_mean(normalized_counts: &RowMajorMatrix<f64>) -> Result<Vec<f64>, DeseqError> {
    normalized_counts.validate_finite("normalized counts")?;
    let mut means = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        means.push(checked::try_mean(row).nonfinite_at("baseMean", Some(gene))?);
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
    let mut means = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let weight_row = weights.row(gene)?;
        let mean = checked_weighted_mean(row, weight_row, "weighted baseMean", gene)?;
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
        let mean = checked::try_mean(row).nonfinite_at("baseVar mean", Some(gene))?;
        variances.push(checked_sample_variance(row, mean).map_err(|error| {
            DeseqError::NonFiniteValue {
                context: "baseVar".to_string(),
                index: Some(gene),
                value: error.value,
            }
        })?);
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
        let mean = checked::try_mean(&weighted_values)
            .nonfinite_at("weighted baseVar mean", Some(gene))?;
        variances.push(
            checked_sample_variance(&weighted_values, mean).map_err(|error| {
                DeseqError::NonFiniteValue {
                    context: "weighted baseVar".to_string(),
                    index: Some(gene),
                    value: error.value,
                }
            })?,
        );
    }
    Ok(variances)
}

fn checked_weighted_mean(
    values: &[f64],
    weights: &[f64],
    context: &str,
    gene: usize,
) -> Result<f64, DeseqError> {
    let weighted_values = checked_weighted_values(values, weights, context, gene)?;
    checked::try_mean(&weighted_values).nonfinite_at(context, Some(gene))
}

fn checked_weighted_values(
    values: &[f64],
    weights: &[f64],
    context: &str,
    gene: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut out = Vec::with_capacity(values.len());
    for (value, weight) in values.iter().copied().zip(weights.iter().copied()) {
        let weighted = checked::try_mul(value, weight).nonfinite_at(context, Some(gene))?;
        out.push(weighted);
    }
    Ok(out)
}

fn checked_sample_variance(values: &[f64], mean: f64) -> checked::Checked {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        let centered = checked::try_sub(value, mean)?;
        scale = scale.max(centered.abs());
    }
    if scale == 0.0 {
        return Ok(0.0);
    }
    let mut normalized_sum_squares = 0.0;
    for value in values.iter().copied() {
        let centered = checked::try_div(checked::try_sub(value, mean)?, scale)?;
        let square = checked::try_mul(centered, centered)?;
        normalized_sum_squares = checked::try_add(normalized_sum_squares, square)?;
    }
    let denominator = checked::try_sub(values.len() as f64, 1.0)?;
    let normalized_variance = checked::try_div(normalized_sum_squares, denominator)?;
    let scaled_deviation = checked::try_mul(normalized_variance.sqrt(), scale)?;
    checked::try_mul(scaled_deviation, scaled_deviation)
}

fn checked_compensated_sum(values: impl IntoIterator<Item = f64>) -> checked::Checked {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        if !value.is_finite() {
            return checked::reject(value);
        }
        let next = checked::try_add(sum, value)?;
        let roundoff = if sum.abs() >= value.abs() {
            checked::try_add(sum - next, value)?
        } else {
            checked::try_add(value - next, sum)?
        };
        correction = checked::try_add(correction, roundoff)?;
        sum = next;
    }
    checked::try_add(sum, correction)
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
                // R's rowMeans() accumulates with extended precision. A
                // compensated sum preserves the same low-order information in
                // portable f64 arithmetic and prevents near-boundary genes
                // from changing optimizer routes through size-factor drift.
                let sum = checked_compensated_sum(row.iter().map(|count| f64::from(*count).ln_p()))
                    .nonfinite_at("ratio geometric mean log sum", Some(gene))?;
                checked::try_div(sum, counts.n_samples() as f64).map_err(|error| {
                    DeseqError::NonFiniteValue {
                        context: "ratio geometric mean log mean".to_string(),
                        index: Some(gene),
                        value: error.value,
                    }
                })
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
                let sum = checked_compensated_sum(
                    row.iter()
                        .filter(|count| **count > 0)
                        .map(|count| f64::from(*count).ln_p()),
                )
                .nonfinite_at("poscounts geometric mean log sum", Some(gene))?;
                checked::try_div(sum, counts.n_samples() as f64).map_err(|error| {
                    DeseqError::NonFiniteValue {
                        context: "poscounts geometric mean log mean".to_string(),
                        index: Some(gene),
                        value: error.value,
                    }
                })
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
                Ok(value.ln_p())
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
                log_ratios.push(
                    checked::try_sub(f64::from(count).ln_p(), log_geo_mean).map_err(|error| {
                        DeseqError::NonFiniteValue {
                            context: "size-factor log ratio".to_string(),
                            index: Some(*gene),
                            value: error.value,
                        }
                    })?,
                );
            }
        }
        if log_ratios.is_empty() {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("sample {sample} has no usable positive count ratios"),
            });
        }
        size_factors.push(median(&log_ratios)?.exp_p());
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
    let log_sum =
        checked_compensated_sum(size_factors.iter().map(|value| value.ln_p())).map_err(|_| {
            DeseqError::InvalidSizeFactors {
                reason: "cannot stabilize size factors to geometric mean one".to_string(),
            }
        })?;
    let log_mean = checked::div(log_sum, size_factors.len() as f64).ok_or_else(|| {
        DeseqError::InvalidSizeFactors {
            reason: "cannot stabilize size factors to geometric mean one".to_string(),
        }
    })?;
    let scale = log_mean.exp_p();
    if !scale.is_finite() || scale <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "cannot stabilize size factors to geometric mean one".to_string(),
        });
    }
    for value in size_factors {
        *value = checked::div(*value, scale).ok_or_else(|| DeseqError::InvalidSizeFactors {
            reason: "cannot stabilize size factors to geometric mean one".to_string(),
        })?;
    }
    Ok(())
}

pub(crate) fn validate_size_factors(
    size_factors: &[f64],
    n_samples: usize,
) -> Result<(), DeseqError> {
    if size_factors.len() != n_samples {
        return Err(invalid_dimensions(
            "size factors",
            n_samples,
            size_factors.len(),
        ));
    }
    validate_size_factor_values(size_factors)
}

/// Reject the first size factor that is not finite and positive.
pub(crate) fn validate_size_factor_values(size_factors: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in size_factors.iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("size factor at sample {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::checked_compensated_sum;

    #[test]
    fn compensated_sum_preserves_representable_low_order_terms() {
        let values = [1.0e16, 1.0, 1.0];
        let naive = values.into_iter().sum::<f64>();
        let compensated = checked_compensated_sum(values).ok().unwrap();

        assert_eq!(naive, 1.0e16);
        assert_eq!(compensated, 1.0e16 + 2.0);
    }

    #[test]
    fn compensated_sum_rejects_non_finite_values() {
        assert_eq!(checked_compensated_sum([1.0, f64::INFINITY]).ok(), None);
    }
}

/// Normalization offsets for a fit: per-sample size factors, shared by every
/// gene, or per-cell normalization factors.
///
/// `PerSample` lets the GLM fits read the size factors directly instead of
/// expanding them into a genes × samples matrix first.
#[derive(Clone, Copy, Debug)]
pub enum Offsets<'a> {
    /// One factor per sample.
    PerSample(&'a [f64]),
    /// One factor per gene and sample.
    PerCell(&'a RowMajorMatrix<f64>),
}

impl<'a> Offsets<'a> {
    /// Factors for one gene, in sample order.
    #[inline]
    pub(crate) fn row(&self, gene: usize) -> Result<&'a [f64], DeseqError> {
        match *self {
            Self::PerSample(factors) => Ok(factors),
            Self::PerCell(matrix) => matrix.row(gene),
        }
    }
}
