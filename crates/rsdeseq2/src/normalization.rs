use crate::errors::{ErrorKind, Requirement, invalid_parameter};
use crate::matrix::try_alloc_matrix;
use crate::numeric::checked::{self, CheckedExt};
use crate::numeric::math::ProfileMath;
mod offsets;
pub use offsets::normalization_factors_with_offsets;
pub(crate) use offsets::sample_factors_with_offsets;
#[doc(hidden)]
pub use offsets::visit_offset_factors_column_major;

use crate::errors::{DeseqError, invalid_dimensions};
use crate::math::median::median_in_place;
use crate::matrix::{CountMatrix, CountView, MatrixOrder, MatrixView, RowMajorMatrix};
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
    estimate_size_factors_for_counts(counts, method, geo_means, control_genes)
}

/// Estimate median-ratio size factors from validated counts in either layout.
#[doc(hidden)]
pub fn estimate_size_factors_from_view(
    counts: CountView<'_>,
    method: SizeFactorMethod,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    estimate_size_factors_for_counts(counts, method, geo_means, control_genes)
}

fn estimate_size_factors_for_counts(
    counts: impl SizeFactorCounts,
    method: SizeFactorMethod,
    geo_means: Option<&[f64]>,
    control_genes: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    let log_geo_means = match (geo_means, method) {
        (Some(values), _) => supplied_log_geo_means(counts.n_genes(), values)?,
        (None, SizeFactorMethod::Ratio) => ratio_log_geo_means(counts)?,
        (None, SizeFactorMethod::PosCounts) => poscounts_log_geo_means(counts)?,
    };
    estimate_from_log_geo_means(counts, &log_geo_means, control_genes, geo_means.is_some())
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
    estimate_size_factors_with_options(counts, SizeFactorMethod::Ratio, geo_means, control_genes)
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
    estimate_size_factors_with_options(
        counts,
        SizeFactorMethod::PosCounts,
        geo_means,
        control_genes,
    )
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
            let log_geo_means = supplied_log_geo_means(counts.n_genes(), &geo_means)?;
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
    normalized_counts_for_rows(counts, Offsets::PerSample(size_factors), None)
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
    normalized_counts_for_rows(counts, Offsets::PerCell(normalization_factors), None)
}

pub(crate) fn normalized_counts_for_rows(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    selected: Option<&[usize]>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let genes = selected.map_or(counts.n_genes(), <[usize]>::len);
    let samples = counts.n_samples();
    let context = match offsets {
        Offsets::PerSample(factors) => {
            validate_size_factors(factors, samples)?;
            "normalized count"
        }
        Offsets::PerCell(factors) => {
            if let Some(selected) = selected {
                if factors.n_cols() != samples {
                    return Err(invalid_dimensions(
                        "normalization factors",
                        genes * samples,
                        genes.saturating_mul(factors.n_cols()),
                    ));
                }
                for (gene, &source) in selected.iter().enumerate() {
                    for (sample, &factor) in factors.row(source)?.iter().enumerate() {
                        validate_normalization_factor_value(factor, gene * samples + sample)?;
                    }
                }
            } else {
                validate_normalization_factors(counts, factors)?;
            }
            "normalization-factor normalized count"
        }
    };
    let mut values = try_alloc_matrix(genes, samples, "normalized counts")?;
    for gene in 0..genes {
        let source = selected.map_or(gene, |selected| selected[gene]);
        let factors = offsets.row(source)?;
        for (sample, (&count, &factor)) in counts.row(source)?.iter().zip(factors).enumerate() {
            values.push(checked::try_div(f64::from(count), factor).map_err(|error| {
                normalized_count_error(context, gene * samples + sample, error.value)
            })?);
        }
    }
    RowMajorMatrix::from_row_major(genes, samples, values)
}

/// Visit normalized counts in column-major order for foreign output buffers.
/// Cell-specific factors take precedence. The caller discards output on error.
#[doc(hidden)]
pub fn visit_normalized_counts_column_major(
    counts: CountView<'_>,
    size_factors: &[f64],
    normalization_factors: Option<MatrixView<'_, f64>>,
    write: impl FnMut(usize, f64),
) -> Result<(), DeseqError> {
    let context = validate_count_factors(counts, size_factors, normalization_factors)?;
    if let Some(factors) = normalization_factors {
        visit_normalized_columns(
            counts,
            |gene, sample| *factors.at(gene, sample),
            context,
            write,
        )
    } else {
        visit_normalized_columns(counts, |_, sample| size_factors[sample], context, write)
    }
}

fn validate_count_factors(
    counts: CountView<'_>,
    size_factors: &[f64],
    normalization_factors: Option<MatrixView<'_, f64>>,
) -> Result<&'static str, DeseqError> {
    if let Some(factors) = normalization_factors {
        validate_normalization_factor_view(
            factors,
            counts.n_genes(),
            counts.n_samples(),
            "normalization factors",
        )?;
        Ok("normalization-factor normalized count")
    } else {
        validate_size_factors(size_factors, counts.n_samples())?;
        Ok("normalized count")
    }
}

fn visit_normalized_columns(
    counts: CountView<'_>,
    factor: impl Fn(usize, usize) -> f64,
    context: &str,
    mut write: impl FnMut(usize, f64),
) -> Result<(), DeseqError> {
    let genes = counts.n_genes();
    let samples = counts.n_samples();
    let mut first_error = None;
    let mut genes_to_check = genes;
    for sample in 0..samples {
        for gene in 0..genes_to_check {
            match checked::try_div(f64::from(counts.at(gene, sample)), factor(gene, sample)) {
                Ok(value) => write(sample * genes + gene, value),
                Err(error) => {
                    first_error = Some((gene * samples + sample, error.value));
                    // In later columns, only an earlier row can take precedence.
                    genes_to_check = gene;
                    break;
                }
            }
        }
    }
    if let Some((index, value)) = first_error {
        return Err(normalized_count_error(context, index, value));
    }
    Ok(())
}

fn normalized_count_error(context: &str, index: usize, value: f64) -> DeseqError {
    DeseqError::NonFiniteValue {
        context: context.to_string(),
        index: Some(index),
        value,
    }
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
    validate_normalization_factor_view(normalization_factors.into(), n_genes, n_samples, context)
}

fn validate_normalization_factor_view(
    normalization_factors: MatrixView<'_, f64>,
    n_genes: usize,
    n_samples: usize,
    context: &str,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != n_genes || normalization_factors.n_cols() != n_samples {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: n_genes * n_samples,
            actual: normalization_factors.as_slice().len(),
        });
    }
    if let Some((index, value)) = first_invalid_cell(normalization_factors, |value| {
        !value.is_finite() || value <= 0.0
    }) {
        return validate_normalization_factor_value(value, index);
    }
    Ok(())
}

fn first_invalid_cell(
    matrix: MatrixView<'_, f64>,
    invalid: impl Fn(f64) -> bool,
) -> Option<(usize, f64)> {
    // Scan contiguous storage first. Only invalid column-major inputs need
    // a second scan to find the earliest error in logical row order.
    let (storage_index, &value) = matrix
        .as_slice()
        .iter()
        .enumerate()
        .find(|(_, value)| invalid(**value))?;
    if matrix.order() == MatrixOrder::RowMajor {
        return Some((storage_index, value));
    }
    for gene in 0..matrix.n_rows() {
        for sample in 0..matrix.n_cols() {
            let value = *matrix.at(gene, sample);
            if invalid(value) {
                return Some((gene * matrix.n_cols() + sample, value));
            }
        }
    }
    None
}

/// Reject a normalization factor that is not finite and positive.
#[inline]
pub(crate) fn validate_normalization_factor_value(
    value: f64,
    index: usize,
) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Normalization,
            format!("normalization factor at index {index}"),
            Requirement::FiniteAndPositive,
        ));
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
    let mut weighted_values = Vec::with_capacity(normalized_counts.n_cols());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let weight_row = weights.row(gene)?;
        let context = "weighted baseMean";
        checked_weighted_values(row, weight_row, context, gene, &mut weighted_values)?;
        let mean = checked::try_mean(&weighted_values).nonfinite_at(context, Some(gene))?;
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
    let mut weighted_values = Vec::with_capacity(n_samples);
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let weight_row = weights.row(gene)?;
        if n_samples < 2 {
            variances.push(f64::NAN);
            continue;
        }
        checked_weighted_values(
            row,
            weight_row,
            "weighted baseVar",
            gene,
            &mut weighted_values,
        )?;
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

/// Compute base metadata from borrowed counts with one row of scratch.
/// The caller discards emitted rows if any input or calculation fails.
#[doc(hidden)]
pub fn visit_base_metadata(
    counts: CountView<'_>,
    size_factors: &[f64],
    factors: Option<MatrixView<'_, f64>>,
    weights: Option<MatrixView<'_, f64>>,
    mut write: impl FnMut(usize, f64, f64, bool),
) -> Result<(), DeseqError> {
    let context = validate_count_factors(counts, size_factors, factors)?;
    let genes = counts.n_genes();
    let samples = counts.n_samples();
    let mut values = try_alloc_matrix(samples, 1, "base metadata row")?;
    let factor = |gene, sample| match factors {
        Some(factors) => *factors.at(gene, sample),
        None => size_factors[sample],
    };
    // Keep error precedence across rows: normalization, weight validity,
    // base means, and then base variances.
    let mut mean_error =
        weights.and_then(|weights| validate_base_weights(weights, genes, samples).err());
    let mut variance_error = None;
    for gene in 0..genes {
        values.clear();
        let mut all_zero = true;
        for sample in 0..samples {
            let count = counts.at(gene, sample);
            all_zero &= count == 0;
            values.push(
                checked::try_div(f64::from(count), factor(gene, sample)).map_err(|error| {
                    normalized_count_error(context, gene * samples + sample, error.value)
                })?,
            );
        }
        if mean_error.is_some() {
            continue;
        }
        let mean = match base_row_mean(&mut values, weights, gene) {
            Ok(mean) => mean,
            Err(error) => {
                mean_error = Some(error);
                continue;
            }
        };
        if variance_error.is_some() {
            continue;
        }
        let variance = if samples == 1 {
            Ok(f64::NAN)
        } else {
            let context = if weights.is_some() {
                "weighted baseVar"
            } else {
                "baseVar"
            };
            checked_sample_variance(&values, mean).nonfinite_at(context, Some(gene))
        };
        match variance {
            Ok(variance) => write(gene, mean, variance, all_zero),
            Err(error) => variance_error = Some(error),
        }
    }
    mean_error.or(variance_error).map_or(Ok(()), Err)
}

fn base_row_mean(
    values: &mut [f64],
    weights: Option<MatrixView<'_, f64>>,
    gene: usize,
) -> Result<f64, DeseqError> {
    let context = if let Some(weights) = weights {
        for (sample, value) in values.iter_mut().enumerate() {
            *value = checked::try_mul(*value, *weights.at(gene, sample))
                .nonfinite_at("weighted baseMean", Some(gene))?;
        }
        "weighted baseMean"
    } else {
        "baseMean"
    };
    checked::try_mean(values).nonfinite_at(context, Some(gene))
}

fn checked_weighted_values(
    values: &[f64],
    weights: &[f64],
    context: &str,
    gene: usize,
    out: &mut Vec<f64>,
) -> Result<(), DeseqError> {
    out.clear();
    for (value, weight) in values.iter().copied().zip(weights.iter().copied()) {
        let weighted = checked::try_mul(value, weight).nonfinite_at(context, Some(gene))?;
        out.push(weighted);
    }
    Ok(())
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
    validate_base_weights(
        weights.into(),
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
    )
}

fn validate_base_weights(
    weights: MatrixView<'_, f64>,
    genes: usize,
    samples: usize,
) -> Result<(), DeseqError> {
    if weights.n_rows() != genes || weights.n_cols() != samples {
        return Err(DeseqError::InvalidDimensions {
            context: "base metadata weights".to_string(),
            expected: genes * samples,
            actual: weights.as_slice().len(),
        });
    }
    if let Some((index, value)) =
        first_invalid_cell(weights, |value| !value.is_finite() || value < 0.0)
    {
        return Err(DeseqError::NonFiniteValue {
            context: "base metadata weight".to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

// Keep native row slices available to the owned-count path while sharing
// geometric means and the median pass with strided foreign input.
trait SizeFactorCounts: Copy {
    fn n_genes(self) -> usize;
    fn n_samples(self) -> usize;
    fn row(self, gene: usize) -> impl Iterator<Item = u32>;
    fn at(self, gene: usize, sample: usize) -> u32;

    fn row_contains_zero(self, gene: usize) -> bool {
        self.row(gene).any(|count| count == 0)
    }
}

impl SizeFactorCounts for &CountMatrix {
    fn n_genes(self) -> usize {
        CountMatrix::n_genes(self)
    }
    fn n_samples(self) -> usize {
        CountMatrix::n_samples(self)
    }
    fn row(self, gene: usize) -> impl Iterator<Item = u32> {
        self.row_values(gene).iter().copied()
    }
    fn at(self, gene: usize, sample: usize) -> u32 {
        self.row_values(gene)[sample]
    }
    fn row_contains_zero(self, gene: usize) -> bool {
        self.row_values(gene).contains(&0)
    }
}

impl SizeFactorCounts for CountView<'_> {
    fn n_genes(self) -> usize {
        CountView::n_genes(&self)
    }
    fn n_samples(self) -> usize {
        CountView::n_samples(&self)
    }
    fn row(self, gene: usize) -> impl Iterator<Item = u32> {
        (0..self.n_samples()).map(move |sample| self.at(gene, sample))
    }
    fn at(self, gene: usize, sample: usize) -> u32 {
        CountView::at(&self, gene, sample)
    }
}

fn ratio_log_geo_means(counts: impl SizeFactorCounts) -> Result<Vec<f64>, DeseqError> {
    let samples = counts.n_samples();
    (0..counts.n_genes())
        .map(|gene| {
            if counts.row_contains_zero(gene) {
                Ok(f64::NEG_INFINITY)
            } else {
                // R's rowMeans() accumulates with extended precision. A
                // compensated sum preserves the same low-order information in
                // portable f64 arithmetic and prevents near-boundary genes
                // from changing optimizer routes through size-factor drift.
                let sum =
                    checked_compensated_sum(counts.row(gene).map(|count| f64::from(count).ln_p()))
                        .nonfinite_at("ratio geometric mean log sum", Some(gene))?;
                checked::try_div(sum, samples as f64).map_err(|error| DeseqError::NonFiniteValue {
                    context: "ratio geometric mean log mean".to_string(),
                    index: Some(gene),
                    value: error.value,
                })
            }
        })
        .collect()
}

fn poscounts_log_geo_means(counts: impl SizeFactorCounts) -> Result<Vec<f64>, DeseqError> {
    let samples = counts.n_samples();
    (0..counts.n_genes())
        .map(|gene| {
            if counts.row(gene).all(|count| count == 0) {
                Ok(f64::NEG_INFINITY)
            } else {
                let sum = checked_compensated_sum(
                    counts
                        .row(gene)
                        .filter(|&count| count > 0)
                        .map(|count| f64::from(count).ln_p()),
                )
                .nonfinite_at("poscounts geometric mean log sum", Some(gene))?;
                checked::try_div(sum, samples as f64).map_err(|error| DeseqError::NonFiniteValue {
                    context: "poscounts geometric mean log mean".to_string(),
                    index: Some(gene),
                    value: error.value,
                })
            }
        })
        .collect()
}

fn supplied_log_geo_means(genes: usize, geo_means: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if geo_means.len() != genes {
        return Err(invalid_dimensions(
            "geometric means",
            genes,
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
    counts: impl SizeFactorCounts,
    log_geo_means: &[f64],
    control_genes: Option<&[usize]>,
    stabilize_to_geometric_mean_one: bool,
) -> Result<Vec<f64>, DeseqError> {
    if log_geo_means.iter().all(|value| value.is_infinite()) {
        return Err(DeseqError::NoUsableGenesForSizeFactors);
    }
    let mut genes = usable_gene_indices(counts.n_genes(), control_genes)?;
    genes.retain(|&gene| log_geo_means[gene].is_finite());
    let mut size_factors = Vec::with_capacity(counts.n_samples());
    let mut log_ratios = try_alloc_matrix(genes.len(), 1, "size-factor log ratios")?;
    for sample in 0..counts.n_samples() {
        log_ratios.clear();
        for gene in &genes {
            let log_geo_mean = log_geo_means[*gene];
            let count = counts.at(*gene, sample);
            if count > 0 {
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
        size_factors.push(median_in_place(&mut log_ratios)?.exp_p());
    }
    validate_size_factors(&size_factors, counts.n_samples())?;
    if stabilize_to_geometric_mean_one {
        stabilize_size_factors(&mut size_factors)?;
    }
    Ok(size_factors)
}

fn usable_gene_indices(
    genes: usize,
    control_genes: Option<&[usize]>,
) -> Result<Vec<usize>, DeseqError> {
    match control_genes {
        Some(indices) => {
            if indices.is_empty() {
                return Err(DeseqError::NoUsableGenesForSizeFactors);
            }
            for index in indices {
                if *index >= genes {
                    return Err(DeseqError::InvalidDimensions {
                        context: "control gene index".to_string(),
                        expected: genes.saturating_sub(1),
                        actual: *index,
                    });
                }
            }
            Ok(indices.to_vec())
        }
        None => Ok((0..genes).collect()),
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
            return Err(invalid_parameter(
                ErrorKind::Normalization,
                format!("size factor at sample {idx}"),
                Requirement::FiniteAndPositive,
            ));
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
