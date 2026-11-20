use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::errors::{ErrorKind, missing_prerequisite};
use crate::matrix::CountMatrix;
use crate::matrix::RowMajorMatrix;
use crate::matrix::{try_alloc_matrix, try_filled_vec};
use crate::normalization::{
    base_mean, base_variance, normalized_counts_with_factors, validate_normalization_factors,
    validate_size_factors,
};
use crate::numeric::checked::{self, CheckedExt};
use rayon::prelude::*;

const MIN_COOKS_DISPERSION: f64 = 0.04;

/// Cook's distance outputs matching DESeq2's Wald/LRT diagnostic shape.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksOutput {
    /// Per-gene, per-sample Cook's distances.
    pub cooks: RowMajorMatrix<f64>,
    /// Per-gene maximum Cook's distance over eligible samples.
    pub max_cooks: Vec<Option<f64>>,
    /// Robust method-of-moments dispersions used in the Cook's variance.
    pub robust_dispersion: Vec<f64>,
    /// Samples considered when recording `maxCooks`.
    pub samples_for_cooks: Vec<bool>,
}

/// Options for DESeq2-style Cook's outlier count replacement.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksReplacementOptions {
    /// Fraction to trim from each side of normalized counts before taking the mean.
    pub trim: f64,
    /// Cook's distance cutoff. DESeq2 defaults to `qf(.99, p, m - p)`.
    pub cooks_cutoff: f64,
    /// Minimum model-matrix cell size for a sample to be replaceable.
    pub min_replicates: usize,
    /// Optional explicit replacement-eligible samples, matching `whichSamples`.
    pub which_samples: Option<Vec<bool>>,
}

impl CooksReplacementOptions {
    /// Construct replacement options with DESeq2's default trim and replicate count.
    pub fn new(cooks_cutoff: f64) -> Self {
        Self {
            trim: 0.2,
            cooks_cutoff,
            min_replicates: 7,
            which_samples: None,
        }
    }
}

/// Output from DESeq2-style Cook's outlier count replacement.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksReplacementOutput {
    /// Counts after applying replacements in replaceable outlier cells.
    pub replaced_counts: CountMatrix,
    /// Candidate replacement counts for every gene/sample cell.
    pub candidate_replacement_counts: CountMatrix,
    /// Per-cell finite Cook's outlier flags, `cooks > cutoff`.
    pub outlier_cells: RowMajorMatrix<bool>,
    /// Per-gene `replace` flags matching R's `any(cooks > cutoff)` NA shape.
    pub replace: Vec<Option<bool>>,
    /// Samples eligible for replacement.
    pub replaceable_samples: Vec<bool>,
}

/// Refit-planning metadata for DESeq2-style Cook's outlier replacement.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksRefitPlan {
    /// Primitive replacement output from `replaceOutliers`.
    pub replacement: CooksReplacementOutput,
    /// Normalized counts after applying replaceable Cook's outlier replacements.
    pub replaced_normalized_counts: RowMajorMatrix<f64>,
    /// `baseMean` recomputed on replacement counts.
    pub replaced_base_mean: Vec<f64>,
    /// `baseVar` recomputed on replacement counts.
    pub replaced_base_var: Vec<f64>,
    /// All-zero flags recomputed on replacement counts.
    pub replaced_all_zero: Vec<bool>,
    /// Number of genes with `replace == TRUE`, matching `sum(replace, na.rm=TRUE)`.
    pub n_refit: usize,
    /// Candidate genes for refitting, before the workflow eligibility guard.
    pub refit_rows: Vec<usize>,
    /// Replaced genes that became all-zero after replacement.
    pub new_all_zero_rows: Vec<usize>,
    /// Whether candidate rows need refitting, before sample eligibility.
    /// Use [`Self::requires_refit`] to choose workflow execution.
    pub should_refit: bool,
    /// `maxCooks` after replacement refit, using original Cook's distances with
    /// replaceable sample columns ignored by zeroing, matching DESeq2.
    pub post_refit_max_cooks: Vec<Option<f64>>,
}

/// Compact metadata summary for a Cook's replacement/refit plan.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CooksReplacementMetadata {
    /// Number of genes with `replace == TRUE`, matching `sum(replace, na.rm=TRUE)`.
    pub n_refit: usize,
    /// Number of replacement-marked candidate genes for GLM refitting.
    pub n_refit_rows: usize,
    /// Number of replacement-marked genes that became all-zero.
    pub n_new_all_zero: usize,
    /// Number of finite Cook's outlier cells before sample eligibility filtering.
    pub n_outlier_cells: usize,
    /// Number of Cook's outlier cells in replaceable samples.
    pub n_replaced_cells: usize,
    /// Number of samples eligible for Cook's outlier replacement.
    pub n_replaceable_samples: usize,
    /// Candidate refit condition before sample eligibility is applied.
    pub should_refit: bool,
}

/// Name/value metadata entry for Cook's replacement/refit summaries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CooksReplacementMetadataEntry {
    /// Entry name.
    pub name: String,
    /// Entry value, formatted as text.
    pub value: String,
}

impl CooksReplacementMetadata {
    /// Return stable scalar metadata entries suitable for TSV/object metadata export.
    pub fn scalar_metadata(&self) -> Vec<CooksReplacementMetadataEntry> {
        vec![
            CooksReplacementMetadataEntry::new("nRefit", self.n_refit),
            CooksReplacementMetadataEntry::new("nRefitRows", self.n_refit_rows),
            CooksReplacementMetadataEntry::new("nNewAllZero", self.n_new_all_zero),
            CooksReplacementMetadataEntry::new("nOutlierCells", self.n_outlier_cells),
            CooksReplacementMetadataEntry::new("nReplacedCells", self.n_replaced_cells),
            CooksReplacementMetadataEntry::new("nReplaceableSamples", self.n_replaceable_samples),
            CooksReplacementMetadataEntry::new("shouldRefit", self.should_refit),
        ]
    }
}

impl CooksReplacementMetadataEntry {
    fn new(name: impl Into<String>, value: impl ToString) -> Self {
        Self {
            name: name.into(),
            value: value.to_string(),
        }
    }
}

impl CooksRefitPlan {
    /// Whether a public fitting workflow has any samples eligible for replacement.
    pub fn has_replaceable_samples(&self) -> bool {
        self.replacement
            .replaceable_samples
            .iter()
            .any(|&value| value)
    }

    /// Whether a fitting workflow should execute the candidate refit.
    /// The raw planner's `should_refit` and row metadata remain available even
    /// when no sample is eligible; those candidates must not trigger fitting.
    pub fn requires_refit(&self) -> bool {
        self.should_refit && self.has_replaceable_samples()
    }

    /// Summarize the replacement/refit branch in a compact, object-friendly form.
    pub fn metadata(&self) -> CooksReplacementMetadata {
        let n_samples = self.replacement.outlier_cells.n_cols();
        let mut n_outlier_cells = 0;
        let mut n_replaced_cells = 0;
        for (idx, is_outlier) in self
            .replacement
            .outlier_cells
            .as_slice()
            .iter()
            .copied()
            .enumerate()
        {
            if !is_outlier {
                continue;
            }
            n_outlier_cells += 1;
            let sample = idx % n_samples;
            if self.replacement.replaceable_samples[sample] {
                n_replaced_cells += 1;
            }
        }
        CooksReplacementMetadata {
            n_refit: self.n_refit,
            n_refit_rows: self.refit_rows.len(),
            n_new_all_zero: self.new_all_zero_rows.len(),
            n_outlier_cells,
            n_replaced_cells,
            n_replaceable_samples: self
                .replacement
                .replaceable_samples
                .iter()
                .filter(|value| **value)
                .count(),
            should_refit: self.should_refit,
        }
    }

    /// Return stable scalar replacement/refit metadata entries.
    pub fn scalar_metadata(&self) -> Vec<CooksReplacementMetadataEntry> {
        self.metadata().scalar_metadata()
    }
}

/// Calculate DESeq2-style Cook's distances.
///
/// This mirrors the formula in DESeq2's `calculateCooksDistance`: Pearson
/// residuals are scaled by the hat diagonal and divided by the number of model
/// coefficients. The variance uses DESeq2's robust method-of-moments
/// dispersion estimate rather than the final fitted dispersion.
pub fn calculate_cooks_distance(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    mu: &RowMajorMatrix<f64>,
    hat_diagonal: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
) -> Result<CooksOutput, DeseqError> {
    validate_cooks_inputs(counts, normalized_counts, mu, hat_diagonal, model_matrix)?;
    let groups = model_matrix_groups(model_matrix)?;
    let robust_dispersion = robust_dispersion_for_groups(normalized_counts, &groups)?;
    let p = model_matrix.n_coefficients() as f64;

    // Genes are independent: fill each gene's row in parallel and report the
    // lowest failing gene's error, as the sequential loop did.
    let n_samples = counts.n_samples();
    let mut values = try_alloc_matrix(counts.n_genes(), n_samples, "Cook's distances")?;
    values.resize(counts.n_genes() * n_samples, f64::NAN);
    let first_error = values
        .par_chunks_mut(n_samples.max(1))
        .zip(robust_dispersion.par_iter().copied())
        .enumerate()
        .map(|(gene, (out, dispersion))| -> Result<(), DeseqError> {
            let count_row = counts.row(gene)?;
            let mu_row = mu.row(gene)?;
            let h_row = hat_diagonal.row(gene)?;
            for sample in 0..counts.n_samples() {
                if dispersion.is_nan() {
                    out[sample] = f64::NAN;
                    continue;
                }
                let mean = mu_row[sample];
                let h = h_row[sample];
                if mean.is_nan() || h.is_nan() {
                    out[sample] = f64::NAN;
                    continue;
                }
                let Ok(dispersion_mean) = checked::try_mul(dispersion, mean) else {
                    out[sample] = f64::NAN;
                    continue;
                };
                let Ok(variance_factor) = checked::try_add(1.0, dispersion_mean) else {
                    out[sample] = f64::NAN;
                    continue;
                };
                let Ok(variance) = checked::try_mul(mean, variance_factor) else {
                    out[sample] = f64::NAN;
                    continue;
                };
                let residual = f64::from(count_row[sample]) - mean;
                let one_minus_h = 1.0 - h;
                let Ok(residual_sq) = checked::try_mul(residual, residual) else {
                    out[sample] = f64::NAN;
                    continue;
                };
                let pearson_sq = residual_sq / variance;
                let Ok(denominator) = checked::try_mul(one_minus_h, one_minus_h) else {
                    out[sample] = f64::NAN;
                    continue;
                };
                let value = pearson_sq / p * h / denominator;
                if !value.is_finite() {
                    out[sample] = f64::NAN;
                    continue;
                }
                out[sample] = value;
            }
            Ok(())
        })
        .find_map_first(Result::err);
    if let Some(error) = first_error {
        return Err(error);
    }

    let cooks = RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)?;
    let samples_for_cooks = eligible_samples(&groups, model_matrix.n_samples(), 3)?;
    let max_cooks = record_max_cooks(&cooks, model_matrix, &samples_for_cooks)?;
    Ok(CooksOutput {
        cooks,
        max_cooks,
        robust_dispersion,
        samples_for_cooks,
    })
}

// Replacement and refit-plan construction.
mod replacement;
pub use replacement::{
    prepare_cooks_replacement_refit, prepare_cooks_replacement_refit_with_weights,
    replace_outlier_counts,
};

/// Calculate DESeq2-style post-refit `maxCooks` after outlier replacement.
///
/// DESeq2 preserves original Cook's distances for diagnostics, but after a
/// replacement-triggered refit it excludes replaceable sample columns from
/// `maxCooks` by setting those Cook's distances to zero before recording the
/// maximum. When all samples are replaceable, the post-refit `maxCooks` field
/// is missing for every gene.
pub fn max_cooks_after_replacement_refit(
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    replaceable_samples: &[bool],
) -> Result<Vec<Option<f64>>, DeseqError> {
    if cooks.n_cols() != model_matrix.n_samples() {
        return Err(invalid_dimensions(
            "replacement-refit Cook's samples",
            model_matrix.n_samples(),
            cooks.n_cols(),
        ));
    }
    if replaceable_samples.len() != model_matrix.n_samples() {
        return Err(invalid_dimensions(
            "replacement-refit replaceable samples",
            model_matrix.n_samples(),
            replaceable_samples.len(),
        ));
    }
    if replaceable_samples.iter().all(|value| *value) {
        return try_filled_vec(cooks.n_rows(), None, "post-refit maxCooks");
    }

    let samples = samples_for_cooks(model_matrix, 3)?;
    record_max_cooks_masked(cooks, model_matrix, &samples, Some(replaceable_samples))
}

/// Robust method-of-moments dispersion used by DESeq2 for Cook's distances.
pub fn robust_method_of_moments_dispersion(
    normalized_counts: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
) -> Result<Vec<f64>, DeseqError> {
    if normalized_counts.n_cols() != model_matrix.n_samples() {
        return Err(invalid_dimensions(
            "Cook's robust dispersion samples",
            model_matrix.n_samples(),
            normalized_counts.n_cols(),
        ));
    }

    let groups = model_matrix_groups(model_matrix)?;
    robust_dispersion_for_groups(normalized_counts, &groups)
}

fn robust_dispersion_for_groups(
    normalized_counts: &RowMajorMatrix<f64>,
    groups: &[Vec<usize>],
) -> Result<Vec<f64>, DeseqError> {
    let mut dispersions = trimmed_variances(normalized_counts, groups)?;
    for (gene, variance) in dispersions.iter_mut().enumerate() {
        let row = normalized_counts.row(gene)?;
        let mean = checked::try_mean(row).nonfinite_at("Cook's robust dispersion mean", None)?;
        let alpha = if mean > 0.0 {
            let inv_mean = mean.recip();
            let centered = *variance - mean;
            if !centered.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "Cook's robust dispersion centered variance".to_string(),
                    index: Some(gene),
                    value: centered,
                });
            }
            checked::try_mul(
                centered,
                checked::try_mul(inv_mean, inv_mean)
                    .nonfinite_at("Cook's robust dispersion inverse mean square", Some(gene))?,
            )
            .nonfinite_at("Cook's robust dispersion", Some(gene))?
        } else {
            f64::NAN
        };
        *variance = if alpha.is_nan() {
            f64::NAN
        } else {
            alpha.max(MIN_COOKS_DISPERSION)
        };
    }
    Ok(dispersions)
}

/// Return samples belonging to model-matrix cells with at least `n` replicates.
pub fn samples_for_cooks(model_matrix: &DesignMatrix, n: usize) -> Result<Vec<bool>, DeseqError> {
    let groups = model_matrix_groups(model_matrix)?;
    eligible_samples(&groups, model_matrix.n_samples(), n)
}

fn eligible_samples(
    groups: &[Vec<usize>],
    samples: usize,
    minimum: usize,
) -> Result<Vec<bool>, DeseqError> {
    let mut eligible = try_filled_vec(samples, false, "Cook's eligible samples")?;
    for group in groups.iter().filter(|group| group.len() >= minimum) {
        for &sample in group {
            eligible[sample] = true;
        }
    }
    Ok(eligible)
}

/// Record maximum Cook's distance over eligible samples.
pub fn record_max_cooks(
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    samples_for_cooks: &[bool],
) -> Result<Vec<Option<f64>>, DeseqError> {
    if cooks.n_cols() != model_matrix.n_samples() {
        return Err(invalid_dimensions(
            "maxCooks samples",
            model_matrix.n_samples(),
            cooks.n_cols(),
        ));
    }
    if samples_for_cooks.len() != model_matrix.n_samples() {
        return Err(invalid_dimensions(
            "samplesForCooks",
            model_matrix.n_samples(),
            samples_for_cooks.len(),
        ));
    }
    record_max_cooks_masked(cooks, model_matrix, samples_for_cooks, None)
}

fn record_max_cooks_masked(
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    samples_for_cooks: &[bool],
    zero_samples: Option<&[bool]>,
) -> Result<Vec<Option<f64>>, DeseqError> {
    if model_matrix.n_samples() <= model_matrix.n_coefficients()
        || !samples_for_cooks.iter().any(|value| *value)
    {
        return try_filled_vec(cooks.n_rows(), None, "maxCooks");
    }

    // Row maxima are independent; rows are contiguous in the row-major matrix.
    let mut max_values = try_alloc_matrix(cooks.n_rows(), 1, "maxCooks")?;
    cooks
        .as_slice()
        .par_chunks(cooks.n_cols())
        .map(|row| {
            let mut row_max: Option<f64> = None;
            for (sample, use_sample) in samples_for_cooks.iter().copied().enumerate() {
                if !use_sample {
                    continue;
                }
                let value = if zero_samples.is_some_and(|mask| mask[sample]) {
                    0.0
                } else {
                    row[sample]
                };
                if value.is_nan() {
                    row_max = None;
                    break;
                }
                row_max = Some(row_max.map_or(value, |current| current.max(value)));
            }
            row_max
        })
        .collect_into_vec(&mut max_values);
    Ok(max_values)
}

fn validate_cooks_inputs(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    mu: &RowMajorMatrix<f64>,
    hat_diagonal: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
) -> Result<(), DeseqError> {
    validate_gene_sample_matrix(
        "normalized counts",
        normalized_counts,
        counts.n_genes(),
        counts.n_samples(),
    )?;
    validate_gene_sample_matrix("Cook's mu", mu, counts.n_genes(), counts.n_samples())?;
    validate_gene_sample_matrix(
        "Cook's hat diagonals",
        hat_diagonal,
        counts.n_genes(),
        counts.n_samples(),
    )?;
    if model_matrix.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "Cook's model matrix samples",
            counts.n_samples(),
            model_matrix.n_samples(),
        ));
    }
    Ok(())
}

fn validate_replacement_inputs(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    options: &CooksReplacementOptions,
) -> Result<(), DeseqError> {
    validate_gene_sample_matrix(
        "replacement normalized counts",
        normalized_counts,
        counts.n_genes(),
        counts.n_samples(),
    )?;
    normalized_counts.validate_finite("replacement normalized counts")?;
    validate_gene_sample_matrix(
        "replacement cooks",
        cooks,
        counts.n_genes(),
        counts.n_samples(),
    )?;
    if model_matrix.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "replacement model matrix samples",
            counts.n_samples(),
            model_matrix.n_samples(),
        ));
    }
    if !options.cooks_cutoff.is_finite() || options.cooks_cutoff < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "replacement Cook's cutoff".to_string(),
            index: None,
            value: options.cooks_cutoff,
        });
    }
    if !options.trim.is_finite() || !(0.0..=0.5).contains(&options.trim) {
        return Err(DeseqError::InvalidOptions {
            reason: "replacement trim must be finite and in [0, 0.5]".to_string(),
        });
    }
    if options.min_replicates < 3 {
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "outlier replacement",
            "at least 3 replicates",
        ));
    }
    if let Some(samples) = &options.which_samples
        && samples.len() != counts.n_samples()
    {
        return Err(invalid_dimensions(
            "replacement whichSamples",
            counts.n_samples(),
            samples.len(),
        ));
    }
    Ok(())
}

fn validate_gene_sample_matrix(
    context: &str,
    matrix: &RowMajorMatrix<f64>,
    n_genes: usize,
    n_samples: usize,
) -> Result<(), DeseqError> {
    if matrix.n_rows() != n_genes {
        return Err(invalid_dimensions(
            format!("{context} rows"),
            n_genes,
            matrix.n_rows(),
        ));
    }
    if matrix.n_cols() != n_samples {
        return Err(invalid_dimensions(
            format!("{context} columns"),
            n_samples,
            matrix.n_cols(),
        ));
    }
    Ok(())
}

// Robust dispersion statistics.
mod trimmed_variance;
use trimmed_variance::trimmed_variances;

fn model_matrix_groups(model_matrix: &DesignMatrix) -> Result<Vec<Vec<usize>>, DeseqError> {
    let samples = model_matrix.n_samples();
    let columns = model_matrix.n_coefficients();
    let values = model_matrix.matrix().as_slice();
    let mut groups: Vec<Vec<usize>> = try_alloc_matrix(samples, 1, "Cook's model matrix cells")?;
    for (sample, row) in values.chunks_exact(columns).enumerate() {
        if let Some(group) = groups.iter_mut().find(|group| {
            let start = group[0] * columns;
            &values[start..start + columns] == row
        }) {
            group.push(sample);
        } else {
            groups.push(vec![sample]);
        }
    }
    Ok(groups)
}

fn cooks_row_buffer(
    scratch: &mut Option<Vec<f64>>,
    samples: usize,
) -> Result<&mut Vec<f64>, DeseqError> {
    match scratch {
        Some(values) => Ok(values),
        None => Ok(scratch.insert(try_alloc_matrix(samples, 1, "Cook's row scratch")?)),
    }
}

/// R `mean(x, trim)`: median for `trim >= 0.5`, otherwise a plain sum of the
/// central values divided by their count; NaN for non-finite input.
fn r_trimmed_mean(values: &mut [f64], trim: f64) -> f64 {
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return f64::NAN;
    }
    values.sort_unstable_by(f64::total_cmp);
    if trim >= 0.5 {
        return if values.len().is_multiple_of(2) {
            let upper = values.len() / 2;
            0.5 * (values[upper - 1] + values[upper])
        } else {
            values[values.len() / 2]
        };
    }
    let trim_count = (values.len() as f64 * trim).floor() as usize;
    let end = values.len().saturating_sub(trim_count);
    if trim_count >= end {
        return f64::NAN;
    }
    checked::try_sum(values[trim_count..end].iter().copied())
        .nonfinite("replacement trimmed mean")
        .map(|sum| sum / (end - trim_count) as f64)
        .unwrap_or(f64::NAN)
}

fn replacement_count_from_scaled_mean(value: f64) -> Result<u32, DeseqError> {
    if !value.is_finite() || value < 0.0 || value > f64::from(u32::MAX) {
        return Err(DeseqError::InvalidCounts {
            reason: format!("replacement count must be finite and fit in u32, got {value}"),
        });
    }
    Ok(value.trunc() as u32)
}

fn count_matrix_like(counts: &CountMatrix, values: Vec<u32>) -> Result<CountMatrix, DeseqError> {
    CountMatrix::from_row_major_u32_with_names(
        counts.n_genes(),
        counts.n_samples(),
        values,
        counts.gene_names().map(<[String]>::to_vec),
        counts.sample_names().map(<[String]>::to_vec),
    )
}

fn no_replacement_output(counts: &CountMatrix) -> Result<CooksReplacementOutput, DeseqError> {
    Ok(CooksReplacementOutput {
        replaced_counts: count_matrix_like(counts, counts.as_slice().to_vec())?,
        candidate_replacement_counts: count_matrix_like(counts, counts.as_slice().to_vec())?,
        outlier_cells: RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), {
            let mut values =
                try_alloc_matrix(counts.n_genes(), counts.n_samples(), "Cook's outlier mask")?;
            values.resize(counts.n_genes() * counts.n_samples(), false);
            values
        })?,
        replace: vec![Some(false); counts.n_genes()],
        replaceable_samples: vec![false; counts.n_samples()],
    })
}
