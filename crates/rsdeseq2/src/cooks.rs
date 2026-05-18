use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{
    base_mean, base_variance, normalization_factors_from_size_factors,
    normalized_counts_with_factors, validate_normalization_factors,
};

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
    /// Genes that should be refit after replacement.
    pub refit_rows: Vec<usize>,
    /// Replaced genes that became all-zero after replacement.
    pub new_all_zero_rows: Vec<usize>,
    /// Whether DESeq2 would enter the replacement-refit branch.
    pub should_refit: bool,
    /// `maxCooks` after replacement refit, using original Cook's distances with
    /// replaceable sample columns ignored by zeroing, matching DESeq2.
    pub post_refit_max_cooks: Vec<Option<f64>>,
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
    let robust_dispersion = robust_method_of_moments_dispersion(normalized_counts, model_matrix)?;
    let p = model_matrix.n_coefficients() as f64;

    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for (gene, dispersion) in robust_dispersion.iter().copied().enumerate() {
        let count_row = counts.row(gene)?;
        let mu_row = mu.row(gene)?;
        let h_row = hat_diagonal.row(gene)?;
        for sample in 0..counts.n_samples() {
            let mean = mu_row[sample];
            let h = h_row[sample];
            let variance = mean + dispersion * mean.powi(2);
            let pearson_sq = (f64::from(count_row[sample]) - mean).powi(2) / variance;
            values.push(pearson_sq / p * h / (1.0 - h).powi(2));
        }
    }

    let cooks = RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)?;
    let samples_for_cooks = samples_for_cooks(model_matrix, 3)?;
    let max_cooks = record_max_cooks(&cooks, model_matrix, &samples_for_cooks)?;
    Ok(CooksOutput {
        cooks,
        max_cooks,
        robust_dispersion,
        samples_for_cooks,
    })
}

/// Replace Cook's outlier counts with trimmed-mean predictions.
///
/// This mirrors the primitive count-transformation part of DESeq2
/// `replaceOutliers`: compute a trimmed mean over normalized counts for each
/// gene, rescale by size factors or gene/sample normalization factors, coerce to
/// integer counts by truncation, and replace only cells with `cooks > cutoff` in
/// samples eligible by `whichSamples` or `minReplicates`.
///
/// This function does not perform the later DESeq2 refit cycle. It returns the
/// transformed counts and metadata needed by a future refit stage.
pub fn replace_outlier_counts(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    size_factors: &[f64],
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    options: &CooksReplacementOptions,
) -> Result<CooksReplacementOutput, DeseqError> {
    validate_replacement_inputs(counts, normalized_counts, cooks, model_matrix, options)?;
    if counts.n_samples() <= model_matrix.n_coefficients() {
        return no_replacement_output(counts);
    }
    let scale_factors = match normalization_factors {
        Some(factors) => {
            validate_normalization_factors(counts, factors)?;
            factors.clone()
        }
        None => normalization_factors_from_size_factors(counts, size_factors)?,
    };
    let replaceable_samples = match &options.which_samples {
        Some(samples) => samples.clone(),
        None => samples_for_cooks(model_matrix, options.min_replicates)?,
    };

    let mut candidate_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let trim_mean = r_trimmed_mean(normalized_counts.row(gene)?.to_vec(), options.trim);
        let scale_row = scale_factors.row(gene)?;
        for scale in scale_row {
            candidate_values.push(replacement_count_from_scaled_mean(trim_mean * *scale)?);
        }
    }

    let mut outlier_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut replace = Vec::with_capacity(counts.n_genes());
    let mut replaced_values = counts.as_slice().to_vec();
    for gene in 0..counts.n_genes() {
        let cooks_row = cooks.row(gene)?;
        let mut row_has_outlier = false;
        let mut row_has_missing = false;
        for sample in 0..counts.n_samples() {
            let cook = cooks_row[sample];
            let is_outlier = cook.is_finite() && cook > options.cooks_cutoff;
            if !cook.is_finite() {
                row_has_missing = true;
            }
            if is_outlier {
                row_has_outlier = true;
                if replaceable_samples[sample] {
                    let idx = gene * counts.n_samples() + sample;
                    replaced_values[idx] = candidate_values[idx];
                }
            }
            outlier_values.push(is_outlier);
        }
        replace.push(if row_has_outlier {
            Some(true)
        } else if row_has_missing {
            None
        } else {
            Some(false)
        });
    }

    Ok(CooksReplacementOutput {
        replaced_counts: count_matrix_like(counts, replaced_values)?,
        candidate_replacement_counts: count_matrix_like(counts, candidate_values)?,
        outlier_cells: RowMajorMatrix::from_row_major(
            counts.n_genes(),
            counts.n_samples(),
            outlier_values,
        )?,
        replace,
        replaceable_samples,
    })
}

/// Prepare DESeq2-style Cook's outlier replacement metadata for a later refit.
///
/// This mirrors the bookkeeping in DESeq2 `refitWithoutOutliers` after
/// `replaceOutliers`: recompute base metadata on replacement counts, identify
/// rows marked for refitting, separate replacement rows that became all-zero,
/// and calculate the post-refit `maxCooks` masking rule. It deliberately does
/// not estimate dispersions or refit the GLM.
pub fn prepare_cooks_replacement_refit(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    size_factors: &[f64],
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    options: &CooksReplacementOptions,
) -> Result<CooksRefitPlan, DeseqError> {
    let replacement = replace_outlier_counts(
        counts,
        normalized_counts,
        size_factors,
        normalization_factors,
        cooks,
        model_matrix,
        options,
    )?;
    let scale_factors = match normalization_factors {
        Some(factors) => {
            validate_normalization_factors(counts, factors)?;
            factors.clone()
        }
        None => normalization_factors_from_size_factors(counts, size_factors)?,
    };
    let replaced_normalized_counts =
        normalized_counts_with_factors(&replacement.replaced_counts, &scale_factors)?;
    let replaced_base_mean = base_mean(&replaced_normalized_counts)?;
    let replaced_base_var = base_variance(&replaced_normalized_counts)?;
    let replaced_all_zero = replacement.replaced_counts.all_zero_flags();

    let mut n_refit = 0;
    let mut refit_rows = Vec::new();
    let mut new_all_zero_rows = Vec::new();
    for (gene, replace) in replacement.replace.iter().copied().enumerate() {
        if replace != Some(true) {
            continue;
        }
        n_refit += 1;
        if replaced_all_zero[gene] {
            new_all_zero_rows.push(gene);
        } else {
            refit_rows.push(gene);
        }
    }
    let post_refit_max_cooks =
        max_cooks_after_replacement_refit(cooks, model_matrix, &replacement.replaceable_samples)?;

    Ok(CooksRefitPlan {
        replacement,
        replaced_normalized_counts,
        replaced_base_mean,
        replaced_base_var,
        replaced_all_zero,
        n_refit,
        should_refit: n_refit > new_all_zero_rows.len(),
        refit_rows,
        new_all_zero_rows,
        post_refit_max_cooks,
    })
}

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
        return Ok(vec![None; cooks.n_rows()]);
    }

    let mut values = cooks.as_slice().to_vec();
    for gene in 0..cooks.n_rows() {
        let row_start = gene * cooks.n_cols();
        for (sample, replaceable) in replaceable_samples.iter().copied().enumerate() {
            if replaceable {
                values[row_start + sample] = 0.0;
            }
        }
    }
    let replace_cooks = RowMajorMatrix::from_row_major(cooks.n_rows(), cooks.n_cols(), values)?;
    let samples = samples_for_cooks(model_matrix, 3)?;
    record_max_cooks(&replace_cooks, model_matrix, &samples)
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

    let samples_with_three_or_more = samples_for_cooks(model_matrix, 3)?;
    let variance = if samples_with_three_or_more.iter().any(|value| *value) {
        trimmed_cell_variance(normalized_counts, model_matrix)?
    } else {
        trimmed_variance(normalized_counts)?
    };

    let mut dispersions = Vec::with_capacity(normalized_counts.n_rows());
    for (gene, variance) in variance.iter().copied().enumerate() {
        let row = normalized_counts.row(gene)?;
        let mean = row.iter().sum::<f64>() / row.len() as f64;
        let alpha = if mean > 0.0 {
            (variance - mean) / mean.powi(2)
        } else {
            f64::NAN
        };
        dispersions.push(if alpha.is_nan() {
            f64::NAN
        } else {
            alpha.max(MIN_COOKS_DISPERSION)
        });
    }
    Ok(dispersions)
}

/// Return samples belonging to model-matrix cells with at least `n` replicates.
pub fn samples_for_cooks(model_matrix: &DesignMatrix, n: usize) -> Result<Vec<bool>, DeseqError> {
    let groups = model_matrix_groups(model_matrix)?;
    Ok(groups
        .iter()
        .map(|group| groups.iter().filter(|other| *other == group).count() >= n)
        .collect())
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
    if model_matrix.n_samples() <= model_matrix.n_coefficients()
        || !samples_for_cooks.iter().any(|value| *value)
    {
        return Ok(vec![None; cooks.n_rows()]);
    }

    let mut max_values = Vec::with_capacity(cooks.n_rows());
    for gene in 0..cooks.n_rows() {
        let row = cooks.row(gene)?;
        let mut row_max: Option<f64> = None;
        for (sample, use_sample) in samples_for_cooks.iter().copied().enumerate() {
            if !use_sample {
                continue;
            }
            let value = row[sample];
            if value.is_nan() {
                row_max = None;
                break;
            }
            row_max = Some(row_max.map_or(value, |current| current.max(value)));
        }
        max_values.push(row_max);
    }
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
        return Err(DeseqError::InvalidOptions {
            reason: "at least 3 replicates are required for outlier replacement".to_string(),
        });
    }
    if let Some(samples) = &options.which_samples {
        if samples.len() != counts.n_samples() {
            return Err(invalid_dimensions(
                "replacement whichSamples",
                counts.n_samples(),
                samples.len(),
            ));
        }
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

fn trimmed_cell_variance(
    normalized_counts: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
) -> Result<Vec<f64>, DeseqError> {
    let groups = model_matrix_groups(model_matrix)?;
    let mut group_ids = unique_group_ids(&groups);
    group_ids.retain(|group| {
        groups
            .iter()
            .filter(|candidate| *candidate == group)
            .count()
            >= 3
    });

    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let mut group_variances = Vec::with_capacity(group_ids.len());
        for group in &group_ids {
            let values = groups
                .iter()
                .enumerate()
                .filter_map(|(sample, candidate)| (candidate == group).then_some(row[sample]))
                .collect::<Vec<_>>();
            let bin = trim_bin(values.len());
            let mean = trimmed_mean(values.clone(), trim_ratio(bin));
            let sq_errors = values
                .into_iter()
                .map(|value| (value - mean).powi(2))
                .collect::<Vec<_>>();
            group_variances.push(trim_scale(bin) * trimmed_mean(sq_errors, trim_ratio(bin)));
        }
        variances.push(row_max(&group_variances));
    }
    Ok(variances)
}

fn trimmed_variance(normalized_counts: &RowMajorMatrix<f64>) -> Result<Vec<f64>, DeseqError> {
    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let mean = trimmed_mean(row.to_vec(), 1.0 / 8.0);
        let sq_errors = row
            .iter()
            .copied()
            .map(|value| (value - mean).powi(2))
            .collect::<Vec<_>>();
        variances.push(1.51 * trimmed_mean(sq_errors, 1.0 / 8.0));
    }
    Ok(variances)
}

fn model_matrix_groups(model_matrix: &DesignMatrix) -> Result<Vec<Vec<f64>>, DeseqError> {
    let mut groups = Vec::with_capacity(model_matrix.n_samples());
    for sample in 0..model_matrix.n_samples() {
        groups.push(model_matrix.matrix().row(sample)?.to_vec());
    }
    Ok(groups)
}

fn unique_group_ids(groups: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut unique = Vec::new();
    for group in groups {
        if !unique.iter().any(|candidate| candidate == group) {
            unique.push(group.clone());
        }
    }
    unique
}

fn trim_bin(n: usize) -> usize {
    if n <= 3 {
        0
    } else if n <= 23 {
        1
    } else {
        2
    }
}

fn trim_ratio(bin: usize) -> f64 {
    [1.0 / 3.0, 1.0 / 4.0, 1.0 / 8.0][bin]
}

fn trim_scale(bin: usize) -> f64 {
    [2.04, 1.86, 1.51][bin]
}

fn trimmed_mean(mut values: Vec<f64>, trim: f64) -> f64 {
    if values.is_empty() || values.iter().any(|value| value.is_nan()) {
        return f64::NAN;
    }
    values.sort_by(f64::total_cmp);
    let trim_count = (values.len() as f64 * trim).floor() as usize;
    let end = values.len().saturating_sub(trim_count);
    if trim_count >= end {
        return f64::NAN;
    }
    let kept = &values[trim_count..end];
    kept.iter().sum::<f64>() / kept.len() as f64
}

fn r_trimmed_mean(mut values: Vec<f64>, trim: f64) -> f64 {
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return f64::NAN;
    }
    values.sort_by(f64::total_cmp);
    if trim >= 0.5 {
        return if values.len() % 2 == 0 {
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
    values[trim_count..end].iter().sum::<f64>() / (end - trim_count) as f64
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
        outlier_cells: RowMajorMatrix::from_row_major(
            counts.n_genes(),
            counts.n_samples(),
            vec![false; counts.n_genes() * counts.n_samples()],
        )?,
        replace: vec![Some(false); counts.n_genes()],
        replaceable_samples: vec![false; counts.n_samples()],
    })
}

fn row_max(values: &[f64]) -> f64 {
    let mut max_value = f64::NEG_INFINITY;
    for value in values {
        if value.is_nan() {
            return f64::NAN;
        }
        max_value = max_value.max(*value);
    }
    max_value
}
