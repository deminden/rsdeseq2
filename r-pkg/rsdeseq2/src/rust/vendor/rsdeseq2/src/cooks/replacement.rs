use super::*;

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
    match normalization_factors {
        Some(factors) => {
            validate_normalization_factors(counts, factors)?;
        }
        None => validate_size_factors(size_factors, counts.n_samples())?,
    }
    let replaceable_samples = match &options.which_samples {
        Some(samples) => samples.clone(),
        None => samples_for_cooks(model_matrix, options.min_replicates)?,
    };

    let mut candidate_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let trim_mean = r_trimmed_mean(normalized_counts.row(gene)?.to_vec(), options.trim);
        let scale_row = match normalization_factors {
            Some(factors) => factors.row(gene)?,
            None => size_factors,
        };
        for (sample, scale) in scale_row.iter().copied().enumerate() {
            let scaled_mean = checked_product2(
                trim_mean,
                scale,
                Some(gene * counts.n_samples() + sample),
                "replacement scaled mean",
            )?;
            candidate_values.push(replacement_count_from_scaled_mean(scaled_mean)?);
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
    let replaced_normalized_counts = match normalization_factors {
        Some(factors) => {
            validate_normalization_factors(counts, factors)?;
            normalized_counts_with_factors(&replacement.replaced_counts, factors)?
        }
        None => {
            crate::normalization::normalized_counts(&replacement.replaced_counts, size_factors)?
        }
    };
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
