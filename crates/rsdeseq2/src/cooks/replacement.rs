use super::*;
use crate::matrix::try_alloc_matrix;
use crate::numeric::checked::{self, CheckedExt};

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

    let genes = counts.n_genes();
    let samples = counts.n_samples();
    let mut candidate_values = try_alloc_matrix(genes, samples, "Cook's replacement")?;
    let mut outlier_values = try_alloc_matrix(genes, samples, "Cook's replacement")?;
    let mut replaced_values = try_alloc_matrix(genes, samples, "Cook's replacement")?;
    candidate_values.resize(counts.as_slice().len(), 0);
    outlier_values.resize(counts.as_slice().len(), false);
    replaced_values.extend_from_slice(counts.as_slice());

    let replace = candidate_values
        .par_chunks_mut(samples)
        .zip(outlier_values.par_chunks_mut(samples))
        .zip(replaced_values.par_chunks_mut(samples))
        .enumerate()
        .map(|(gene, ((candidate, outlier), replaced))| {
            let trim_mean = r_trimmed_mean(normalized_counts.row(gene)?.to_vec(), options.trim);
            let scales = match normalization_factors {
                Some(factors) => factors.row(gene)?,
                None => size_factors,
            };
            let cooks_row = cooks.row(gene)?;
            let mut has_outlier = false;
            let mut has_missing = false;
            for sample in 0..samples {
                let scaled_mean = checked::try_mul(trim_mean, scales[sample])
                    .nonfinite_at("replacement scaled mean", Some(gene * samples + sample))?;
                candidate[sample] = replacement_count_from_scaled_mean(scaled_mean)?;

                let cook = cooks_row[sample];
                let finite = cook.is_finite();
                outlier[sample] = finite && cook > options.cooks_cutoff;
                has_outlier |= outlier[sample];
                has_missing |= !finite;
                if outlier[sample] && replaceable_samples[sample] {
                    replaced[sample] = candidate[sample];
                }
            }
            Ok(if has_outlier {
                Some(true)
            } else if has_missing {
                None
            } else {
                Some(false)
            })
        })
        .collect::<Vec<_>>();
    // Collect in row order so concurrent failures report the same first gene
    // and sample as the serial implementation.
    let replace = replace
        .into_iter()
        .collect::<Result<Vec<_>, DeseqError>>()?;

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
    prepare_cooks_replacement_refit_with_weights(
        counts,
        normalized_counts,
        size_factors,
        normalization_factors,
        cooks,
        model_matrix,
        options,
        None,
    )
}

/// Prepare replacement/refit metadata using the original observation weights.
/// Weights affect base means/variances, not trimmed replacement counts. Pass raw
/// assay weights, before row normalization used by the GLM score.
#[allow(clippy::too_many_arguments)]
pub fn prepare_cooks_replacement_refit_with_weights(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    size_factors: &[f64],
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    options: &CooksReplacementOptions,
    observation_weights: Option<&RowMajorMatrix<f64>>,
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
    let (replaced_base_mean, replaced_base_var) = match observation_weights {
        Some(weights) => (
            crate::normalization::base_mean_with_weights(&replaced_normalized_counts, weights)?,
            crate::normalization::base_variance_with_weights(&replaced_normalized_counts, weights)?,
        ),
        None => (
            base_mean(&replaced_normalized_counts)?,
            base_variance(&replaced_normalized_counts)?,
        ),
    };
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
