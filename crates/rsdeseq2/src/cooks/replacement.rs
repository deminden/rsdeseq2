use super::*;
use crate::matrix::try_alloc_matrix;
use crate::normalization::{Offsets, normalized_counts_for_rows, renormalize_changed_rows};
use crate::numeric::checked::{self, CheckedExt};
use std::borrow::Cow;

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
    let mut changed_rows = try_filled_vec(genes, false, "Cook's changed count rows")?;
    let mut replace = try_filled_vec(genes, None, "Cook's replacement flags")?;
    let mut candidates = try_alloc_matrix(genes, samples, "Cook's replacement candidates")?;
    let mut outliers = try_alloc_matrix(genes, samples, "Cook's outlier mask")?;
    candidates.resize(counts.as_slice().len(), 0);
    outliers.resize(counts.as_slice().len(), false);

    let first_error = candidates
        .par_chunks_mut(samples)
        .zip(outliers.par_chunks_mut(samples))
        .zip(replace.par_iter_mut())
        .zip(changed_rows.par_iter_mut())
        .enumerate()
        .map_init(
            || None,
            |scratch,
             (gene, (((candidate, outlier), replace), changed))|
             -> Result<(), DeseqError> {
                let values = cooks_row_buffer(scratch, samples)?;
                values.clear();
                values.extend_from_slice(normalized_counts.row(gene)?);
                let trim_mean = r_trimmed_mean(values, options.trim);
                let scales = match normalization_factors {
                    Some(factors) => factors.row(gene)?,
                    None => size_factors,
                };
                let counts_row = counts.row_values(gene);
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
                    if outlier[sample]
                        && replaceable_samples[sample]
                        && counts_row[sample] != candidate[sample]
                    {
                        *changed = true;
                    }
                }
                *replace = if has_outlier {
                    Some(true)
                } else if has_missing {
                    None
                } else {
                    Some(false)
                };
                Ok(())
            },
        )
        .find_map_first(Result::err);
    if let Some(error) = first_error {
        return Err(error);
    }

    let replaced_counts = if changed_rows.iter().any(|&changed| changed) {
        let mut replaced = try_alloc_matrix(genes, samples, "Cook's replaced counts")?;
        replaced.extend_from_slice(counts.as_slice());
        replaced
            .par_chunks_mut(samples)
            .zip(candidates.par_chunks(samples))
            .zip(outliers.par_chunks(samples))
            .zip(changed_rows.par_iter())
            .filter(|&(_, &changed)| changed)
            .for_each(|(((replaced, candidates), outliers), _)| {
                for (sample, &replaceable) in replaceable_samples.iter().enumerate() {
                    if replaceable && outliers[sample] {
                        replaced[sample] = candidates[sample];
                    }
                }
            });
        counts.with_counts(replaced)?
    } else {
        counts.clone()
    };

    Ok(CooksReplacementOutput {
        replaced_counts,
        candidate_replacement_counts: counts.with_counts(candidates)?,
        outlier_cells: RowMajorMatrix::from_row_major(
            counts.n_genes(),
            counts.n_samples(),
            outliers,
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
    let offsets = normalization_factors.map_or(Offsets::PerSample(size_factors), Offsets::PerCell);
    prepare_cooks_replacement_refit_from_normalized(
        counts,
        Cow::Borrowed(normalized_counts),
        offsets,
        cooks,
        model_matrix,
        options,
        observation_weights,
    )
}

/// Owned input was normalized from these counts with the same offsets.
pub(crate) fn prepare_cooks_replacement_refit_from_normalized(
    counts: &CountMatrix,
    normalized: Cow<'_, RowMajorMatrix<f64>>,
    offsets: Offsets<'_>,
    cooks: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
    options: &CooksReplacementOptions,
    observation_weights: Option<&RowMajorMatrix<f64>>,
) -> Result<CooksRefitPlan, DeseqError> {
    let (size_factors, normalization_factors) = match offsets {
        Offsets::PerSample(factors) => (factors, None),
        Offsets::PerCell(factors) => (&[][..], Some(factors)),
    };
    let replacement = replace_outlier_counts(
        counts,
        &normalized,
        size_factors,
        normalization_factors,
        cooks,
        model_matrix,
        options,
    )?;
    let replaced_normalized_counts = match normalized {
        Cow::Owned(mut values) => {
            if !counts.shares_count_storage(&replacement.replaced_counts) {
                let marked_rows = replacement
                    .replace
                    .iter()
                    .enumerate()
                    .filter_map(|(gene, &flag)| (flag == Some(true)).then_some(gene));
                renormalize_changed_rows(
                    &mut values,
                    counts,
                    &replacement.replaced_counts,
                    offsets,
                    marked_rows,
                )?;
            }
            values
        }
        Cow::Borrowed(_) => {
            normalized_counts_for_rows(&replacement.replaced_counts, offsets, None)?
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn owned_normalization_reuses_storage_for_weighted_replacement_metadata() {
        let counts =
            CountMatrix::from_row_major_u32(2, 4, vec![10, 200, 10, 10, 5, 5, 5, 5]).unwrap();
        let weights =
            RowMajorMatrix::from_row_major(2, 4, vec![1.0, 0.5, 0.8, 1.0, 1.0, 1.0, 0.5, 0.3])
                .unwrap();
        let size_factors = [1.5, 1.0, 1.25, 0.8];
        let cell_factors =
            RowMajorMatrix::from_row_major(2, 4, vec![0.8, 1.3, 1.2, 1.1, 1.4, 1.0, 0.9, 1.1])
                .unwrap();
        let cooks =
            RowMajorMatrix::from_row_major(2, 4, vec![0., 2., 0., 0., 0., 0., 0., 0.]).unwrap();
        let design = DesignMatrix::from_row_major(4, 1, vec![1.; 4], None).unwrap();
        let mut options = CooksReplacementOptions::new(1.0);
        options.which_samples = Some(vec![false, true, false, false]);

        for factors in [None, Some(&cell_factors)] {
            let offsets = factors.map_or(Offsets::PerSample(&size_factors), Offsets::PerCell);
            let normalized = normalized_counts_for_rows(&counts, offsets, None).unwrap();
            let reference = prepare_cooks_replacement_refit_with_weights(
                &counts,
                &normalized,
                &size_factors,
                factors,
                &cooks,
                &design,
                &options,
                Some(&weights),
            )
            .unwrap();
            let original_buffer = normalized.as_slice().as_ptr();
            let reused = prepare_cooks_replacement_refit_from_normalized(
                &counts,
                Cow::Owned(normalized),
                offsets,
                &cooks,
                &design,
                &options,
                Some(&weights),
            )
            .unwrap();
            assert_eq!(reused, reference);
            assert_ne!(reused.replacement.replaced_counts.row(0).unwrap()[1], 200);
            assert_eq!(
                reused.replaced_normalized_counts.as_slice().as_ptr(),
                original_buffer,
            );
        }
    }

    #[test]
    fn public_planner_recomputes_normalization_when_seed_counts_are_caller_supplied() {
        let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 20, 30, 40]).unwrap();
        let seed = RowMajorMatrix::from_elem(1, 4, 3.0).unwrap();
        let size_factors = [1.0, 2.0, 3.0, 4.0];
        let cooks = RowMajorMatrix::from_elem(1, 4, 0.0).unwrap();
        let design = DesignMatrix::from_row_major(4, 1, vec![1.; 4], None).unwrap();
        let plan = prepare_cooks_replacement_refit(
            &counts,
            &seed,
            &size_factors,
            None,
            &cooks,
            &design,
            &CooksReplacementOptions::new(1.0),
        )
        .unwrap();
        assert_eq!(plan.replaced_normalized_counts.as_slice(), &[10.0; 4]);
        assert_eq!(
            plan.replacement.candidate_replacement_counts.as_slice(),
            &[3, 6, 9, 12]
        );
    }
}
