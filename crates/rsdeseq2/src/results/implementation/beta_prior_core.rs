/// Return the core DESeq2 `results()` column names currently emitted by Rust.
pub fn deseq2_result_core_column_names() -> &'static [&'static str] {
    &DESEQ2_RESULT_CORE_COLUMNS
}

/// Return optional diagnostic result-column names used by Rust result rows.
pub fn rsdeseq2_result_diagnostic_column_names() -> &'static [&'static str] {
    &RSDESEQ2_RESULT_DIAGNOSTIC_COLUMNS
}

/// Build DESeq2-shaped Wald result rows for one coefficient.
///
/// This mirrors the non-contrast, no-independent-filtering result assembly:
/// `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`.
pub fn build_wald_results(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let wald = wald_test_coefficient(fit, coefficient)?;
    build_wald_results_from_wald(base_mean, fit, coefficient, gene_names, dispersions, &wald)
}

/// Collapse an expanded-model fit and build DESeq2-shaped Wald results.
///
/// This is a primitive result-table companion for the beta-prior expanded
/// model workflow. It performs grouped coefficient/covariance collapse, then
/// reports the requested standard-design coefficient with ordinary Wald
/// statistics and BH adjustment.
pub fn build_wald_results_from_expanded_model_fit(
    base_mean: &[f64],
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let collapsed = collapse_expanded_model_fit(expanded_fit, standard_design, coefficient_groups)?;
    build_wald_results(base_mean, &collapsed, coefficient, gene_names, dispersions)
}

/// Collapse an expanded-model fit and build DESeq2-shaped Wald contrast rows.
///
/// This is the contrast companion to
/// [`build_wald_results_from_expanded_model_fit`]. The supplied contrast is on
/// the collapsed standard-design coefficient scale; the helper propagates the
/// expanded covariance through the grouped coefficient average before computing
/// `c' beta` and `sqrt(c' Sigma c)`.
pub fn build_wald_contrast_results_from_expanded_model_fit(
    base_mean: &[f64],
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
    contrast: &[f64],
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let collapsed = collapse_expanded_model_fit(expanded_fit, standard_design, coefficient_groups)?;
    let contrast = wald_test_contrast(&collapsed, contrast)?;
    build_wald_contrast_results(base_mean, &collapsed, &contrast, gene_names, dispersions)
}

/// Build DESeq2-shaped Wald rows from an expanded beta-prior refit output.
///
/// The helper reports the already-collapsed standard-design prior fit stored in
/// [`ExpandedModelBetaPriorGlmFit`], so callers that use the expanded beta-prior
/// workflow do not need to manually pass the collapsed fit to result assembly.
pub fn build_wald_results_from_expanded_beta_prior_fit(
    base_mean: &[f64],
    fit: &ExpandedModelBetaPriorGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    build_wald_results(
        base_mean,
        &fit.prior_fit,
        coefficient,
        gene_names,
        dispersions,
    )
}

/// Build DESeq2-shaped Wald contrast rows from an expanded beta-prior refit output.
///
/// The supplied contrast is on the collapsed standard-design coefficient scale.
pub fn build_wald_contrast_results_from_expanded_beta_prior_fit(
    base_mean: &[f64],
    fit: &ExpandedModelBetaPriorGlmFit,
    contrast: &[f64],
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let contrast = wald_test_contrast(&fit.prior_fit, contrast)?;
    build_wald_contrast_results(
        base_mean,
        &fit.prior_fit,
        &contrast,
        gene_names,
        dispersions,
    )
}

/// Fit an expanded beta-prior model and assemble Wald rows for one coefficient.
///
/// This is a primitive all-Rust companion for callers that already provide the
/// expanded design, standard design, and coefficient groups.
pub fn fit_expanded_beta_prior_wald_results(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit = fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(
        input.counts,
        input.design,
        BetaPriorSizeFactorWeightInput {
            size_factors: input.size_factors,
            weights: input.weights,
        },
        input.dispersions,
        input.base_mean,
        input.disp_fit,
        input.options,
    )?;
    let results = build_wald_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        coefficient,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior model and assemble Wald rows for a numeric contrast.
///
/// The contrast is on the collapsed standard-design coefficient scale.
pub fn fit_expanded_beta_prior_wald_contrast_results(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit = fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(
        input.counts,
        input.design,
        BetaPriorSizeFactorWeightInput {
            size_factors: input.size_factors,
            weights: input.weights,
        },
        input.dispersions,
        input.base_mean,
        input.disp_fit,
        input.options,
    )?;
    let results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        contrast,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior Wald coefficient workflow with Cook's replacement refit.
///
/// Cook's distances are calculated from the collapsed prior fit on the reported
/// standard-design result. Replacement counts are then refit through the same
/// expanded beta-prior workflow with the original size factors and supplied
/// dispersions.
pub fn fit_expanded_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_results(input.clone(), coefficient)?;
    let cooks = beta_prior_cooks_output(input.counts, input.size_factors, &original.fit)?;
    let normalized = normalized_counts(input.counts, input.size_factors)?;
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        input.size_factors,
        None,
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: &refit_plan.replacement.replaced_counts,
                design: input.design,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: &refit_plan.replaced_base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?)
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit an expanded beta-prior Wald contrast workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_contrast_results(input.clone(), contrast)?;
    let cooks = beta_prior_cooks_output(input.counts, input.size_factors, &original.fit)?;
    let normalized = normalized_counts(input.counts, input.size_factors)?;
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        input.size_factors,
        None,
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: &refit_plan.replacement.replaced_counts,
                design: input.design,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: &refit_plan.replaced_base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?)
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit an expanded beta-prior model with normalization factors and assemble Wald rows.
pub fn fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            input.counts,
            input.design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: input.normalization_factors,
                weights: input.weights,
            },
            input.dispersions,
            input.base_mean,
            input.disp_fit,
            input.options,
        )?;
    let results = build_wald_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        coefficient,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior model with normalization factors and assemble contrast rows.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            input.counts,
            input.design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: input.normalization_factors,
                weights: input.weights,
            },
            input.dispersions,
            input.base_mean,
            input.disp_fit,
            input.options,
        )?;
    let results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        contrast,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit a normalization-factor expanded beta-prior Wald coefficient workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        input.clone(),
        coefficient,
    )?;
    let cooks = beta_prior_normalized_cooks_output(
        input.counts,
        input.normalization_factors,
        &original.fit,
    )?;
    let normalized = normalized_counts_with_factors(input.counts, input.normalization_factors)?;
    let replacement_size_factors = vec![1.0; input.counts.n_samples()];
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        &replacement_size_factors,
        Some(input.normalization_factors),
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(
            fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: &refit_plan.replacement.replaced_counts,
                    design: input.design,
                    normalization_factors: input.normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: &refit_plan.replaced_base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                coefficient,
            )?,
        )
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit a normalization-factor expanded beta-prior Wald contrast workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            input.clone(),
            contrast,
        )?;
    let cooks = beta_prior_normalized_cooks_output(
        input.counts,
        input.normalization_factors,
        &original.fit,
    )?;
    let normalized = normalized_counts_with_factors(input.counts, input.normalization_factors)?;
    let replacement_size_factors = vec![1.0; input.counts.n_samples()];
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        &replacement_size_factors,
        Some(input.normalization_factors),
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(
            fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: &refit_plan.replacement.replaced_counts,
                    design: input.design,
                    normalization_factors: input.normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: &refit_plan.replaced_base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                contrast,
            )?,
        )
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}
