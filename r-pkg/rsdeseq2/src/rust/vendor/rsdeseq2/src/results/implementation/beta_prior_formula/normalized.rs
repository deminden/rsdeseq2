use super::*;

pub(super) fn fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: FormulaPriorNormalizedInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let design = formula_design.design;
    let fit_and_results =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

pub(super) fn fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: FormulaPriorNormalizedInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let design = formula_design.design;
    let fit_and_results =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

pub(super) fn fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: FormulaPriorNormalizedInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let design = formula_design.design;
    let replacement =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

pub(super) fn fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: FormulaPriorNormalizedInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let design = formula_design.design;
    let replacement =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}
