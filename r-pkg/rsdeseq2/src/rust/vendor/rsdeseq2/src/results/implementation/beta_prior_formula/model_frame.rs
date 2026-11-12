use super::*;

/// Parse a primitive formula from owned model-frame metadata, fit the expanded
/// beta-prior model, and assemble coefficient Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_results(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_results(formula_design, input.into(), coefficient)
}

/// Parse a primitive formula from owned model-frame metadata, fit the expanded
/// beta-prior model, and assemble contrast Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_contrast_results(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results(
        formula_design,
        input.into(),
        contrast,
    )
}

/// Parse a primitive formula from owned model-frame metadata and run
/// coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_results_with_cooks_replacement(
        formula_design,
        input.into(),
        coefficient,
        replacement_options,
    )
}

/// Parse a primitive formula from owned model-frame metadata and run contrast
/// beta-prior Wald replacement refit.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_cooks_replacement(
        formula_design,
        input.into(),
        contrast,
        replacement_options,
    )
}

/// Parse a primitive formula from owned model-frame metadata, use
/// normalization factors, and assemble coefficient Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights(
        formula_design,
        input.into(),
        coefficient,
    )
}

/// Parse a primitive formula from owned model-frame metadata, use
/// normalization factors, and assemble contrast Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
        formula_design,
        input.into(),
        contrast,
    )
}

/// Parse a primitive formula from owned model-frame metadata, use
/// normalization factors, and run coefficient beta-prior Wald replacement
/// refit.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
        formula_design,
        input.into(),
        coefficient,
        replacement_options,
    )
}

/// Parse a primitive formula from owned model-frame metadata, use
/// normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
        formula_design,
        input.into(),
        contrast,
        replacement_options,
    )
}
