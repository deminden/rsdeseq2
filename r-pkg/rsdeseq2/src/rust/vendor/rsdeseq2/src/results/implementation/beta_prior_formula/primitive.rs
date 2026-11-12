use super::*;

pub fn fit_expanded_formula_beta_prior_wald_results(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_results(formula_design, input.into(), coefficient)
}

/// Parse a primitive formula, fit the expanded beta-prior model, and assemble contrast rows.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results(
        formula_design,
        input.into(),
        contrast,
    )
}

/// Parse a primitive formula and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_results_with_cooks_replacement(
        formula_design,
        input.into(),
        coefficient,
        replacement_options,
    )
}

/// Parse a primitive formula and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_cooks_replacement(
        formula_design,
        input.into(),
        contrast,
        replacement_options,
    )
}

/// Parse a primitive formula, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights(
        formula_design,
        input.into(),
        coefficient,
    )
}

/// Parse a primitive formula, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(formula_design, input.into(), contrast)
}

/// Parse a primitive formula, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(formula_design, input.into(), coefficient, replacement_options)
}

/// Parse a primitive formula, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(formula_design, input.into(), contrast, replacement_options)
}
