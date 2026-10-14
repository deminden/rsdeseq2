pub fn fit_expanded_additive_beta_prior_wald_results(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, fit the beta-prior model, and assemble contrast rows.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
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

/// Build an additive-factor expanded design and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement = fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
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

/// Build an additive-factor expanded design, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
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

/// Build an additive-factor expanded design, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
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
