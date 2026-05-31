pub fn fit_expanded_factor_beta_prior_wald_results(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, fit the beta-prior model, and assemble contrast rows.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
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
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}
