pub fn fit_expanded_formula_beta_prior_wald_results(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
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
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
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
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
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
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
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
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
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
                normalization_factors,
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
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
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
                normalization_factors,
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
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
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
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
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

/// Parse a primitive formula from owned model-frame metadata, fit the expanded
/// beta-prior model, and assemble coefficient Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_results(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_results(formula_design, input, coefficient)
}

/// Parse a primitive formula from owned model-frame metadata, fit the expanded
/// beta-prior model, and assemble contrast Wald rows.
pub fn fit_expanded_formula_model_frame_beta_prior_wald_contrast_results(
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design =
        expanded_formula_design_with_offsets_from_model_frame(input.formula, input.model_frame)?;
    fit_expanded_formula_design_beta_prior_wald_contrast_results(formula_design, input, contrast)
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
        input,
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
        input,
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
        input,
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
        input,
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
        input,
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
        input,
        contrast,
        replacement_options,
    )
}

fn fit_expanded_formula_design_beta_prior_wald_results(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    let design = formula_design.design;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
    let fit_and_results = fit_expanded_beta_prior_wald_results(
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
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

fn fit_expanded_formula_design_beta_prior_wald_contrast_results(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    let design = formula_design.design;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
    let fit_and_results = fit_expanded_beta_prior_wald_contrast_results(
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
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

fn fit_expanded_formula_design_beta_prior_wald_results_with_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    let design = formula_design.design;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
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

fn fit_expanded_formula_design_beta_prior_wald_contrast_results_with_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    let design = formula_design.design;
    if let Some(normalization_factors) = offset_factors.as_ref() {
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
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
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

fn fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
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

fn fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
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

fn fit_expanded_formula_design_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
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

fn fit_expanded_formula_design_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    formula_design: ExpandedFormulaDesignWithOffsets,
    input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'_>,
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
