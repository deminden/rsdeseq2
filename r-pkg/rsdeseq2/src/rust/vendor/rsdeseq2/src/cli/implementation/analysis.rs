fn cli_normalized_counts(
    counts: &crate::core::CountMatrix,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
    method: SizeFactorMethodArg,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
) -> Result<crate::matrix::RowMajorMatrix<f64>, DeseqError> {
    if let Some(path) = normalization_factors {
        if size_factors.is_some() {
            return Err(cli_conflicting_normalization_inputs());
        }
        let factors = read_cli_normalization_factors(path, counts)?;
        normalized_counts_with_factors(counts, &factors)
    } else if let Some(path) = size_factors {
        let size_factors = read_cli_size_factors(path, counts)?;
        normalized_counts(counts, &size_factors)
    } else {
        let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
        let size_factors = estimate_size_factors_with_options(
            counts,
            method.into(),
            geometric_means.as_deref(),
            control_genes.as_deref(),
        )?;
        normalized_counts(counts, &size_factors)
    }
}

fn cli_fit_output((fit, results): (DeseqFit, DeseqResults)) -> CliAnalysisOutput {
    let cooks = fit.cooks.clone();
    CliAnalysisOutput {
        results,
        fit: Some(fit),
        refit: None,
        cooks,
        refit_plan: None,
    }
}

fn cli_wald_replacement_output(output: CooksReplacementWaldOutput) -> CliAnalysisOutput {
    let cooks = output.original_fit.cooks.clone();
    CliAnalysisOutput {
        results: output.results,
        fit: Some(output.original_fit),
        refit: output.refit_fit,
        cooks,
        refit_plan: Some(output.refit_plan),
    }
}

fn cli_lrt_replacement_output(output: CooksReplacementLrtOutput) -> CliAnalysisOutput {
    let cooks = output.original_fit.cooks.clone();
    CliAnalysisOutput {
        results: output.results,
        fit: Some(output.original_fit),
        refit: output.refit_fit,
        cooks,
        refit_plan: Some(output.refit_plan),
    }
}

fn cli_expanded_beta_prior_output(output: ExpandedBetaPriorWaldResults) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

fn cli_expanded_beta_prior_replacement_output(
    output: ExpandedBetaPriorWaldReplacementResults,
) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: Some(output.cooks.cooks),
        refit_plan: Some(output.refit_plan),
    }
}

fn cli_factor_beta_prior_output(output: ExpandedFactorBetaPriorWaldResults) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

fn cli_factor_beta_prior_replacement_output(
    output: ExpandedFactorBetaPriorWaldReplacementResults,
) -> CliAnalysisOutput {
    cli_expanded_beta_prior_replacement_output(output.replacement)
}

fn cli_additive_beta_prior_output(
    output: ExpandedAdditiveBetaPriorWaldResults,
) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

fn cli_additive_beta_prior_replacement_output(
    output: ExpandedAdditiveBetaPriorWaldReplacementResults,
) -> CliAnalysisOutput {
    cli_expanded_beta_prior_replacement_output(output.replacement)
}

#[allow(clippy::too_many_arguments)]
fn cli_expanded_beta_prior_wald_analysis(
    counts: &crate::core::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    expanded_design: PathBuf,
    coefficient_groups: &str,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
    observation_weights: Option<PathBuf>,
    method: SizeFactorMethodArg,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
    coefficient: Option<usize>,
    coefficient_name: Option<String>,
    contrast: Option<Vec<f64>>,
    contrast_name: Option<String>,
    contrast_positive: Option<Vec<String>>,
    contrast_negative: Option<Vec<String>>,
    contrast_positive_weight: f64,
    contrast_negative_weight: f64,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    if normalization_factors.is_some() && size_factors.is_some() {
        return Err(cli_conflicting_normalization_inputs());
    }
    let expanded_design = read_cli_design_matrix(expanded_design, counts)?;
    let coefficient_groups = parse_cli_coefficient_groups(
        coefficient_groups,
        standard_design.n_coefficients(),
        expanded_design.n_coefficients(),
    )?;
    let design = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded_design,
        standard_design,
        coefficient_groups: &coefficient_groups,
    };
    let dispersions = read_cli_gene_numeric(dispersions, counts, "beta-prior dispersion")?;
    let base_mean = read_cli_gene_numeric(base_mean, counts, "beta-prior baseMean")?;
    let disp_fit = read_cli_gene_numeric(disp_fit, counts, "beta-prior dispFit")?;
    let weights = observation_weights
        .map(|path| read_cli_observation_weights(path, counts))
        .transpose()?;
    let options = BetaPriorRefitOptions::default();
    let replacement_options = cutoff.map(CooksReplacementOptions::new);

    let numeric_contrast = cli_beta_prior_numeric_contrast(
        standard_design,
        contrast,
        contrast_name,
        contrast_positive,
        contrast_negative,
        contrast_positive_weight,
        contrast_negative_weight,
    )?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(standard_design, &name)?),
        (None, None, false) => Some(default_cli_coefficient(standard_design)?),
        (None, None, true) => None,
        _ => unreachable!("checked above"),
    };

    match normalization_factors {
        Some(path) => {
            let normalization_factors = read_cli_normalization_factors(path, counts)?;
            let input = ExpandedBetaPriorWaldNormalizedResultsInput {
                counts,
                design,
                normalization_factors: &normalization_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_expanded_beta_prior_replacement_output(
                        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                        input,
                        &contrast,
                    )?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_expanded_beta_prior_replacement_output(
                        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        coefficient.unwrap(),
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_with_options(
                    counts,
                    method.into(),
                    geometric_means.as_deref(),
                    control_genes.as_deref(),
                )?
            };
            let input = ExpandedBetaPriorWaldResultsInput {
                counts,
                design,
                size_factors: &size_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_expanded_beta_prior_replacement_output(
                        fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_contrast_results(input, &contrast)?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_expanded_beta_prior_replacement_output(
                        fit_expanded_beta_prior_wald_results_with_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_results(input, coefficient.unwrap())?,
                )),
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn cli_factor_beta_prior_wald_analysis(
    counts: &crate::core::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    factor: String,
    reference: String,
    sample_levels: PathBuf,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
    observation_weights: Option<PathBuf>,
    method: SizeFactorMethodArg,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
    coefficient: Option<usize>,
    coefficient_name: Option<String>,
    contrast: Option<Vec<f64>>,
    contrast_name: Option<String>,
    contrast_positive: Option<Vec<String>>,
    contrast_negative: Option<Vec<String>>,
    contrast_positive_weight: f64,
    contrast_negative_weight: f64,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    if normalization_factors.is_some() && size_factors.is_some() {
        return Err(cli_conflicting_normalization_inputs());
    }
    let sample_levels = align_sample_levels_to_samples(
        &read_sample_levels_tsv(sample_levels)?,
        counts
            .sample_names()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "count sample names are required to align beta-prior sample levels"
                    .to_string(),
            })?,
    )?;
    let generated_design = expanded_factor_design(&factor, &sample_levels, &reference)?;
    if &generated_design.standard_design != standard_design {
        return Err(DeseqError::InvalidOptions {
            reason:
                "reported --design does not match the beta-prior factor design generated from sample levels"
                    .to_string(),
        });
    }
    let dispersions = read_cli_gene_numeric(dispersions, counts, "beta-prior dispersion")?;
    let base_mean = read_cli_gene_numeric(base_mean, counts, "beta-prior baseMean")?;
    let disp_fit = read_cli_gene_numeric(disp_fit, counts, "beta-prior dispFit")?;
    let weights = observation_weights
        .map(|path| read_cli_observation_weights(path, counts))
        .transpose()?;
    let options = BetaPriorRefitOptions::default();
    let replacement_options = cutoff.map(CooksReplacementOptions::new);
    let numeric_contrast = cli_beta_prior_numeric_contrast(
        &generated_design.standard_design,
        contrast,
        contrast_name,
        contrast_positive,
        contrast_negative,
        contrast_positive_weight,
        contrast_negative_weight,
    )?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(
            &generated_design.standard_design,
            &name,
        )?),
        (None, None, false) => Some(default_cli_coefficient(&generated_design.standard_design)?),
        (None, None, true) => None,
        _ => unreachable!("checked above"),
    };

    match normalization_factors {
        Some(path) => {
            let normalization_factors = read_cli_normalization_factors(path, counts)?;
            let input = ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts,
                factor: &factor,
                sample_levels: &sample_levels,
                reference: &reference,
                normalization_factors: &normalization_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_factor_beta_prior_replacement_output(
                        fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                        input,
                        &contrast,
                    )?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_factor_beta_prior_replacement_output(
                        fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        coefficient.unwrap(),
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_with_options(
                    counts,
                    method.into(),
                    geometric_means.as_deref(),
                    control_genes.as_deref(),
                )?
            };
            let input = ExpandedFactorBetaPriorWaldResultsInput {
                counts,
                factor: &factor,
                sample_levels: &sample_levels,
                reference: &reference,
                size_factors: &size_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_factor_beta_prior_replacement_output(
                        fit_expanded_factor_beta_prior_wald_contrast_results_with_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_contrast_results(input, &contrast)?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_factor_beta_prior_replacement_output(
                        fit_expanded_factor_beta_prior_wald_results_with_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_results(input, coefficient.unwrap())?,
                )),
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn cli_additive_beta_prior_wald_analysis(
    counts: &crate::core::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    factor_names: Vec<String>,
    references: Vec<String>,
    sample_level_paths: Vec<PathBuf>,
    numeric_names: Vec<String>,
    numeric_value_paths: Vec<PathBuf>,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
    observation_weights: Option<PathBuf>,
    method: SizeFactorMethodArg,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
    coefficient: Option<usize>,
    coefficient_name: Option<String>,
    contrast: Option<Vec<f64>>,
    contrast_name: Option<String>,
    contrast_positive: Option<Vec<String>>,
    contrast_negative: Option<Vec<String>>,
    contrast_positive_weight: f64,
    contrast_negative_weight: f64,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    if normalization_factors.is_some() && size_factors.is_some() {
        return Err(cli_conflicting_normalization_inputs());
    }
    if factor_names.len() != references.len() || factor_names.len() != sample_level_paths.len() {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior additive factor inputs".to_string(),
            expected: factor_names.len(),
            actual: references.len().max(sample_level_paths.len()),
        });
    }
    if numeric_names.len() != numeric_value_paths.len() {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior additive numeric inputs".to_string(),
            expected: numeric_names.len(),
            actual: numeric_value_paths.len(),
        });
    }
    let sample_names = counts
        .sample_names()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "count sample names are required to align beta-prior additive sample levels"
                .to_string(),
        })?;
    let sample_levels = sample_level_paths
        .iter()
        .map(|path| align_sample_levels_to_samples(&read_sample_levels_tsv(path)?, sample_names))
        .collect::<Result<Vec<_>, _>>()?;
    let factors = factor_names
        .iter()
        .zip(sample_levels.iter())
        .zip(references.iter())
        .map(|((factor, sample_levels), reference)| ExpandedFactorSpec {
            factor,
            sample_levels,
            reference,
            levels: None,
        })
        .collect::<Vec<_>>();
    let numeric_values = numeric_value_paths
        .iter()
        .zip(numeric_names.iter())
        .map(|(path, name)| {
            align_sample_numeric_values_to_samples(
                &read_labeled_sample_numeric_tsv(
                    path,
                    &format!("beta-prior additive numeric covariate {name}"),
                )?,
                sample_names,
                &format!("beta-prior additive numeric covariate {name}"),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let numeric_covariates = numeric_names
        .iter()
        .zip(numeric_values.iter())
        .map(|(name, values)| ExpandedNumericSpec { name, values })
        .collect::<Vec<_>>();
    let generated_design = if numeric_covariates.is_empty() {
        expanded_additive_factor_design(&factors)?
    } else {
        expanded_additive_design(&factors, &numeric_covariates)?
    };
    if &generated_design.standard_design != standard_design {
        return Err(DeseqError::InvalidOptions {
            reason:
                "reported --design does not match the beta-prior additive design generated from sample levels"
                    .to_string(),
        });
    }
    let dispersions = read_cli_gene_numeric(dispersions, counts, "beta-prior dispersion")?;
    let base_mean = read_cli_gene_numeric(base_mean, counts, "beta-prior baseMean")?;
    let disp_fit = read_cli_gene_numeric(disp_fit, counts, "beta-prior dispFit")?;
    let weights = observation_weights
        .map(|path| read_cli_observation_weights(path, counts))
        .transpose()?;
    let options = BetaPriorRefitOptions::default();
    let replacement_options = cutoff.map(CooksReplacementOptions::new);
    let numeric_contrast = cli_beta_prior_numeric_contrast(
        &generated_design.standard_design,
        contrast,
        contrast_name,
        contrast_positive,
        contrast_negative,
        contrast_positive_weight,
        contrast_negative_weight,
    )?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(
            &generated_design.standard_design,
            &name,
        )?),
        (None, None, false) => Some(default_cli_coefficient(&generated_design.standard_design)?),
        (None, None, true) => None,
        _ => unreachable!("checked above"),
    };

    match normalization_factors {
        Some(path) => {
            let normalization_factors = read_cli_normalization_factors(path, counts)?;
            let input = ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_additive_beta_prior_replacement_output(
                        fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                        input,
                        &contrast,
                    )?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_additive_beta_prior_replacement_output(
                        fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        coefficient.unwrap(),
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_with_options(
                    counts,
                    method.into(),
                    geometric_means.as_deref(),
                    control_genes.as_deref(),
                )?
            };
            let input = ExpandedAdditiveBetaPriorWaldResultsInput {
                counts,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                size_factors: &size_factors,
                weights: weights.as_ref(),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: counts.gene_names(),
                options,
            };
            match (numeric_contrast, replacement_options) {
                (Some(contrast), Some(replacement_options)) => {
                    Ok(cli_additive_beta_prior_replacement_output(
                        fit_expanded_additive_beta_prior_wald_contrast_results_with_cooks_replacement(
                            input,
                            &contrast,
                            &replacement_options,
                        )?,
                    ))
                }
                (Some(contrast), None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_contrast_results(input, &contrast)?,
                )),
                (None, Some(replacement_options)) => {
                    Ok(cli_additive_beta_prior_replacement_output(
                        fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement(
                            input,
                            coefficient.unwrap(),
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_results(input, coefficient.unwrap())?,
                )),
            }
        }
    }
}

fn cli_beta_prior_numeric_contrast(
    design: &crate::design::DesignMatrix,
    contrast: Option<Vec<f64>>,
    contrast_name: Option<String>,
    contrast_positive: Option<Vec<String>>,
    contrast_negative: Option<Vec<String>>,
    contrast_positive_weight: f64,
    contrast_negative_weight: f64,
) -> Result<Option<Vec<f64>>, DeseqError> {
    if let Some(contrast) = contrast {
        return Ok(Some(contrast));
    }
    if let Some(contrast_name) = contrast_name {
        return resolve_contrast(design, &ContrastSpec::coefficient_name(contrast_name)).map(Some);
    }
    if contrast_positive.is_some() || contrast_negative.is_some() {
        let contrast = ContrastSpec::list_with_values(
            contrast_positive.unwrap_or_default(),
            contrast_negative.unwrap_or_default(),
            contrast_positive_weight,
            contrast_negative_weight,
        );
        return resolve_contrast(design, &contrast).map(Some);
    }
    Ok(None)
}

fn parse_cli_coefficient_groups(
    raw: &str,
    n_coefficients: usize,
    n_expanded_coefficients: usize,
) -> Result<Vec<Vec<usize>>, DeseqError> {
    let groups = raw
        .split('|')
        .map(|group| {
            let indices = group
                .split(',')
                .map(|value| {
                    let value = value.trim();
                    value
                        .parse::<usize>()
                        .map_err(|_| DeseqError::InvalidOptions {
                            reason: format!("invalid coefficient group index '{value}'"),
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            if indices.is_empty() {
                return Err(DeseqError::InvalidOptions {
                    reason: "coefficient groups must not contain empty groups".to_string(),
                });
            }
            if let Some(index) = indices
                .iter()
                .copied()
                .find(|index| *index >= n_expanded_coefficients)
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "coefficient group index {index} is outside the expanded design columns"
                    ),
                });
            }
            Ok(indices)
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;
    if groups.len() != n_coefficients {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior coefficient groups".to_string(),
            expected: n_coefficients,
            actual: groups.len(),
        });
    }
    let mut seen = HashSet::new();
    if let Some(index) = groups
        .iter()
        .flatten()
        .copied()
        .find(|index| !seen.insert(*index))
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("coefficient group index {index} appears more than once"),
        });
    }
    Ok(groups)
}
