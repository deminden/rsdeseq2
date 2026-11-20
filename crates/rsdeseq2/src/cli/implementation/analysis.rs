use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};

pub(super) fn cli_normalized_counts(
    counts: &crate::matrix::CountMatrix,
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
        let size_factors = estimate_size_factors_for_dataset(
            counts,
            method.into(),
            geometric_means.as_deref(),
            control_genes.as_deref(),
        )?;
        normalized_counts(counts, &size_factors)
    }
}

pub(super) fn cli_fit_output((fit, results): (DeseqFit, DeseqResults)) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results,
        fit: Some(fit),
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

pub(super) fn cli_wald_replacement_output(output: CooksReplacementWaldOutput) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: Some(output.original_fit),
        refit: output.refit_fit,
        cooks: None,
        refit_plan: Some(output.refit_plan),
    }
}

pub(super) fn cli_lrt_replacement_output(output: CooksReplacementLrtOutput) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: Some(output.original_fit),
        refit: output.refit_fit,
        cooks: None,
        refit_plan: Some(output.refit_plan),
    }
}

pub(super) fn cli_expanded_beta_prior_output(
    output: ExpandedBetaPriorWaldResults,
) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

pub(super) fn cli_expanded_beta_prior_replacement_output(
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

pub(super) fn cli_factor_beta_prior_output(
    output: ExpandedFactorBetaPriorWaldResults,
) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        fit: None,
        refit: None,
        cooks: None,
        refit_plan: None,
    }
}

pub(super) fn cli_factor_beta_prior_replacement_output(
    output: ExpandedFactorBetaPriorWaldReplacementResults,
) -> CliAnalysisOutput {
    cli_expanded_beta_prior_replacement_output(output.replacement)
}

pub(super) fn cli_additive_beta_prior_output(
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

pub(super) fn cli_additive_beta_prior_replacement_output(
    output: ExpandedAdditiveBetaPriorWaldReplacementResults,
) -> CliAnalysisOutput {
    cli_expanded_beta_prior_replacement_output(output.replacement)
}

pub(super) fn run_cli_model(
    builder: &DeseqBuilder,
    counts: &crate::matrix::CountMatrix,
    design: &crate::design::DesignMatrix,
    reduced_design: Option<&crate::design::DesignMatrix>,
    target: ResultTargetArgs,
    contrast: ContrastArgs,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    let factor_contrast = contrast.has_factor();
    let ResultTargetArgs {
        coefficient,
        coefficient_name,
        contrast: numeric_contrast,
    } = target;
    let ContrastArgs {
        contrast_name,
        contrast_positive,
        contrast_negative,
        contrast_positive_weight,
        contrast_negative_weight,
        contrast_factor,
        contrast_numerator,
        contrast_denominator,
        contrast_reference,
        contrast_sample_levels,
    } = contrast;
    let spec;
    let character;
    let levels;
    let column = if let Some(contrast) = &numeric_contrast {
        ResultColumn::Numeric(contrast)
    } else if let Some(contrast_name) = contrast_name {
        spec = ContrastSpec::coefficient_name(contrast_name);
        ResultColumn::Spec(&spec)
    } else if contrast_positive.is_some() || contrast_negative.is_some() {
        spec = ContrastSpec::list_with_values(
            contrast_positive.unwrap_or_default(),
            contrast_negative.unwrap_or_default(),
            contrast_positive_weight,
            contrast_negative_weight,
        );
        ResultColumn::Spec(&spec)
    } else if factor_contrast {
        let factor_inputs = usize::from(contrast_factor.is_some())
            + usize::from(contrast_numerator.is_some())
            + usize::from(contrast_denominator.is_some());
        let factor = contrast_factor.ok_or_else(|| DeseqError::InvalidDimensions {
            context: "factor-level contrast inputs".to_string(),
            expected: 3,
            actual: factor_inputs,
        })?;
        let numerator = contrast_numerator.ok_or_else(|| DeseqError::InvalidDimensions {
            context: "factor-level contrast inputs".to_string(),
            expected: 3,
            actual: factor_inputs,
        })?;
        let denominator = contrast_denominator.ok_or_else(|| DeseqError::InvalidDimensions {
            context: "factor-level contrast inputs".to_string(),
            expected: 3,
            actual: factor_inputs,
        })?;
        let path = contrast_sample_levels.ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "factor-level results contrast",
                "--contrast-sample-levels",
            )
        })?;
        levels = align_sample_levels_to_samples(
            &read_sample_levels_tsv(path)?,
            counts.sample_names().ok_or_else(|| {
                missing_prerequisite(
                    ErrorKind::Options,
                    "aligning sample levels",
                    "count sample names",
                )
            })?,
        )?;
        character = match contrast_reference {
            Some(reference) => {
                ResultsContrast::character_with_reference(factor, numerator, denominator, reference)
            }
            None => ResultsContrast::character(factor, numerator, denominator),
        };
        ResultColumn::Results {
            contrast: &character,
            sample_levels: Some(&levels),
            model_frame: None,
        }
    } else {
        ResultColumn::Coefficient(match (coefficient, coefficient_name) {
            (Some(coefficient), None) => coefficient,
            (None, Some(name)) => resolve_coefficient_index(design, &name)?,
            (None, None) => default_cli_coefficient(design)?,
            (Some(_), Some(_)) => return Err(invariant("checked above")),
        })
    };
    let replacement = cutoff.map(CooksReplacementOptions::new);
    let mut request = Analysis::new(counts, design)
        .means(MeanModel::GlmMu)
        .result(column);
    if let Some(reduced_design) = reduced_design {
        request = request.lrt(reduced_design);
    }
    if let Some(options) = &replacement {
        request = request.cooks_replacement(options);
    }
    match (builder.run(&request)?, reduced_design.is_some()) {
        (RunOutput::Results { fit, results }, _) => Ok(cli_fit_output((fit, results))),
        (RunOutput::Replacement(CooksReplacementTestOutput::Wald(output)), false) => {
            Ok(cli_wald_replacement_output(output))
        }
        (RunOutput::Replacement(CooksReplacementTestOutput::Lrt(output)), true) => {
            Ok(cli_lrt_replacement_output(output))
        }
        (_, false) => Err(invariant("Wald CLI run should return Wald results")),
        (_, true) => Err(invariant("LRT CLI run should return LRT results")),
    }
}

pub(super) fn validate_cli_result_target(
    kind: &str,
    target: &ResultTargetArgs,
    contrast: &ContrastArgs,
) -> Result<(), DeseqError> {
    let provided = usize::from(target.coefficient.is_some())
        + usize::from(target.coefficient_name.is_some())
        + usize::from(target.contrast.is_some())
        + usize::from(contrast.contrast_name.is_some())
        + usize::from(contrast.contrast_positive.is_some() || contrast.contrast_negative.is_some())
        + usize::from(contrast.has_factor());
    if provided > 1 {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{kind} coefficient and contrast inputs"),
            expected: 1,
            actual: provided,
        });
    }
    Ok(())
}

pub(super) fn configure_cli_results(
    mut builder: DeseqBuilder,
    counts: &crate::matrix::CountMatrix,
    normalization: NormalizationArgs,
    filtering: ResultFilteringArgs,
) -> Result<DeseqBuilder, DeseqError> {
    builder = apply_cli_result_options(
        builder,
        filtering.disable_cooks_cutoff,
        filtering.cooks_cutoff,
        filtering.disable_independent_filtering,
        filtering.independent_filtering_alpha,
        filtering.independent_filtering_theta,
    )?;
    builder = apply_cli_normalization_inputs(
        builder,
        counts,
        normalization.normalization_factors,
        normalization.size_factors,
    )?;
    builder = apply_cli_size_factor_controls(
        builder,
        counts,
        normalization.geometric_means,
        normalization.control_genes,
    )?;
    if let Some(path) = normalization.observation_weights {
        builder = builder.observation_weights(read_cli_observation_weights(path, counts)?);
    }
    Ok(builder)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cli_expanded_beta_prior_wald_analysis(
    counts: &crate::matrix::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    expanded_design: PathBuf,
    coefficient_groups: &str,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization: NormalizationArgs,
    target: ResultTargetArgs,
    contrast: ContrastArgs,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    let NormalizationArgs {
        normalization_factors,
        size_factors,
        observation_weights,
        method,
        geometric_means,
        control_genes,
        ..
    } = normalization;
    let ResultTargetArgs {
        coefficient,
        coefficient_name,
        contrast: numeric_contrast,
    } = target;
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

    let numeric_contrast =
        cli_beta_prior_numeric_contrast(standard_design, numeric_contrast, contrast)?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(standard_design, &name)?),
        (None, None, false) => Some(default_cli_coefficient(standard_design)?),
        (None, None, true) => None,
        _ => return Err(invariant("checked above")),
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        required(coefficient, "--coefficient")?,
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_for_dataset(
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_expanded_beta_prior_output(
                    fit_expanded_beta_prior_wald_results(
                        input,
                        required(coefficient, "--coefficient")?,
                    )?,
                )),
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cli_factor_beta_prior_wald_analysis(
    counts: &crate::matrix::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    factor: String,
    reference: String,
    sample_levels: PathBuf,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization: NormalizationArgs,
    target: ResultTargetArgs,
    contrast: ContrastArgs,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    let NormalizationArgs {
        normalization_factors,
        size_factors,
        observation_weights,
        method,
        geometric_means,
        control_genes,
        ..
    } = normalization;
    let ResultTargetArgs {
        coefficient,
        coefficient_name,
        contrast: numeric_contrast,
    } = target;
    if normalization_factors.is_some() && size_factors.is_some() {
        return Err(cli_conflicting_normalization_inputs());
    }
    let sample_levels = align_sample_levels_to_samples(
        &read_sample_levels_tsv(sample_levels)?,
        counts.sample_names().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "aligning beta-prior sample levels",
                "count sample names",
            )
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
        numeric_contrast,
        contrast,
    )?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(
            &generated_design.standard_design,
            &name,
        )?),
        (None, None, false) => Some(default_cli_coefficient(&generated_design.standard_design)?),
        (None, None, true) => None,
        _ => return Err(invariant("checked above")),
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        required(coefficient, "--coefficient")?,
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_for_dataset(
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_factor_beta_prior_output(
                    fit_expanded_factor_beta_prior_wald_results(input, required(coefficient, "--coefficient")?)?,
                )),
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cli_additive_beta_prior_wald_analysis(
    counts: &crate::matrix::CountMatrix,
    standard_design: &crate::design::DesignMatrix,
    factor_names: Vec<String>,
    references: Vec<String>,
    sample_level_paths: Vec<PathBuf>,
    numeric_names: Vec<String>,
    numeric_value_paths: Vec<PathBuf>,
    dispersions: PathBuf,
    base_mean: PathBuf,
    disp_fit: PathBuf,
    normalization: NormalizationArgs,
    target: ResultTargetArgs,
    contrast: ContrastArgs,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    let NormalizationArgs {
        normalization_factors,
        size_factors,
        observation_weights,
        method,
        geometric_means,
        control_genes,
        ..
    } = normalization;
    let ResultTargetArgs {
        coefficient,
        coefficient_name,
        contrast: numeric_contrast,
    } = target;
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
    let sample_names = counts.sample_names().ok_or_else(|| {
        missing_prerequisite(
            ErrorKind::Options,
            "aligning beta-prior additive sample levels",
            "count sample names",
        )
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
        numeric_contrast,
        contrast,
    )?;
    let coefficient = match (coefficient, coefficient_name, numeric_contrast.is_some()) {
        (Some(coefficient), None, false) => Some(coefficient),
        (None, Some(name), false) => Some(resolve_coefficient_index(
            &generated_design.standard_design,
            &name,
        )?),
        (None, None, false) => Some(default_cli_coefficient(&generated_design.standard_design)?),
        (None, None, true) => None,
        _ => return Err(invariant("checked above")),
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
                        input,
                        required(coefficient, "--coefficient")?,
                    )?,
                )),
            }
        }
        None => {
            let size_factors = if let Some(path) = size_factors {
                read_cli_size_factors(path, counts)?
            } else {
                let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
                estimate_size_factors_for_dataset(
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
                            required(coefficient, "--coefficient")?,
                            &replacement_options,
                        )?,
                    ))
                }
                (None, None) => Ok(cli_additive_beta_prior_output(
                    fit_expanded_additive_beta_prior_wald_results(input, required(coefficient, "--coefficient")?)?,
                )),
            }
        }
    }
}

pub(super) fn cli_beta_prior_numeric_contrast(
    design: &crate::design::DesignMatrix,
    contrast: Option<Vec<f64>>,
    options: ContrastArgs,
) -> Result<Option<Vec<f64>>, DeseqError> {
    let ContrastArgs {
        contrast_name,
        contrast_positive,
        contrast_negative,
        contrast_positive_weight,
        contrast_negative_weight,
        ..
    } = options;
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

pub(super) fn parse_cli_coefficient_groups(
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
