use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};

// Execute the LRT command after output staging.
pub(super) fn run_lrt(command: Commands) -> Result<(), DeseqError> {
    match command {
        Commands::Lrt {
            counts,
            design,
            reduced_design,
            normalization:
                NormalizationArgs {
                    normalization_factors,
                    size_factors,
                    observation_weights,
                    method,
                    geometric_means,
                    control_genes,
                    fit_type,
                },
            coefficient,
            coefficient_name,
            contrast,
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
            filtering:
                ResultFilteringArgs {
                    disable_cooks_cutoff,
                    cooks_cutoff,
                    disable_independent_filtering,
                    independent_filtering_alpha,
                    independent_filtering_theta,
                },
            outputs:
                DiagnosticOutputArgs {
                    result_column_metadata_output,
                    result_table_metadata_output,
                    independent_filter_metadata_output,
                    independent_filter_num_rej_output,
                    independent_filter_lowess_output,
                    fit_diagnostics_output,
                    refit_diagnostics_output,
                    fit_beta_output,
                    fit_beta_se_output,
                    fit_beta_optim_start_output,
                    refit_beta_output,
                    refit_beta_se_output,
                    refit_beta_optim_start_output,
                    cooks_distance_output,
                    cooks_replacement_metadata_output,
                    cooks_replacement_row_metadata_output,
                    cooks_replaced_counts_output,
                    cooks_candidate_replacement_counts_output,
                    cooks_outlier_cells_output,
                    output,
                },
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_cli_design_matrix(design, &counts)?;
            let reduced_design = read_cli_design_matrix(reduced_design, &counts)?;
            let mut builder = DeseqBuilder::new()
                .automatic_dispersion_means()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
            builder = apply_cli_result_options(
                builder,
                disable_cooks_cutoff,
                cooks_cutoff,
                disable_independent_filtering,
                independent_filtering_alpha,
                independent_filtering_theta,
            )?;
            builder = apply_cli_normalization_inputs(
                builder,
                &counts,
                normalization_factors,
                size_factors,
            )?;
            builder =
                apply_cli_size_factor_controls(builder, &counts, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_cli_observation_weights(path, &counts)?);
            }
            let contrast_inputs = usize::from(coefficient.is_some())
                + usize::from(coefficient_name.is_some())
                + usize::from(contrast.is_some())
                + usize::from(contrast_name.is_some())
                + usize::from(contrast_positive.is_some() || contrast_negative.is_some())
                + usize::from(
                    contrast_factor.is_some()
                        || contrast_numerator.is_some()
                        || contrast_denominator.is_some()
                        || contrast_reference.is_some()
                        || contrast_sample_levels.is_some(),
                );
            if contrast_inputs > 1 {
                return Err(DeseqError::InvalidDimensions {
                    context: "LRT coefficient and contrast inputs".to_string(),
                    expected: 1,
                    actual: contrast_inputs,
                });
            }
            let cutoff = resolve_cooks_cutoff(
                builder.current_cooks_cutoff(),
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let analysis = {
                // One run request: the column comes from whichever contrast
                // or coefficient option was given.
                let spec;
                let character;
                let levels;
                let column = if let Some(contrast) = &contrast {
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
                } else if contrast_factor.is_some()
                    || contrast_numerator.is_some()
                    || contrast_denominator.is_some()
                    || contrast_reference.is_some()
                    || contrast_sample_levels.is_some()
                {
                    let factor_inputs = usize::from(contrast_factor.is_some())
                        + usize::from(contrast_numerator.is_some())
                        + usize::from(contrast_denominator.is_some());
                    let factor = contrast_factor.ok_or_else(|| DeseqError::InvalidDimensions {
                        context: "factor-level contrast inputs".to_string(),
                        expected: 3,
                        actual: factor_inputs,
                    })?;
                    let numerator =
                        contrast_numerator.ok_or_else(|| DeseqError::InvalidDimensions {
                            context: "factor-level contrast inputs".to_string(),
                            expected: 3,
                            actual: factor_inputs,
                        })?;
                    let denominator =
                        contrast_denominator.ok_or_else(|| DeseqError::InvalidDimensions {
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
                        Some(reference) => ResultsContrast::character_with_reference(
                            factor,
                            numerator,
                            denominator,
                            reference,
                        ),
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
                        (None, Some(name)) => resolve_coefficient_index(&design, &name)?,
                        (None, None) => default_cli_coefficient(&design)?,
                        _ => return Err(invariant("checked above")),
                    })
                };
                // The builder selects dispersion means automatically, so the
                // GLM-mu request is DESeq2's default path.
                let replacement = cutoff.map(CooksReplacementOptions::new);
                let mut request = Analysis::new(&counts, &design)
                    .lrt(&reduced_design)
                    .means(MeanModel::GlmMu)
                    .result(column);
                if let Some(options) = &replacement {
                    request = request.cooks_replacement(options);
                }
                match builder.run(&request)? {
                    RunOutput::Results { fit, results } => cli_fit_output((fit, results)),
                    RunOutput::Replacement(CooksReplacementTestOutput::Lrt(output)) => {
                        cli_lrt_replacement_output(output)
                    }
                    _ => return Err(invariant("LRT CLI run should return LRT results")),
                }
            };
            let sidecars = CliCooksOutputPaths {
                cooks_distance: cooks_distance_output,
                replacement_metadata: cooks_replacement_metadata_output,
                replacement_row_metadata: cooks_replacement_row_metadata_output,
                replaced_counts: cooks_replaced_counts_output,
                candidate_replacement_counts: cooks_candidate_replacement_counts_output,
                outlier_cells: cooks_outlier_cells_output,
            };
            let result_sidecars = CliResultSidecarPaths {
                column_metadata: result_column_metadata_output,
                table_metadata: result_table_metadata_output,
                independent_filter_metadata: independent_filter_metadata_output,
                independent_filter_num_rej: independent_filter_num_rej_output,
                independent_filter_lowess: independent_filter_lowess_output,
                fit_diagnostics: fit_diagnostics_output,
                refit_diagnostics: refit_diagnostics_output,
                fit_beta: fit_beta_output,
                fit_beta_se: fit_beta_se_output,
                fit_beta_optim_start: fit_beta_optim_start_output,
                refit_beta: refit_beta_output,
                refit_beta_se: refit_beta_se_output,
                refit_beta_optim_start: refit_beta_optim_start_output,
            };
            write_cli_cooks_outputs(
                &sidecars,
                counts.gene_names(),
                counts.sample_names(),
                &analysis,
            )?;
            write_cli_result_sidecars(&result_sidecars, counts.gene_names(), &analysis)?;
            write_deseq_results_tsv(output, &analysis.results)
        }
        _ => Err(invariant("command dispatch selects the matching variant")),
    }
}
