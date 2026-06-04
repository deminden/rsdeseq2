/// Parse process arguments and run the CLI.
pub fn run_cli() -> Result<(), DeseqError> {
    run(Cli::parse())
}

fn run(cli: Cli) -> Result<(), DeseqError> {
    match cli.command {
        Commands::SizeFactors {
            counts,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let geometric_means = read_cli_geometric_means(geometric_means, &counts)?;
            let size_factors = estimate_size_factors_with_options(
                &counts,
                method.into(),
                geometric_means.as_deref(),
                control_genes.as_deref(),
            )?;
            write_size_factors_tsv(output, counts.sample_names(), &size_factors)
        }
        Commands::BaseMean {
            counts,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let normalized = cli_normalized_counts(
                &counts,
                normalization_factors,
                size_factors,
                method,
                geometric_means,
                control_genes,
            )?;
            let base_mean = if let Some(path) = observation_weights {
                let weights = read_cli_observation_weights(path, &counts)?;
                base_mean_with_weights(&normalized, &weights)?
            } else {
                base_mean(&normalized)?
            };
            write_base_mean_tsv(output, counts.gene_names(), &base_mean)
        }
        Commands::NormalizedCounts {
            counts,
            normalization_factors,
            size_factors,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let normalized = cli_normalized_counts(
                &counts,
                normalization_factors,
                size_factors,
                method,
                geometric_means,
                control_genes,
            )?;
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &normalized,
            )
        }
        Commands::Vst {
            counts,
            design,
            blind,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            nsub,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
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
            let transformed = if blind {
                builder.blind_vst_glm_mu_auto(&counts, nsub)?.transformed
            } else {
                let design = design.ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "VST design path".to_string(),
                    expected: 1,
                    actual: 0,
                })?;
                let design = read_cli_design_matrix(design, &counts)?;
                builder.vst_glm_mu_auto(&counts, &design, nsub)?.transformed
            };
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &transformed,
            )
        }
        Commands::Rlog {
            counts,
            design,
            blind,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            frozen_intercept,
            rlog_prior_variance,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
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
            let frozen_intercept = read_cli_frozen_intercept(frozen_intercept, &counts)?;
            let transformed = if blind {
                if let Some(frozen_intercept) = frozen_intercept {
                    let prior = required_cli_rlog_prior_variance(rlog_prior_variance)?;
                    let fit = builder.fit_map_dispersions_glm_mu(
                        &counts,
                        &crate::design::DesignMatrix::intercept_only(counts.n_samples())?,
                    )?;
                    fit.frozen_rlog(&counts, &frozen_intercept, prior)?
                        .transformed
                } else {
                    if rlog_prior_variance.is_some() {
                        return Err(cli_rlog_prior_without_frozen_intercept());
                    }
                    builder.blind_rlog_glm_mu(&counts)?.transformed
                }
            } else {
                let design = design.ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "rlog design path".to_string(),
                    expected: 1,
                    actual: 0,
                })?;
                let design = read_cli_design_matrix(design, &counts)?;
                if let Some(frozen_intercept) = frozen_intercept {
                    let prior = required_cli_rlog_prior_variance(rlog_prior_variance)?;
                    let fit = builder.fit_map_dispersions_glm_mu(&counts, &design)?;
                    fit.frozen_rlog(&counts, &frozen_intercept, prior)?
                        .transformed
                } else {
                    if rlog_prior_variance.is_some() {
                        return Err(cli_rlog_prior_without_frozen_intercept());
                    }
                    builder.rlog_glm_mu(&counts, &design)?.transformed
                }
            };
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &transformed,
            )
        }
        Commands::Wald {
            counts,
            design,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            coefficient,
            coefficient_name,
            contrast,
            beta_prior_expanded_design,
            beta_prior_coefficient_groups,
            beta_prior_dispersions,
            beta_prior_base_mean,
            beta_prior_disp_fit,
            beta_prior_factor,
            beta_prior_reference,
            beta_prior_sample_levels,
            beta_prior_additive_factors,
            beta_prior_additive_references,
            beta_prior_additive_sample_levels,
            beta_prior_additive_numeric,
            beta_prior_additive_numeric_values,
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
            lfc_threshold,
            alternative,
            use_t,
            t_degrees_of_freedom,
            t_degrees_of_freedom_file,
            disable_cooks_cutoff,
            cooks_cutoff,
            disable_independent_filtering,
            independent_filtering_alpha,
            independent_filtering_theta,
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
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_cli_design_matrix(design, &counts)?;
            let beta_prior_normalization_factors = normalization_factors.clone();
            let beta_prior_size_factors = size_factors.clone();
            let beta_prior_observation_weights = observation_weights.clone();
            let beta_prior_geometric_means = geometric_means.clone();
            let beta_prior_control_genes = control_genes.clone();
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into())
                .wald_lfc_threshold(lfc_threshold, alternative.into());
            builder = apply_cli_wald_t_options(
                builder,
                &counts,
                use_t,
                t_degrees_of_freedom,
                t_degrees_of_freedom_file,
            )?;
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
                    context: "Wald coefficient and contrast inputs".to_string(),
                    expected: 1,
                    actual: contrast_inputs,
                });
            }
            let cutoff = resolve_cooks_cutoff(
                builder.current_cooks_cutoff(),
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let beta_prior_inputs = [
                beta_prior_expanded_design.is_some(),
                beta_prior_coefficient_groups.is_some(),
                beta_prior_dispersions.is_some(),
                beta_prior_base_mean.is_some(),
                beta_prior_disp_fit.is_some(),
                beta_prior_factor.is_some(),
                beta_prior_reference.is_some(),
                beta_prior_sample_levels.is_some(),
                beta_prior_additive_factors.is_some(),
                beta_prior_additive_references.is_some(),
                beta_prior_additive_sample_levels.is_some(),
                beta_prior_additive_numeric.is_some(),
                beta_prior_additive_numeric_values.is_some(),
            ]
            .into_iter()
            .filter(|present| *present)
            .count();
            let analysis = if beta_prior_inputs > 0 {
                let beta_prior_matrix_inputs =
                    beta_prior_expanded_design.is_some() || beta_prior_coefficient_groups.is_some();
                let beta_prior_factor_inputs = beta_prior_factor.is_some()
                    || beta_prior_reference.is_some()
                    || beta_prior_sample_levels.is_some();
                let beta_prior_additive_inputs = beta_prior_additive_factors.is_some()
                    || beta_prior_additive_references.is_some()
                    || beta_prior_additive_sample_levels.is_some()
                    || beta_prior_additive_numeric.is_some()
                    || beta_prior_additive_numeric_values.is_some();
                let beta_prior_design_routes = [
                    beta_prior_matrix_inputs,
                    beta_prior_factor_inputs,
                    beta_prior_additive_inputs,
                ]
                .into_iter()
                .filter(|present| *present)
                .count();
                if beta_prior_design_routes > 1 {
                    return Err(DeseqError::InvalidOptions {
                        reason:
                            "beta-prior expanded matrix, one-factor, and additive-factor inputs are mutually exclusive"
                                .to_string(),
                    });
                }
                if beta_prior_matrix_inputs && beta_prior_inputs != 5 {
                    return Err(DeseqError::InvalidDimensions {
                        context: "beta-prior expanded Wald inputs".to_string(),
                        expected: 5,
                        actual: beta_prior_inputs,
                    });
                }
                if beta_prior_factor_inputs && beta_prior_inputs != 6 {
                    return Err(DeseqError::InvalidDimensions {
                        context: "beta-prior factor Wald inputs".to_string(),
                        expected: 6,
                        actual: beta_prior_inputs,
                    });
                }
                if beta_prior_additive_inputs {
                    let additive_common_inputs = [
                        beta_prior_dispersions.is_some(),
                        beta_prior_base_mean.is_some(),
                        beta_prior_disp_fit.is_some(),
                    ]
                    .into_iter()
                    .filter(|present| *present)
                    .count();
                    let additive_factor_inputs = [
                        beta_prior_additive_factors.is_some(),
                        beta_prior_additive_references.is_some(),
                        beta_prior_additive_sample_levels.is_some(),
                    ]
                    .into_iter()
                    .filter(|present| *present)
                    .count();
                    let additive_numeric_inputs = [
                        beta_prior_additive_numeric.is_some(),
                        beta_prior_additive_numeric_values.is_some(),
                    ]
                    .into_iter()
                    .filter(|present| *present)
                    .count();
                    if additive_common_inputs != 3
                        || matches!(additive_factor_inputs, 1 | 2)
                        || additive_numeric_inputs == 1
                        || additive_factor_inputs + additive_numeric_inputs == 0
                    {
                        return Err(DeseqError::InvalidDimensions {
                            context: "beta-prior additive-factor Wald inputs".to_string(),
                            expected: 3,
                            actual: additive_common_inputs,
                        });
                    }
                }
                if beta_prior_additive_inputs && beta_prior_inputs < 5 {
                    return Err(DeseqError::InvalidDimensions {
                        context: "beta-prior additive-factor Wald inputs".to_string(),
                        expected: 5,
                        actual: beta_prior_inputs,
                    });
                }
                if contrast_factor.is_some()
                    || contrast_numerator.is_some()
                    || contrast_denominator.is_some()
                    || contrast_reference.is_some()
                    || contrast_sample_levels.is_some()
                {
                    return Err(DeseqError::InvalidOptions {
                        reason:
                            "beta-prior expanded Wald CLI currently accepts coefficient, named, list, or numeric contrasts"
                                .to_string(),
                    });
                }
                if beta_prior_factor_inputs {
                    cli_factor_beta_prior_wald_analysis(
                        &counts,
                        &design,
                        beta_prior_factor.unwrap(),
                        beta_prior_reference.unwrap(),
                        beta_prior_sample_levels.unwrap(),
                        beta_prior_dispersions.unwrap(),
                        beta_prior_base_mean.unwrap(),
                        beta_prior_disp_fit.unwrap(),
                        beta_prior_normalization_factors,
                        beta_prior_size_factors,
                        beta_prior_observation_weights,
                        method,
                        beta_prior_geometric_means,
                        beta_prior_control_genes,
                        coefficient,
                        coefficient_name,
                        contrast,
                        contrast_name,
                        contrast_positive,
                        contrast_negative,
                        contrast_positive_weight,
                        contrast_negative_weight,
                        cutoff,
                    )?
                } else if beta_prior_additive_inputs {
                    cli_additive_beta_prior_wald_analysis(
                        &counts,
                        &design,
                        beta_prior_additive_factors.unwrap_or_default(),
                        beta_prior_additive_references.unwrap_or_default(),
                        beta_prior_additive_sample_levels.unwrap_or_default(),
                        beta_prior_additive_numeric.unwrap_or_default(),
                        beta_prior_additive_numeric_values.unwrap_or_default(),
                        beta_prior_dispersions.unwrap(),
                        beta_prior_base_mean.unwrap(),
                        beta_prior_disp_fit.unwrap(),
                        beta_prior_normalization_factors,
                        beta_prior_size_factors,
                        beta_prior_observation_weights,
                        method,
                        beta_prior_geometric_means,
                        beta_prior_control_genes,
                        coefficient,
                        coefficient_name,
                        contrast,
                        contrast_name,
                        contrast_positive,
                        contrast_negative,
                        contrast_positive_weight,
                        contrast_negative_weight,
                        cutoff,
                    )?
                } else {
                    cli_expanded_beta_prior_wald_analysis(
                        &counts,
                        &design,
                        beta_prior_expanded_design.unwrap(),
                        &beta_prior_coefficient_groups.unwrap(),
                        beta_prior_dispersions.unwrap(),
                        beta_prior_base_mean.unwrap(),
                        beta_prior_disp_fit.unwrap(),
                        beta_prior_normalization_factors,
                        beta_prior_size_factors,
                        beta_prior_observation_weights,
                        method,
                        beta_prior_geometric_means,
                        beta_prior_control_genes,
                        coefficient,
                        coefficient_name,
                        contrast,
                        contrast_name,
                        contrast_positive,
                        contrast_negative,
                        contrast_positive_weight,
                        contrast_negative_weight,
                        cutoff,
                    )?
                }
            } else if let Some(contrast) = contrast {
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_wald_glm_mu_contrast(&counts, &design, &contrast)?)
                }
            } else if let Some(contrast_name) = contrast_name {
                let contrast = ContrastSpec::coefficient_name(contrast_name);
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(
                        builder.fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?,
                    )
                }
            } else if contrast_positive.is_some() || contrast_negative.is_some() {
                let contrast = ContrastSpec::list_with_values(
                    contrast_positive.unwrap_or_default(),
                    contrast_negative.unwrap_or_default(),
                    contrast_positive_weight,
                    contrast_negative_weight,
                );
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(
                        builder.fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?,
                    )
                }
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
                let numerator = contrast_numerator.ok_or_else(|| DeseqError::InvalidDimensions {
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
                let path = contrast_sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "factor-level results contrast requires --contrast-sample-levels"
                        .to_string(),
                })?;
                let levels = align_sample_levels_to_samples(
                    &read_sample_levels_tsv(path)?,
                    counts
                        .sample_names()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "count sample names are required to align sample levels"
                                .to_string(),
                        })?,
                )?;
                let contrast = match contrast_reference {
                    Some(reference) => ResultsContrast::character_with_reference(
                        factor,
                        numerator,
                        denominator,
                        reference,
                    ),
                    None => ResultsContrast::character(factor, numerator, denominator),
                };
                if let Some(cutoff) = cutoff {
                    let output = builder
                        .fit_with_test_results_contrast_request_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            Some(&levels),
                            &CooksReplacementOptions::new(cutoff),
                        )?;
                    let CooksReplacementTestOutput::Wald(output) = output else {
                        unreachable!("Wald CLI branch should route to Wald");
                    };
                    cli_wald_replacement_output(output)
                } else {
                    cli_fit_output(builder.fit_with_results_contrast_request(
                        &counts,
                        &design,
                        &contrast,
                        Some(&levels),
                    )?)
                }
            } else {
                let coefficient = match (coefficient, coefficient_name) {
                    (Some(coefficient), None) => coefficient,
                    (None, Some(name)) => resolve_coefficient_index(&design, &name)?,
                    (None, None) => default_cli_coefficient(&design)?,
                    (Some(_), Some(_)) => unreachable!("checked above"),
                };
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(builder.fit_wald_glm_mu_with_cooks_replacement(
                        &counts,
                        &design,
                        coefficient,
                        &CooksReplacementOptions::new(cutoff),
                    )?)
                } else {
                    cli_fit_output(builder.fit_wald_glm_mu(&counts, &design, coefficient)?)
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
        Commands::Lrt {
            counts,
            design,
            reduced_design,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
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
            disable_cooks_cutoff,
            cooks_cutoff,
            disable_independent_filtering,
            independent_filtering_alpha,
            independent_filtering_theta,
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
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_cli_design_matrix(design, &counts)?;
            let reduced_design = read_cli_design_matrix(reduced_design, &counts)?;
            let mut builder = DeseqBuilder::new()
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
            let analysis = if let Some(contrast) = contrast {
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else if let Some(contrast_name) = contrast_name {
                let contrast = ContrastSpec::coefficient_name(contrast_name);
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast_spec(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else if contrast_positive.is_some() || contrast_negative.is_some() {
                let contrast = ContrastSpec::list_with_values(
                    contrast_positive.unwrap_or_default(),
                    contrast_negative.unwrap_or_default(),
                    contrast_positive_weight,
                    contrast_negative_weight,
                );
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast_spec(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
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
                let numerator = contrast_numerator.ok_or_else(|| DeseqError::InvalidDimensions {
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
                let path = contrast_sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "factor-level results contrast requires --contrast-sample-levels"
                        .to_string(),
                })?;
                let levels = align_sample_levels_to_samples(
                    &read_sample_levels_tsv(path)?,
                    counts
                        .sample_names()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "count sample names are required to align sample levels"
                                .to_string(),
                        })?,
                )?;
                let contrast = match contrast_reference {
                    Some(reference) => ResultsContrast::character_with_reference(
                        factor,
                        numerator,
                        denominator,
                        reference,
                    ),
                    None => ResultsContrast::character(factor, numerator, denominator),
                };
                if let Some(cutoff) = cutoff {
                    let output = builder
                        .clone()
                        .test(TestType::Lrt)
                        .reduced_design(reduced_design.clone())
                        .fit_with_test_results_contrast_request_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            Some(&levels),
                            &CooksReplacementOptions::new(cutoff),
                        )?;
                    let CooksReplacementTestOutput::Lrt(output) = output else {
                        unreachable!("LRT CLI branch should route to LRT");
                    };
                    cli_lrt_replacement_output(output)
                } else {
                    cli_fit_output(builder.fit_lrt_with_results_contrast_request(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                        Some(&levels),
                    )?)
                }
            } else {
                let coefficient = match (coefficient, coefficient_name) {
                    (Some(coefficient), None) => coefficient,
                    (None, Some(name)) => resolve_coefficient_index(&design, &name)?,
                    (None, None) => default_cli_coefficient(&design)?,
                    _ => unreachable!("checked above"),
                };
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(builder.fit_lrt_glm_mu_with_cooks_replacement(
                        &counts,
                        &design,
                        &reduced_design,
                        coefficient,
                        &CooksReplacementOptions::new(cutoff),
                    )?)
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu(
                        &counts,
                        &design,
                        &reduced_design,
                        coefficient,
                    )?)
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
    }
}
