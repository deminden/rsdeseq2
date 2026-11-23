use super::*;
impl DeseqBuilder {
    pub(super) fn leaf_wald_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            validate_pipeline_wald_coefficient(design, coefficient)?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) =
                raw_builder.leaf_wald_glm_mu(counts, design, coefficient)?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let (fit, results) = refit_builder.attach_native_wald(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    coefficient,
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                design.n_samples(),
                design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                design,
                self.model_frame.as_ref(),
                None,
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            Ok(CooksReplacementWaldOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            })
        })
    }

    pub(super) fn leaf_wald_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) =
                raw_builder.leaf_wald_glm_mu_contrast(counts, design, contrast)?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let (fit, results) = refit_builder.attach_native_wald_contrast(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    contrast,
                    None,
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                design.n_samples(),
                design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                design,
                self.model_frame.as_ref(),
                None,
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            Ok(CooksReplacementWaldOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            })
        })
    }

    pub(super) fn leaf_wald_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let mut output = self.leaf_wald_glm_mu_contrast_with_cooks_replacement(
                counts,
                design,
                &numeric_contrast,
                replacement_options,
            )?;
            apply_contrast_metadata_to_replacement_output(
                &mut output,
                contrast.result_name(),
                contrast.comparison(),
                Some(&numeric_contrast),
            );
            Ok(output)
        })
    }

    pub(super) fn leaf_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let metadata_contrast = resolve_contrast(design, &contrast_spec)?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) =
                raw_builder.leaf_wald_glm_mu_factor_level_contrast(counts, design, contrast)?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let contrast_spec = match contrast.reference {
                    Some(reference) => ContrastSpec::factor_level_with_reference(
                        contrast.factor,
                        contrast.numerator,
                        contrast.denominator,
                        reference,
                    ),
                    None => ContrastSpec::factor_level(
                        contrast.factor,
                        contrast.numerator,
                        contrast.denominator,
                    ),
                };
                let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
                let contrast_all_zero = contrast_all_zero_factor_levels(
                    &refit_plan.replacement.replaced_counts,
                    contrast.sample_levels,
                    contrast.numerator,
                    contrast.denominator,
                )?;
                let (fit, results) = refit_builder.attach_native_wald_contrast(
                    &refit_plan.replacement.replaced_counts,
                    design,
                    &numeric_contrast,
                    Some(&contrast_all_zero),
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                design.n_samples(),
                design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                design,
                self.model_frame.as_ref(),
                Some(contrast),
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            let (result_name, comparison) = factor_level_result_metadata(contrast);
            let mut output = CooksReplacementWaldOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            };
            apply_contrast_metadata_to_replacement_output(
                &mut output,
                result_name,
                comparison,
                Some(&metadata_contrast),
            );
            Ok(output)
        })
    }

    pub(super) fn leaf_lrt_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) =
                raw_builder.leaf_lrt_glm_mu(counts, full_design, reduced_design, coefficient)?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                full_design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let (fit, results) = refit_builder.attach_native_lrt(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    reduced_design,
                    coefficient,
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                full_design.n_samples(),
                full_design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                full_design,
                self.model_frame.as_ref(),
                None,
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            Ok(CooksReplacementLrtOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            })
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder.leaf_lrt_glm_mu_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
            )?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                full_design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let (fit, results) = refit_builder.attach_native_lrt_contrast(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    reduced_design,
                    contrast,
                    None,
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                full_design.n_samples(),
                full_design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                full_design,
                self.model_frame.as_ref(),
                None,
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            Ok(CooksReplacementLrtOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            })
        })
    }

    pub(super) fn leaf_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let metadata_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder
                .leaf_lrt_glm_mu_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                )?;
            let refit_plan = replacement_refit_plan_from_original(
                counts,
                full_design,
                &original_fit,
                replacement_options,
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (trend, disp_prior_var, var_log_disp_estimates) =
                    replacement_dispersion_inputs(&original_fit)?;
                let fit = refit_builder.fit_map_dispersions_glm_mu_with_dispersion_function(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    trend,
                    disp_prior_var,
                    var_log_disp_estimates,
                )?;
                let contrast_spec = match contrast.reference {
                    Some(reference) => ContrastSpec::factor_level_with_reference(
                        contrast.factor,
                        contrast.numerator,
                        contrast.denominator,
                        reference,
                    ),
                    None => ContrastSpec::factor_level(
                        contrast.factor,
                        contrast.numerator,
                        contrast.denominator,
                    ),
                };
                let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
                let contrast_all_zero = contrast_all_zero_factor_levels(
                    &refit_plan.replacement.replaced_counts,
                    contrast.sample_levels,
                    contrast.numerator,
                    contrast.denominator,
                )?;
                let (fit, results) = refit_builder.attach_native_lrt_contrast(
                    &refit_plan.replacement.replaced_counts,
                    full_design,
                    reduced_design,
                    &numeric_contrast,
                    Some(&contrast_all_zero),
                    fit,
                )?;
                (Some(fit), Some(results))
            } else {
                (None, None)
            };

            let mut results = merge_replacement_refit_results(
                &original_results,
                refit_results.as_ref(),
                &refit_plan,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                full_design.n_samples(),
                full_design.n_coefficients(),
            )?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                original_fit.cooks.as_ref(),
                full_design,
                self.model_frame.as_ref(),
                Some(contrast),
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

            let (result_name, comparison) = factor_level_result_metadata(contrast);
            let mut output = CooksReplacementLrtOutput {
                original_fit,
                original_results,
                refit_plan,
                refit_fit,
                refit_results,
                results,
            };
            apply_lrt_contrast_metadata_to_replacement_output(
                &mut output,
                result_name,
                comparison,
                Some(&metadata_contrast),
            );
            Ok(output)
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let mut output = self.leaf_lrt_glm_mu_contrast_with_cooks_replacement(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
                replacement_options,
            )?;
            apply_lrt_contrast_metadata_to_replacement_output(
                &mut output,
                contrast.result_name(),
                contrast.comparison(),
                Some(&numeric_contrast),
            );
            Ok(output)
        })
    }
}
