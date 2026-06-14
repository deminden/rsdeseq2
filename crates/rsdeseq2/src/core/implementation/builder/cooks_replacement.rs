impl DeseqBuilder {
    /// Run the current GLM-mu native Wald path with limited Cook's replacement refit.
    ///
    /// This is an explicitly scoped analogue of the `replaceOutliers` /
    /// `refitWithoutOutliers` part of DESeq2 for the currently implemented
    /// native GLM-mu Wald branch. It preserves the original size factors,
    /// builds replacement counts from original Cook's distances, reruns the
    /// implemented GLM-mu dispersion/MAP/Wald path on replacement counts, and
    /// merges only rows marked for refit. It does not yet implement beta
    /// priors or Bioconductor object slots.
    pub fn fit_wald_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_wald_glm_mu(counts, design, coefficient)?;
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        let normalized = match original_fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &original_fit.size_factors)?,
        };
        let refit_plan = prepare_cooks_replacement_refit(
            counts,
            &normalized,
            &original_fit.size_factors,
            original_fit.normalization_factors.as_ref(),
            original_cooks,
            design,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        if let Some(contrast) =
            self.model_frame_factor_level_contrast_for_coefficient(design, coefficient)?
        {
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                counts,
                original_cooks,
                contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }

    /// Run the current GLM-mu native Wald contrast path with limited Cook's replacement refit.
    ///
    /// This mirrors [`Self::fit_wald_glm_mu_with_cooks_replacement`] for
    /// primitive numeric contrasts: original and replacement-count refits both
    /// use the native GLM-mu contrast path, then only rows marked by the
    /// replacement plan are merged into the final result table.
    pub fn fit_wald_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_wald_glm_mu_contrast(counts, design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        if let Some(factor_contrast) =
            self.model_frame_factor_level_contrast_for_numeric_contrast(design, contrast)?
        {
            let original_cooks =
                original_fit
                    .cooks
                    .as_ref()
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "Cook's distances are required before replacement refit"
                            .to_string(),
                    })?;
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                counts,
                original_cooks,
                factor_contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }

    /// Run native GLM-mu Wald replacement refit for a named primitive contrast specification.
    pub fn fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let mut output = self.fit_wald_glm_mu_contrast_with_cooks_replacement(
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
    }

    /// Run native GLM-mu Wald replacement refit for a caller-supplied factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
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
            raw_builder.fit_wald_glm_mu_factor_level_contrast(counts, design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut results,
            cooks_cutoff,
            counts,
            original_cooks,
            contrast,
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
    }

    /// Run the current GLM-mu native LRT path with limited Cook's replacement refit.
    ///
    /// This mirrors the scoped Wald replacement-refit path for the implemented
    /// native GLM-mu LRT branch: first fit on original counts with Cook's
    /// filtering disabled, replace eligible Cook's outlier counts, rerun the
    /// native GLM-mu dispersion/MAP/LRT path on replacement counts with the
    /// original size factors, and merge only rows marked by the refit plan.
    /// Broader DESeq2 behavior for contrasts, beta priors, and Bioconductor
    /// object metadata remains future work.
    pub fn fit_lrt_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)?;
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        let normalized = match original_fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &original_fit.size_factors)?,
        };
        let refit_plan = prepare_cooks_replacement_refit(
            counts,
            &normalized,
            &original_fit.size_factors,
            original_fit.normalization_factors.as_ref(),
            original_cooks,
            full_design,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        if let Some(contrast) =
            self.model_frame_factor_level_contrast_for_coefficient(full_design, coefficient)?
        {
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                counts,
                original_cooks,
                contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }

    /// Run the current GLM-mu native LRT contrast path with limited Cook's replacement refit.
    pub fn fit_lrt_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            full_design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        if let Some(factor_contrast) =
            self.model_frame_factor_level_contrast_for_numeric_contrast(full_design, contrast)?
        {
            let original_cooks =
                original_fit
                    .cooks
                    .as_ref()
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "Cook's distances are required before replacement refit"
                            .to_string(),
                    })?;
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                counts,
                original_cooks,
                factor_contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }

    /// Run native GLM-mu LRT replacement refit for a caller-supplied factor-level contrast.
    pub fn fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
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
        let (original_fit, original_results) = raw_builder.fit_lrt_glm_mu_factor_level_contrast(
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
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
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
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut results,
            cooks_cutoff,
            counts,
            original_cooks,
            contrast,
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
    }

    /// Run native GLM-mu LRT replacement refit for a named primitive contrast specification.
    pub fn fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let mut output = self.fit_lrt_glm_mu_contrast_with_cooks_replacement(
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
    }
}
