impl DeseqBuilder {
    /// Run a supplied-dispersion Wald pipeline for one coefficient.
    ///
    /// This is the current closest analogue to the core `nbinomWaldTest` path,
    /// but it requires caller-supplied dispersions and does not yet implement
    /// contrasts or beta priors.
    /// All-zero rows are skipped during GLM fitting and expanded back as
    /// missing numeric values, matching DESeq2's `buildMatrixWithNARows`
    /// pattern.
    pub fn fit_fixed_dispersion_wald(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let stages = self.normalization_stages_for_design(counts, design)?;
        let wald_output = self.fixed_dispersion_wald_components(WaldPipelineInput {
            counts,
            design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient,
        })?;
        let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald);
        Ok((fit, wald_output.results))
    }

    /// Run a supplied-dispersion Wald pipeline for a primitive numeric contrast.
    ///
    /// The contrast vector must contain one finite value per design
    /// coefficient. This is the Rust primitive analogue of DESeq2's numeric
    /// contrast path after R has resolved formula terms and coefficient names.
    pub fn fit_fixed_dispersion_wald_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let stages = self.normalization_stages_for_design(counts, design)?;
        let wald_output = self.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                observation_weights: stages.observation_weights.as_ref(),
                all_zero: &stages.all_zero,
                dispersions,
            },
            &stages.normalized,
            &stages.base_mean,
            contrast,
            None,
        )?;
        let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        Ok((fit, wald_output.results))
    }

    /// Run a supplied-dispersion Wald pipeline for a DESeq2 `results(contrast=...)` request.
    ///
    /// Character triplet contrasts require one sample level per count-matrix
    /// column so the Rust core can apply DESeq2's character contrast all-zero
    /// handling. List and numeric contrasts ignore `sample_levels` and use the
    /// numeric all-zero rule.
    pub fn fit_fixed_dispersion_wald_results_contrast<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_fixed_dispersion_wald_factor_level_contrast(
                    counts,
                    design,
                    dispersions,
                    factor_contrast,
                );
            }
        match contrast {
            ResultsContrast::Character {
                factor,
                numerator,
                denominator,
                reference,
            } => {
                let levels = sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "character results contrast requires sample levels for contrastAllZero"
                        .to_string(),
                })?;
                let levels = levels
                    .iter()
                    .map(|level| level.as_ref().to_string())
                    .collect::<Vec<_>>();
                let contrast = factor_level_contrast_from_sample_levels(
                    factor,
                    numerator,
                    denominator,
                    reference.as_deref(),
                    &levels,
                )?;
                self.fit_fixed_dispersion_wald_factor_level_contrast(
                    counts,
                    design,
                    dispersions,
                    contrast,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let resolved = resolve_results_contrast(design, contrast)?;
                let stages = self.normalization_stages_for_design(counts, design)?;
                let wald_output = self.fixed_dispersion_wald_contrast_components(
                    FixedDispersionGlmInput {
                        counts,
                        design,
                        size_factors: &stages.size_factors,
                        normalization_factors: stages.normalization_factors.as_ref(),
                        observation_weights: stages.observation_weights.as_ref(),
                        all_zero: &stages.all_zero,
                        dispersions,
                    },
                    &stages.normalized,
                    &stages.base_mean,
                    &resolved.numeric,
                    None,
                )?;
                let mut fit =
                    self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
                fit.dispersion = Some(wald_output.expanded_dispersions);
                fit.cooks = Some(wald_output.cooks.cooks);
                fit.max_cooks = Some(wald_output.cooks.max_cooks);
                attach_glm_fit(&mut fit, wald_output.glm_fit);
                fit.wald = Some(wald_output.wald_contrast.wald);
                let mut results = wald_output.results;
                results.set_resolved_contrast_metadata(
                    resolved.result_name,
                    resolved.comparison,
                    &resolved.numeric,
                );
                Ok((fit, results))
            }
        }
    }

    /// Run a supplied-dispersion Wald pipeline for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    ///
    /// Character triplet contrasts resolve their factor reference and
    /// per-sample levels from `model_frame`. List and numeric contrasts use the
    /// same numeric all-zero handling as [`Self::fit_fixed_dispersion_wald_results_contrast`].
    pub fn fit_fixed_dispersion_wald_results_contrast_from_model_frame(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_fixed_dispersion_wald_factor_level_contrast(
                counts,
                design,
                dispersions,
                factor_contrast,
            );
        }
        builder.fit_fixed_dispersion_wald_results_contrast::<String>(
            counts,
            design,
            dispersions,
            contrast,
            None,
        )
    }

    /// Run the supplied-dispersion Wald coefficient path with Cook's replacement refit.
    ///
    /// This extends the replacement/refit machinery to the fixed-dispersion
    /// core path. Dispersions remain caller-supplied for both the original fit
    /// and the replacement-count refit; no native dispersion re-estimation or
    /// beta-prior refit is performed here.
    pub fn fit_fixed_dispersion_wald_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_fixed_dispersion_wald(counts, design, dispersions, coefficient)?;
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_wald(
                &refit_plan.replacement.replaced_counts,
                design,
                dispersions,
                coefficient,
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
            self.model_frame_factor_level_contrast_for_coefficient(design, coefficient)?
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

    /// Run the supplied-dispersion Wald numeric-contrast path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder.fit_fixed_dispersion_wald_contrast(
            counts,
            design,
            dispersions,
            contrast,
        )?;
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_wald_contrast(
                &refit_plan.replacement.replaced_counts,
                design,
                dispersions,
                contrast,
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

    /// Run supplied-dispersion Wald replacement refit for a named primitive contrast specification.
    pub fn fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let mut output = self.fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
            counts,
            design,
            dispersions,
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

    /// Run supplied-dispersion Wald replacement refit for a DESeq2 `results(contrast=...)` request.
    pub fn fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    dispersions,
                    factor_contrast,
                    replacement_options,
                );
            }
        match contrast {
            ResultsContrast::Character {
                factor,
                numerator,
                denominator,
                reference,
            } => {
                let levels = sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "character results contrast requires sample levels for contrastAllZero"
                        .to_string(),
                })?;
                let levels = levels
                    .iter()
                    .map(|level| level.as_ref().to_string())
                    .collect::<Vec<_>>();
                let contrast = factor_level_contrast_from_sample_levels(
                    factor,
                    numerator,
                    denominator,
                    reference.as_deref(),
                    &levels,
                )?;
                self.fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    dispersions,
                    contrast,
                    replacement_options,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    dispersions,
                    &contrast_spec,
                    replacement_options,
                )
            }
        }
    }

    /// Run supplied-dispersion Wald replacement refit for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    pub fn fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
                counts,
                design,
                dispersions,
                factor_contrast,
                replacement_options,
            );
        }
        builder.fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement::<String>(
            counts,
            design,
            dispersions,
            contrast,
            None,
            replacement_options,
        )
    }

    /// Run a supplied-dispersion Wald pipeline for a named primitive contrast specification.
    ///
    /// This resolves coefficient names and DESeq2-style positive/negative
    /// coefficient-name lists to a numeric contrast before calling
    /// `fit_fixed_dispersion_wald_contrast`.
    pub fn fit_fixed_dispersion_wald_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) = self.fit_fixed_dispersion_wald_contrast(
            counts,
            design,
            dispersions,
            &numeric_contrast,
        )?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run a supplied-dispersion Wald pipeline for a factor-level contrast.
    ///
    /// This resolves DESeq2-shaped coefficient names from the design matrix and
    /// applies DESeq2-style character contrast all-zero handling using
    /// caller-supplied sample levels. Model-frame callers can use
    /// [`Self::fit_fixed_dispersion_wald_results_contrast_from_model_frame`].
    pub fn fit_fixed_dispersion_wald_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let result_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let mut wald_output = result_builder.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                observation_weights: stages.observation_weights.as_ref(),
                all_zero: &stages.all_zero,
                dispersions,
            },
            &stages.normalized,
            &stages.base_mean,
            &numeric_contrast,
            Some(&contrast_all_zero),
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut wald_output.results,
            cooks_cutoff,
            counts,
            &wald_output.cooks.cooks,
            contrast,
        )?;
        apply_independent_filtering(
            &mut wald_output.results,
            &self.independent_filtering_options,
        )?;
        let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        wald_output.results.set_resolved_contrast_metadata(
            result_name,
            comparison,
            &numeric_contrast,
        );
        Ok((fit, wald_output.results))
    }

    /// Run supplied-dispersion Wald replacement refit for a factor-level contrast.
    pub fn fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
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
        let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder
            .fit_fixed_dispersion_wald_factor_level_contrast(
                counts,
                design,
                dispersions,
                contrast,
            )?;
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_wald_factor_level_contrast(
                &refit_plan.replacement.replaced_counts,
                design,
                dispersions,
                contrast,
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
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }
}
