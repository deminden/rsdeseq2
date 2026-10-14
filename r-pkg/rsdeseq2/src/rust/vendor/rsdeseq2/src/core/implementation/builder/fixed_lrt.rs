impl DeseqBuilder {
    /// Run a supplied-dispersion likelihood-ratio test pipeline.
    ///
    /// This mirrors the DESeq2 `nbinomLRT` shape for primitive matrices when
    /// dispersions are already available. The full-model beta fields are
    /// exposed in result rows alongside the model-level LRT statistic and
    /// p-value.
    pub fn fit_fixed_dispersion_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient,
        })?;

        let mut fit = self.base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    /// Run the supplied-dispersion LRT path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_lrt_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder.fit_fixed_dispersion_lrt(
            counts,
            full_design,
            reduced_design,
            dispersions,
            coefficient,
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_lrt(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
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
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        if let Some(contrast) =
            self.model_frame_factor_level_contrast_for_coefficient(full_design, coefficient)?
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

    /// Run a supplied-dispersion likelihood-ratio test and report a numeric contrast.
    ///
    /// This keeps the LRT model comparison unchanged while reporting contrast
    /// estimates and standard errors from the full model in result rows.
    pub fn fit_fixed_dispersion_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &stages.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, full_design, contrast)?;
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &stages.all_zero,
        )?;
        for (gene, all_zero) in stages.all_zero.iter().copied().enumerate() {
            lrt_output.results.rows[gene].max_cooks = lrt_output.cooks.max_cooks[gene];
            if all_zero {
                lrt_output.results.rows[gene].converged = None;
                lrt_output.results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut lrt_output.results, cooks_cutoff)?;
        apply_independent_filtering(&mut lrt_output.results, &self.independent_filtering_options)?;

        let mut fit = self.base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    /// Run the supplied-dispersion LRT numeric-contrast path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_lrt_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder.fit_fixed_dispersion_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            dispersions,
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_lrt_contrast(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
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

    /// Run a supplied-dispersion LRT and report a named full-model contrast.
    pub fn fit_fixed_dispersion_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_fixed_dispersion_lrt_contrast(
            counts,
            full_design,
            reduced_design,
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

    /// Run a supplied-dispersion LRT for a DESeq2 `results(contrast=...)` request.
    ///
    /// Character triplet contrasts require one sample level per count-matrix
    /// column so the Rust core can apply DESeq2's character contrast all-zero
    /// handling. List and numeric contrasts ignore `sample_levels` and use the
    /// numeric all-zero rule. As in DESeq2 LRT result tables, the selected
    /// contrast controls the reported LFC and SE while the statistic and
    /// p-value remain the LRT model-comparison values.
    pub fn fit_fixed_dispersion_lrt_results_contrast<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_fixed_dispersion_lrt_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
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
                self.fit_fixed_dispersion_lrt_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    dispersions,
                    contrast,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_fixed_dispersion_lrt_contrast_spec(
                    counts,
                    full_design,
                    reduced_design,
                    dispersions,
                    &contrast_spec,
                )
            }
        }
    }

    /// Run a supplied-dispersion LRT for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    ///
    /// Character triplet contrasts resolve their factor reference and
    /// per-sample levels from `model_frame`. List and numeric contrasts use the
    /// same numeric all-zero handling as [`Self::fit_fixed_dispersion_lrt_results_contrast`].
    pub fn fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_fixed_dispersion_lrt_factor_level_contrast(
                counts,
                full_design,
                reduced_design,
                dispersions,
                factor_contrast,
            );
        }
        builder.fit_fixed_dispersion_lrt_results_contrast::<String>(
            counts,
            full_design,
            reduced_design,
            dispersions,
            contrast,
            None,
        )
    }

    /// Run supplied-dispersion LRT replacement refit for a named primitive contrast specification.
    pub fn fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let mut output = self.fit_fixed_dispersion_lrt_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            dispersions,
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

    /// Run supplied-dispersion LRT replacement refit for a DESeq2 `results(contrast=...)` request.
    #[allow(clippy::too_many_arguments)]
    pub fn fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
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
                self.fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    dispersions,
                    contrast,
                    replacement_options,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    dispersions,
                    &contrast_spec,
                    replacement_options,
                )
            }
        }
    }

    /// Run supplied-dispersion LRT replacement refit for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    #[allow(clippy::too_many_arguments)]
    pub fn fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
                counts,
                full_design,
                reduced_design,
                dispersions,
                factor_contrast,
                replacement_options,
            );
        }
        builder.fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement::<String>(
            counts,
            full_design,
            reduced_design,
            dispersions,
            contrast,
            None,
            replacement_options,
        )
    }

    /// Run a supplied-dispersion LRT and report a factor-level full-model contrast.
    ///
    /// This resolves DESeq2-shaped coefficient names from the full design
    /// matrix and applies character-style `contrastAllZero` handling from
    /// caller-supplied sample levels. As in DESeq2 LRT result tables, the
    /// all-zero cleanup only zeroes the displayed LFC; the model-comparison
    /// statistic and p-values remain unchanged.
    pub fn fit_fixed_dispersion_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
        let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
        let contrast_all_zero = contrast_all_zero_factor_levels(
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            &numeric_contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &stages.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &stages.all_zero,
        )?;
        for (gene, all_zero) in stages.all_zero.iter().copied().enumerate() {
            lrt_output.results.rows[gene].max_cooks = lrt_output.cooks.max_cooks[gene];
            if all_zero {
                lrt_output.results.rows[gene].converged = None;
                lrt_output.results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut lrt_output.results,
            cooks_cutoff,
            counts,
            &lrt_output.cooks.cooks,
            contrast,
        )?;
        apply_independent_filtering(&mut lrt_output.results, &self.independent_filtering_options)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        lrt_output.results.set_resolved_contrast_metadata(
            result_name,
            comparison,
            &numeric_contrast,
        );

        let mut fit = self.base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    /// Run supplied-dispersion LRT replacement refit for a factor-level full-model contrast.
    pub fn fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
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
        let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder
            .fit_fixed_dispersion_lrt_factor_level_contrast(
                counts,
                full_design,
                reduced_design,
                dispersions,
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
            let (fit, results) = refit_builder.fit_fixed_dispersion_lrt_factor_level_contrast(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
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
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }
}
