use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};
impl DeseqBuilder {
    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
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
        })
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_contrast`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
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
        })
    }

    /// Run a supplied-dispersion Wald pipeline for a DESeq2 `results(contrast=...)` request.
    ///
    /// Character triplet contrasts require one sample level per count-matrix
    /// column so the Rust core can apply DESeq2's character contrast all-zero
    /// handling. List and numeric contrasts ignore `sample_levels` and use the
    /// numeric all-zero rule.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design).fixed_dispersion(dispersions).result(ResultColumn::Results { contrast, sample_levels, model_frame: None }))"
    )]
    pub fn fit_fixed_dispersion_wald_results_contrast<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                }),
        )?
        .into_results_pair()
    }

    /// Run a supplied-dispersion Wald pipeline for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    ///
    /// Character triplet contrasts resolve their factor reference and
    /// per-sample levels from `model_frame`. List and numeric contrasts use the
    /// same numeric all-zero handling as [`Self::fit_fixed_dispersion_wald_results_contrast`].
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design).fixed_dispersion(dispersions).result(ResultColumn::Results { contrast, sample_levels: None, model_frame: Some(model_frame) }))"
    )]
    pub fn fit_fixed_dispersion_wald_results_contrast_from_model_frame(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                }),
        )?
        .into_results_pair()
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_with_cooks_replacement`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) =
                raw_builder.leaf_fixed_dispersion_wald(counts, design, dispersions, coefficient)?;
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
                let (fit, results) = refit_builder.leaf_fixed_dispersion_wald(
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

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_contrast_with_cooks_replacement`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder
                .leaf_fixed_dispersion_wald_contrast(counts, design, dispersions, contrast)?;
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
                let (fit, results) = refit_builder.leaf_fixed_dispersion_wald_contrast(
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

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let mut output = self.leaf_fixed_dispersion_wald_contrast_with_cooks_replacement(
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
        })
    }

    /// Run supplied-dispersion Wald replacement refit for a DESeq2 `results(contrast=...)` request.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design).fixed_dispersion(dispersions).result(ResultColumn::Results { contrast, sample_levels, model_frame: None }).cooks_replacement(replacement_options))"
    )]
    pub fn fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                })
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run supplied-dispersion Wald replacement refit for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design).fixed_dispersion(dispersions).result(ResultColumn::Results { contrast, sample_levels: None, model_frame: Some(model_frame) }).cooks_replacement(replacement_options))"
    )]
    pub fn fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                })
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_contrast_spec`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let (fit, mut results) = self.leaf_fixed_dispersion_wald_contrast(
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
        })
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_factor_level_contrast`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
            apply_cooks_cutoff_for_design(
                &mut wald_output.results,
                cooks_cutoff,
                counts,
                Some(&wald_output.cooks.cooks),
                design,
                self.model_frame.as_ref(),
                Some(contrast),
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
        })
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
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
            let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder
                .leaf_fixed_dispersion_wald_factor_level_contrast(
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
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (fit, results) = refit_builder
                    .leaf_fixed_dispersion_wald_factor_level_contrast(
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
            let original_cooks = original_fit.cooks.as_ref().ok_or_else(|| {
                missing_prerequisite(ErrorKind::Options, "replacement refit", "Cook's distances")
            })?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                Some(original_cooks),
                design,
                self.model_frame.as_ref(),
                Some(contrast),
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
        })
    }
}
