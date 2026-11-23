use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};
impl DeseqBuilder {
    pub(super) fn leaf_fixed_dispersion_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
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
            fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like);
            fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged);
            fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter);
            fit.reduced_mu = Some(lrt_output.reduced_fit.mu);
            fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal);
            attach_glm_fit(&mut fit, lrt_output.full_fit);
            fit.lrt = Some(lrt_output.lrt);
            Ok((fit, lrt_output.results))
        })
    }

    pub(super) fn leaf_fixed_dispersion_lrt_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder.leaf_fixed_dispersion_lrt(
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
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (fit, results) = refit_builder.leaf_fixed_dispersion_lrt(
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

    pub(super) fn leaf_fixed_dispersion_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
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
            let view = crate::glm::inference::GlmInferenceView::from(&lrt_output.full_fit);
            let contrast_output = view.estimates()?.project_contrast(contrast)?;
            lrt_output.results = build_lrt_projection_results_unadjusted_view(
                &stages.base_mean,
                &view,
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
            apply_cooks_cutoff_for_design(
                &mut lrt_output.results,
                cooks_cutoff,
                counts,
                Some(&lrt_output.cooks.cooks),
                full_design,
                self.model_frame.as_ref(),
                None,
            )?;
            apply_independent_filtering(
                &mut lrt_output.results,
                &self.independent_filtering_options,
            )?;

            let mut fit = self.base_fit(
                counts,
                Some(full_design.clone()),
                stages.into_base_fit_input(),
            );
            fit.reduced_design = Some(reduced_design.clone());
            fit.dispersion = Some(lrt_output.expanded_dispersions);
            fit.cooks = Some(lrt_output.cooks.cooks);
            fit.max_cooks = Some(lrt_output.cooks.max_cooks);
            fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like);
            fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged);
            fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter);
            fit.reduced_mu = Some(lrt_output.reduced_fit.mu);
            fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal);
            attach_glm_fit(&mut fit, lrt_output.full_fit);
            fit.lrt = Some(lrt_output.lrt);
            Ok((fit, lrt_output.results))
        })
    }

    pub(super) fn leaf_fixed_dispersion_lrt_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder.leaf_fixed_dispersion_lrt_contrast(
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
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (fit, results) = refit_builder.leaf_fixed_dispersion_lrt_contrast(
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

    pub(super) fn leaf_fixed_dispersion_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let (fit, mut results) = self.leaf_fixed_dispersion_lrt_contrast(
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
        })
    }

    pub(super) fn leaf_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let mut output = self.leaf_fixed_dispersion_lrt_contrast_with_cooks_replacement(
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
        })
    }

    pub(super) fn leaf_fixed_dispersion_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
            let view = crate::glm::inference::GlmInferenceView::from(&lrt_output.full_fit);
            let contrast_output = view.estimates()?.project_contrast(&numeric_contrast)?;
            lrt_output.results = build_lrt_projection_results_unadjusted_view(
                &stages.base_mean,
                &view,
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
            apply_cooks_cutoff_for_design(
                &mut lrt_output.results,
                cooks_cutoff,
                counts,
                Some(&lrt_output.cooks.cooks),
                full_design,
                self.model_frame.as_ref(),
                Some(contrast),
            )?;
            apply_independent_filtering(
                &mut lrt_output.results,
                &self.independent_filtering_options,
            )?;
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
            fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like);
            fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged);
            fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter);
            fit.reduced_mu = Some(lrt_output.reduced_fit.mu);
            fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal);
            attach_glm_fit(&mut fit, lrt_output.full_fit);
            fit.lrt = Some(lrt_output.lrt);
            Ok((fit, lrt_output.results))
        })
    }

    pub(super) fn leaf_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
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
            let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let (original_fit, original_results) = raw_builder
                .leaf_fixed_dispersion_lrt_factor_level_contrast(
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
                self.observation_weights.as_deref(),
            )?;

            let (refit_fit, refit_results) = if refit_plan.requires_refit() {
                let mut refit_builder = raw_builder.clone();
                refit_builder.size_factor_options.supplied_size_factors =
                    Some(original_fit.size_factors.clone());
                let (fit, results) = refit_builder
                    .leaf_fixed_dispersion_lrt_factor_level_contrast(
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
            let original_cooks = original_fit.cooks.as_ref().ok_or_else(|| {
                missing_prerequisite(ErrorKind::Options, "replacement refit", "Cook's distances")
            })?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                Some(original_cooks),
                full_design,
                self.model_frame.as_ref(),
                Some(contrast),
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
        })
    }
}
