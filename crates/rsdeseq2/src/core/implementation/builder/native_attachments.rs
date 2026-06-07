impl DeseqBuilder {
    /// Run the parametric GLM-mu native dispersion path and then a Wald test.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
    }

    /// Run the parametric GLM-mu native Wald path for a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_glm_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the parametric GLM-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_glm_mu_contrast_parametric(counts, design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the parametric GLM-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.metadata.result_name = Some(result_name);
        results.metadata.comparison = Some(comparison);
        Ok((fit, results))
    }

    fn attach_native_wald(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before Wald fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let wald_output = self.fixed_dispersion_wald_components(WaldPipelineInput {
            counts,
            design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient,
        })?;

        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald);
        Ok((fit, wald_output.results))
    }

    fn attach_native_wald_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before Wald fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let wald_output = self.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &fit.size_factors,
                normalization_factors: fit.normalization_factors.as_ref(),
                observation_weights: fit.observation_weights.as_ref(),
                all_zero: &fit.all_zero,
                dispersions,
            },
            &normalized,
            &fit.base_mean,
            contrast,
            contrast_all_zero_override,
        )?;

        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        Ok((fit, wald_output.results))
    }

    fn attach_native_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before LRT fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient,
        })?;

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

    fn attach_native_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before LRT fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &fit.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        let contrast_all_zero = match contrast_all_zero_override {
            Some(flags) => {
                if flags.len() != counts.n_genes() {
                    return Err(invalid_dimensions(
                        "contrastAllZero rows",
                        counts.n_genes(),
                        flags.len(),
                    ));
                }
                flags.to_vec()
            }
            None => contrast_all_zero_numeric(counts, full_design, contrast)?,
        };
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &fit.all_zero,
        )?;
        for (gene, all_zero) in fit.all_zero.iter().copied().enumerate() {
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
}
