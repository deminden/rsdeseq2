impl DeseqBuilder {
    /// Run the implemented native dispersion path and then a Wald test.
    ///
    /// This is an early, explicitly scoped analogue of `DESeq(..., test="Wald")`:
    /// size factors, base means, linear-mu gene-wise dispersions, selected
    /// trend, deterministic prior variance, MAP dispersions, fixed-dispersion
    /// GLM fitting, Cook's distances, and selected-coefficient Wald results.
    /// Parametric, local, and mean trends are currently implemented. It does not yet
    /// implement DESeq2's general mean/dispersion iteration, beta priors,
    /// contrasts, exact locfit smoothing, or observation weights.
    pub fn fit_wald_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
    }

    /// Run the linear-mu native dispersion path and then a Wald test for a numeric contrast.
    ///
    /// This is the primitive contrast companion to [`Self::fit_wald_linear_mu`].
    /// It reuses the implemented linear-mu dispersion/MAP path, then reports
    /// the requested contrast from the final fixed-dispersion GLM fit.
    pub fn fit_wald_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the linear-mu native Wald path for a named primitive contrast specification.
    pub fn fit_wald_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_linear_mu_contrast(counts, design, &numeric_contrast)?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the linear-mu native Wald path for a caller-supplied factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast(
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
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the parametric native dispersion path and then a Wald test.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_wald_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
    }

    /// Run the parametric linear-mu native Wald path for a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the parametric linear-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_linear_mu_contrast_parametric(counts, design, &numeric_contrast)?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the parametric linear-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast_parametric(
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the GLM-mu native dispersion path and then a Wald test.
    ///
    /// This mirrors `fit_wald_linear_mu` but uses the GLM-mu mean/dispersion
    /// alternation before trend, MAP, fixed-dispersion GLM, Cook's distances,
    /// and selected-coefficient Wald results. Builder-supplied observation
    /// weights are supported for this branch.
    pub fn fit_wald_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
    }

    /// Run the GLM-mu native dispersion path and then a Wald test for a numeric contrast.
    ///
    /// This follows the same native dispersion/MAP and final GLM fitting path
    /// as [`Self::fit_wald_glm_mu`], then reports a primitive numeric
    /// contrast over the fitted coefficient vector. Higher-level formula and
    /// factor handling remains caller or wrapper responsibility.
    pub fn fit_wald_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the GLM-mu native Wald path for a named primitive contrast specification.
    ///
    /// This resolves coefficient names, coefficient-name lists, or supported
    /// factor-level coefficient shapes to a numeric contrast before running
    /// the implemented native GLM-mu Wald contrast path.
    pub fn fit_wald_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_glm_mu_contrast(counts, design, &numeric_contrast)?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the GLM-mu native Wald path for a factor-level contrast.
    ///
    /// In addition to resolving the coefficient contrast, this applies
    /// DESeq2-style character `contrastAllZero` handling using the supplied
    /// per-sample factor levels.
    pub fn fit_wald_glm_mu_factor_level_contrast(
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let fit = raw_builder.fit_map_dispersions_glm_mu(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        let cooks = fit
            .cooks
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "Cook's distances are required before factor-level Cook's filtering"
                    .to_string(),
            })?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut results,
            cooks_cutoff,
            counts,
            cooks,
            contrast,
        )?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the implemented linear-mu native dispersion path and then an LRT.
    ///
    /// This is a limited native analogue of `nbinomLRT`: the full design is
    /// used for the currently implemented linear-mu dispersion/MAP stages, then
    /// full and reduced fixed-dispersion GLMs are fit using those final
    /// dispersions.
    pub fn fit_lrt_linear_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
    }

    /// Run the linear-mu native LRT path and report a full-model numeric contrast.
    ///
    /// The likelihood-ratio statistic and p-values remain the full-vs-reduced
    /// model comparison; the result table's effect-size columns come from the
    /// supplied contrast over the full-model coefficients.
    pub fn fit_lrt_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the linear-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_linear_mu_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the linear-mu native LRT path for a factor-level full-model contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the parametric linear-mu native dispersion path and then an LRT.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_lrt_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
    }

    /// Run the parametric linear-mu native LRT path and report a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the parametric linear-mu native LRT path and report a named contrast.
    pub fn fit_lrt_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_linear_mu_contrast_parametric(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the parametric linear-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the implemented GLM-mu native dispersion path and then an LRT.
    ///
    /// Builder-supplied observation weights are supported for the GLM-mu branch
    /// through the same preprocessing used by native Wald.
    pub fn fit_lrt_glm_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
    }

    /// Run the GLM-mu native LRT path and report a full-model numeric contrast.
    ///
    /// The model comparison, deviance, and p-values come from the LRT. The
    /// result table's effect-size columns use the supplied contrast over the
    /// fitted full-model coefficients, matching DESeq2's result-table shape.
    pub fn fit_lrt_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the GLM-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) =
            self.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, &numeric_contrast)?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the GLM-mu native LRT path for a factor-level full-model contrast.
    ///
    /// In addition to resolving the coefficient contrast, this applies
    /// DESeq2-style character `contrastAllZero` handling using the supplied
    /// per-sample factor levels. For LRT result tables, only the reported LFC
    /// is zeroed; the full-vs-reduced statistic and p-values are preserved.
    pub fn fit_lrt_glm_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let fit = raw_builder.fit_map_dispersions_glm_mu(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        let cooks = fit
            .cooks
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "Cook's distances are required before factor-level Cook's filtering"
                    .to_string(),
            })?;
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut results,
            cooks_cutoff,
            counts,
            cooks,
            contrast,
        )?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }

    /// Run the parametric GLM-mu native dispersion path and then an LRT.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
    }

    /// Run the parametric GLM-mu native LRT path and report a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_glm_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the parametric GLM-mu native LRT path and report a named contrast.
    pub fn fit_lrt_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_glm_mu_contrast_parametric(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.set_resolved_contrast_metadata(
            contrast.result_name(),
            contrast.comparison(),
            &numeric_contrast,
        );
        Ok((fit, results))
    }

    /// Run the parametric GLM-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
        Ok((fit, results))
    }
}
