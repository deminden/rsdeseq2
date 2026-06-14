impl DeseqBuilder {
    /// Run the currently implemented DESeq-like workflow.
    ///
    /// For `test=Wald`, this follows the implemented GLM-mu native path and
    /// reports the last design coefficient, matching DESeq2's default
    /// coefficient selection shape. For `test=Lrt`, callers must first store a
    /// reduced design with [`DeseqBuilder::reduced_design`].
    pub fn fit(&self, counts: &CountMatrix, design: &DesignMatrix) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results(counts, design)
            .map(|(fit, _results)| fit)
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow, returning only the
    /// fit state.
    pub fn fit_formula(
        &self,
        counts: &CountMatrix,
        formula: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_formula_with_results(counts, formula)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented DESeq-like workflow and return result rows.
    ///
    /// This is the result-table-producing companion to [`DeseqBuilder::fit`].
    pub fn fit_with_results(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results(counts, design, reduced_design)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow with result rows.
    ///
    pub fn fit_formula_with_results(
        &self,
        counts: &CountMatrix,
        formula: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results(counts, &formula_design.design.standard_design)
    }

    /// Run the top-level Wald workflow and report a named design coefficient.
    pub fn fit_with_results_name(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = resolve_coefficient_index(design, coefficient_name)?;
                self.fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_name(counts, design, reduced_design, coefficient_name)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient.
    pub fn fit_formula_with_results_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results_name(
                counts,
                &formula_design.design.standard_design,
                coefficient_name,
            )
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient, returning only the fit state.
    pub fn fit_formula_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_formula_with_results_name(counts, formula, coefficient_name)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.fit_wald_glm_mu_with_cooks_replacement(
                    counts,
                    design,
                    coefficient,
                    replacement_options,
                )
            }
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            }),
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_formula_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results_with_cooks_replacement(
                counts,
                &formula_design.design.standard_design,
                replacement_options,
            )
    }

    /// Run the top-level named Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = resolve_coefficient_index(design, coefficient_name)?;
                self.fit_wald_glm_mu_with_cooks_replacement(
                    counts,
                    design,
                    coefficient,
                    replacement_options,
                )
            }
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            }),
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level named Wald workflow with limited Cook's replacement
    /// refit.
    pub fn fit_formula_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results_name_with_cooks_replacement(
                counts,
                &formula_design.design.standard_design,
                coefficient_name,
                replacement_options,
            )
    }

    /// Run the top-level workflow with limited Cook's replacement refit.
    ///
    /// Unlike [`Self::fit_with_results_with_cooks_replacement`], this method
    /// returns an enum so `test=Lrt` can route through a stored reduced design.
    pub fn fit_with_test_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_with_cooks_replacement(counts, design, replacement_options)
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level workflow with limited Cook's replacement refit.
    pub fn fit_formula_with_test_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_test_results_with_cooks_replacement(
                counts,
                &formula_design.design.standard_design,
                replacement_options,
            )
    }

    /// Run the top-level named workflow with limited Cook's replacement refit.
    ///
    /// The returned enum keeps Wald and LRT replacement output types explicit
    /// while allowing `test=Lrt` to use the builder's stored reduced design.
    pub fn fit_with_test_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_name_with_cooks_replacement(
                    counts,
                    design,
                    coefficient_name,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_name_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    coefficient_name,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level named workflow with limited Cook's replacement refit.
    pub fn fit_formula_with_test_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_test_results_name_with_cooks_replacement(
                counts,
                &formula_design.design.standard_design,
                coefficient_name,
                replacement_options,
            )
    }

    /// Run the currently implemented top-level Wald workflow for a numeric contrast.
    ///
    /// This is the primitive contrast companion to [`Self::fit_with_results`].
    /// It follows the implemented GLM-mu native Wald path when `test=Wald`.
    /// Top-level LRT remains explicit because a contrast is not an LRT reduced
    /// model.
    pub fn fit_with_results_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast(counts, design, reduced_design, contrast)
            }
        }
    }

    /// Run the currently implemented top-level Wald workflow for a numeric contrast.
    pub fn fit_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast".to_string(),
            }),
        }
    }

    /// Run the top-level numeric-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_spec(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_spec(counts, design, reduced_design, contrast)
            }
        }
    }

    /// Run the top-level workflow for a DESeq2 `results(contrast=...)` request.
    ///
    /// Character triplet contrasts require one sample level per count-matrix
    /// column so the Rust core can apply DESeq2's character contrast all-zero
    /// handling. List and numeric contrasts ignore `sample_levels` and use the
    /// numeric all-zero rule.
    pub fn fit_with_results_contrast_request<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if sample_levels.is_none() {
            if let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_with_results_factor_level_contrast(counts, design, factor_contrast);
            }
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
                self.fit_with_results_factor_level_contrast(counts, design, contrast)
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_with_results_contrast_spec(counts, design, &contrast_spec)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run a DESeq2 `results(contrast=...)` request.
    pub fn fit_formula_with_results_contrast_request(
        &self,
        counts: &CountMatrix,
        formula: &str,
        contrast: &ResultsContrast,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results_contrast_request::<String>(
                counts,
                &formula_design.design.standard_design,
                contrast,
                None,
            )
    }

    /// Run the top-level workflow for a DESeq2 `results(contrast=...)`
    /// request using formula model-frame metadata.
    ///
    /// Character triplet contrasts resolve their factor reference and
    /// per-sample levels from `model_frame`. List and numeric contrasts use the
    /// same numeric all-zero handling as [`Self::fit_with_results_contrast_request`].
    pub fn fit_with_results_contrast_request_from_model_frame(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_with_results_factor_level_contrast(counts, design, factor_contrast);
        }
        builder.fit_with_results_contrast_request::<String>(counts, design, contrast, None)
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_contrast_spec(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level named Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast specification"
                    .to_string(),
            }),
        }
    }

    /// Run the top-level named-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run a DESeq2 `results(contrast=...)` request with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_request_with_cooks_replacement<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        if sample_levels.is_none() {
            if let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self
                    .fit_with_test_results_factor_level_contrast_with_cooks_replacement(
                        counts,
                        design,
                        factor_contrast,
                        replacement_options,
                    );
            }
        }
        match self.test {
            TestType::Wald => match contrast {
                ResultsContrast::Character {
                    factor,
                    numerator,
                    denominator,
                    reference,
                } => {
                    let levels = sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                        reason:
                            "character results contrast requires sample levels for contrastAllZero"
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
                    self.fit_with_results_factor_level_contrast_with_cooks_replacement(
                        counts,
                        design,
                        contrast,
                        replacement_options,
                    )
                    .map(CooksReplacementTestOutput::Wald)
                }
                ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                    let contrast_spec = contrast.as_contrast_spec();
                    self.fit_with_results_contrast_spec_with_cooks_replacement(
                        counts,
                        design,
                        &contrast_spec,
                        replacement_options,
                    )
                    .map(CooksReplacementTestOutput::Wald)
                }
            },
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_request_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    sample_levels,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run a DESeq2 `results(contrast=...)` request with limited Cook's
    /// replacement refit.
    pub fn fit_formula_with_test_results_contrast_request_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        contrast: &ResultsContrast,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_test_results_contrast_request_with_cooks_replacement::<String>(
                counts,
                &formula_design.design.standard_design,
                contrast,
                None,
                replacement_options,
            )
    }

    /// Run a DESeq2 `results(contrast=...)` request with limited Cook's
    /// replacement refit using formula model-frame metadata.
    pub fn fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        match self.test {
            TestType::Wald => {
                if let Some(factor_contrast) =
                    factor_level_contrast_from_model_frame(contrast, model_frame)?
                {
                    return builder
                        .fit_with_results_factor_level_contrast_with_cooks_replacement(
                            counts,
                            design,
                            factor_contrast,
                            replacement_options,
                        )
                        .map(CooksReplacementTestOutput::Wald);
                }
                builder.fit_with_test_results_contrast_request_with_cooks_replacement::<String>(
                    counts,
                    design,
                    contrast,
                    None,
                    replacement_options,
                )
            }
            TestType::Lrt => {
                let reduced_design = builder.reduced_design_for_top_level_lrt()?;
                builder
                    .fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    model_frame,
                    replacement_options,
                )
                    .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_factor_level_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_factor_level_contrast(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                )
            }
        }
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_factor_level_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level factor-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald factor-level contrast"
                    .to_string(),
            }),
        }
    }

    /// Run the top-level factor-level contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the currently implemented top-level LRT workflow with a reduced design.
    ///
    /// This follows the implemented GLM-mu native LRT path and reports the last
    /// full-design coefficient by default.
    pub fn fit_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results(counts, full_design, reduced_design)
            .map(|(fit, _results)| fit)
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow, returning only
    /// the fit state.
    pub fn fit_lrt_formula(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_formula_with_results(counts, full_formula, reduced_formula)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and return result rows.
    pub fn fit_lrt_with_results(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let coefficient = default_results_coefficient(full_design)?;
        self.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow.
    ///
    pub fn fit_lrt_formula_with_results(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
            )
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    pub fn fit_lrt_with_results_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let coefficient = resolve_coefficient_index(full_design, coefficient_name)?;
        self.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata, run LRT, and report a named full-design
    /// coefficient effect.
    pub fn fit_lrt_formula_with_results_name(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_name(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                coefficient_name,
            )
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT replacement-refit
    /// workflow for the default reported full-design coefficient.
    pub fn fit_lrt_formula_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_with_cooks_replacement(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                replacement_options,
            )
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT replacement-refit
    /// workflow for a named full-design coefficient effect.
    pub fn fit_lrt_formula_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_name_with_cooks_replacement(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                coefficient_name,
                replacement_options,
            )
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    pub fn fit_lrt_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_name(counts, full_design, reduced_design, coefficient_name)
            .map(|(fit, _results)| fit)
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and report a named full-design coefficient,
    /// returning only the fit state.
    pub fn fit_lrt_formula_name(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_formula_with_results_name(
            counts,
            full_formula,
            reduced_formula,
            coefficient_name,
        )
        .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    pub fn fit_lrt_with_results_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    pub fn fit_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_contrast(counts, full_design, reduced_design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    pub fn fit_lrt_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_contrast_spec(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    pub fn fit_lrt_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_factor_level_contrast(counts, full_design, reduced_design, contrast)
    }

    /// Run the top-level LRT workflow for a DESeq2 `results(contrast=...)` request.
    pub fn fit_lrt_with_results_contrast_request<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if sample_levels.is_none() {
            if let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_lrt_with_results_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    factor_contrast,
                );
            }
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
                self.fit_lrt_with_results_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_lrt_with_results_contrast_spec(
                    counts,
                    full_design,
                    reduced_design,
                    &contrast_spec,
                )
            }
        }
    }

    /// Run the top-level LRT workflow for a DESeq2 `results(contrast=...)`
    /// request using formula model-frame metadata.
    pub fn fit_lrt_with_results_contrast_request_from_model_frame(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_lrt_with_results_factor_level_contrast(
                counts,
                full_design,
                reduced_design,
                factor_contrast,
            );
        }
        builder.fit_lrt_with_results_contrast_request::<String>(
            counts,
            full_design,
            reduced_design,
            contrast,
            None,
        )
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow for a DESeq2
    /// `results(contrast=...)` request.
    pub fn fit_lrt_formula_with_results_contrast_request(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        contrast: &ResultsContrast,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_contrast_request::<String>(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                contrast,
                None,
            )
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    pub fn fit_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_contrast_spec(counts, full_design, reduced_design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    pub fn fit_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_factor_level_contrast(
            counts,
            full_design,
            reduced_design,
            contrast,
        )
        .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let coefficient = default_results_coefficient(full_design)?;
        self.fit_lrt_glm_mu_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            coefficient,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level LRT replacement-refit workflow and report a named full-design coefficient.
    pub fn fit_lrt_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let coefficient = resolve_coefficient_index(full_design, coefficient_name)?;
        self.fit_lrt_glm_mu_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            coefficient,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level named LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level factor-level LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level LRT replacement-refit workflow for a DESeq2 `results(contrast=...)` request.
    pub fn fit_lrt_with_results_contrast_request_with_cooks_replacement<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        if sample_levels.is_none() {
            if let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)? {
                return self.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    factor_contrast,
                    replacement_options,
                );
            }
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
                self.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.fit_lrt_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    &contrast_spec,
                    replacement_options,
                )
            }
        }
    }

    /// Run the top-level LRT replacement-refit workflow for a DESeq2
    /// `results(contrast=...)` request using formula model-frame metadata.
    pub fn fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) = factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                counts,
                full_design,
                reduced_design,
                factor_contrast,
                replacement_options,
            );
        }
        builder.fit_lrt_with_results_contrast_request_with_cooks_replacement::<String>(
            counts,
            full_design,
            reduced_design,
            contrast,
            None,
            replacement_options,
        )
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT replacement-refit
    /// workflow for a DESeq2 `results(contrast=...)` request.
    pub fn fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        contrast: &ResultsContrast,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let full_formula_design = self.expanded_formula_design_with_offsets(full_formula)?;
        let reduced_design = self.standard_design_from_formula_without_offsets(reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_contrast_request_with_cooks_replacement::<String>(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                contrast,
                None,
                replacement_options,
            )
    }

    fn reduced_design_for_top_level_lrt(&self) -> Result<&DesignMatrix, DeseqError> {
        self.reduced_design
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "top-level LRT fit without a reduced design".to_string(),
            })
    }
}
