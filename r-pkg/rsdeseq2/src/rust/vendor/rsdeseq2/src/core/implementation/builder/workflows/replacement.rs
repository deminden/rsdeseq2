//! Top-level Wald and selected-test Cook's replacement workflows.

use super::super::{
    ContrastSpec, CooksReplacementOptions, CooksReplacementTestOutput, CooksReplacementWaldOutput,
    CountMatrix, DeseqBuilder, DeseqError, DesignMatrix, FactorLevelContrast, FormulaModelFrame,
    ResultsContrast, TestType, default_results_coefficient, factor_level_contrast_from_model_frame,
    factor_level_contrast_from_sample_levels, resolve_coefficient_index,
};

impl DeseqBuilder {
    /// Run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_with_cooks_replacement(counts, design, replacement_options)
            });
        }

        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.with_automatic_dispersion_means()
                    .fit_wald_glm_mu_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_results_with_cooks_replacement(
                    counts,
                    formula,
                    replacement_options,
                )
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_name_with_cooks_replacement(
                    counts,
                    design,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

        match self.test {
            TestType::Wald => {
                let coefficient = resolve_coefficient_index(design, coefficient_name)?;
                self.with_automatic_dispersion_means()
                    .fit_wald_glm_mu_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_results_name_with_cooks_replacement(
                    counts,
                    formula,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_with_cooks_replacement(
                    counts,
                    design,
                    replacement_options,
                )
            });
        }

        match self.test {
            TestType::Wald => self
                .fit_with_results_with_cooks_replacement(counts, design, replacement_options)
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_test_results_with_cooks_replacement(
                    counts,
                    formula,
                    replacement_options,
                )
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_name_with_cooks_replacement(
                    counts,
                    design,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_name_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_test_results_name_with_cooks_replacement(
                    counts,
                    formula,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_test_results_name_with_cooks_replacement(
                counts,
                &formula_design.design.standard_design,
                coefficient_name,
                replacement_options,
            )
    }

    /// Run the top-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_contrast_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_with_cooks_replacement(
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

    /// Run the top-level named Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_spec_with_cooks_replacement(
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
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)?
        {
            return self.fit_with_test_results_factor_level_contrast_with_cooks_replacement(
                counts,
                design,
                factor_contrast,
                replacement_options,
            );
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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_request_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_test_results_contrast_request_with_cooks_replacement(
                    counts,
                    formula,
                    contrast,
                    replacement_options,
                )
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    model_frame,
                    replacement_options,
                )
            });
        }

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
                    .with_automatic_dispersion_means()
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

    /// Run the top-level factor-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_test_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
            });
        }

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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
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
}
