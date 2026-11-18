//! Top-level Wald and selected-test Cook's replacement workflows.

use super::super::*;

impl DeseqBuilder {
    /// Run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            });
        }
        self.run(&Analysis::new(counts, design).cooks_replacement(replacement_options))?
            .into_wald_replacement()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_formula_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            });
        }
        self.run(&Analysis::formula(counts, formula).cooks_replacement(replacement_options))?
            .into_wald_replacement()
    }

    /// Run the top-level named Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            });
        }
        self.run(
            &Analysis::new(counts, design)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
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
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            });
        }
        self.run(
            &Analysis::formula(counts, formula)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the top-level workflow with limited Cook's replacement refit.
    pub fn fit_formula_with_test_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        formula: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
    }

    /// Run the top-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast".to_string(),
            });
        }
        self.run(
            &Analysis::new(counts, design)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the top-level numeric-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
    }

    /// Run the top-level named Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast specification"
                    .to_string(),
            });
        }
        self.run(
            &Analysis::new(counts, design)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the top-level named-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
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
        let test = self.top_level_test()?;
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                })
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: None,
                })
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                })
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
    }

    /// Run the top-level factor-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        if self.test == TestType::Lrt {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald factor-level contrast"
                    .to_string(),
            });
        }
        self.run(
            &Analysis::new(counts, design)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the top-level factor-level contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_replacement()
    }
}
