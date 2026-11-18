//! Top-level likelihood-ratio fit workflows.

use super::super::*;

impl DeseqBuilder {
    /// Run the currently implemented top-level LRT workflow with a reduced design.
    ///
    /// This follows the implemented GLM-mu native LRT path and reports the last
    /// full-design coefficient by default.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design))"
    )]
    pub fn fit_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(&Analysis::new(counts, full_design).lrt(reduced_design))?
            .into_fit_only()
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow, returning only
    /// the fit state.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula))"
    )]
    pub fn fit_lrt_formula(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula))?
            .into_fit_only()
    }

    /// Run the currently implemented top-level LRT workflow and return result rows.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design))"
    )]
    pub fn fit_lrt_with_results(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(&Analysis::new(counts, full_design).lrt(reduced_design))?
            .into_results_pair()
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow.
    ///
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula))"
    )]
    pub fn fit_lrt_formula_with_results(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula))?
            .into_results_pair()
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Name(coefficient_name)))"
    )]
    pub fn fit_lrt_with_results_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_results_pair()
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata, run LRT, and report a named full-design
    /// coefficient effect.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula).result(ResultColumn::Name(coefficient_name)))"
    )]
    pub fn fit_lrt_formula_with_results_name(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_results_pair()
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Name(coefficient_name)))"
    )]
    pub fn fit_lrt_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_fit_only()
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and report a named full-design coefficient,
    /// returning only the fit state.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula).result(ResultColumn::Name(coefficient_name)))"
    )]
    pub fn fit_lrt_formula_name(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_fit_only()
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Numeric(contrast)))"
    )]
    pub fn fit_lrt_with_results_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Numeric(contrast)))"
    )]
    pub fn fit_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_fit_only()
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Spec(contrast)))"
    )]
    pub fn fit_lrt_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::FactorLevel(contrast)))"
    )]
    pub fn fit_lrt_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the top-level LRT workflow for a DESeq2 `results(contrast=...)` request.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Results { contrast, sample_levels, model_frame: None }))"
    )]
    pub fn fit_lrt_with_results_contrast_request<S: AsRef<str>>(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        sample_levels: Option<&[S]>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                }),
        )?
        .into_results_pair()
    }

    /// Run the top-level LRT workflow for a DESeq2 `results(contrast=...)`
    /// request using formula model-frame metadata.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Results { contrast, sample_levels: None, model_frame: Some(model_frame) }))"
    )]
    pub fn fit_lrt_with_results_contrast_request_from_model_frame(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ResultsContrast,
        model_frame: &FormulaModelFrame,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                }),
        )?
        .into_results_pair()
    }

    /// Build supported full and reduced formula designs from stored
    /// model-frame metadata and run the top-level LRT workflow for a DESeq2
    /// `results(contrast=...)` request.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, full_formula).lrt_formula(reduced_formula).result(ResultColumn::Results { contrast, sample_levels: None, model_frame: None }))"
    )]
    pub fn fit_lrt_formula_with_results_contrast_request(
        &self,
        counts: &CountMatrix,
        full_formula: &str,
        reduced_formula: &str,
        contrast: &ResultsContrast,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: None,
                }),
        )?
        .into_results_pair()
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::Spec(contrast)))"
    )]
    pub fn fit_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_fit_only()
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, full_design).lrt(reduced_design).result(ResultColumn::FactorLevel(contrast)))"
    )]
    pub fn fit_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_fit_only()
    }

    pub(super) fn reduced_design_for_top_level_lrt(&self) -> Result<&DesignMatrix, DeseqError> {
        self.reduced_design
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "top-level LRT fit without a reduced design".to_string(),
            })
    }
}
