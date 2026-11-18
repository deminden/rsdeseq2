//! Top-level likelihood-ratio Cook's replacement workflows.

use super::super::*;

impl DeseqBuilder {
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
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run the currently implemented top-level LRT workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Name(coefficient_name))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                })
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                })
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
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
        self.run(
            &Analysis::formula(counts, full_formula)
                .lrt_formula(reduced_formula)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: None,
                })
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }
}
