//! Top-level Wald contrast workflows.

use super::super::*;

impl DeseqBuilder {
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the currently implemented top-level Wald workflow for a numeric contrast.
    pub fn fit_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_fit_only()
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
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
        let test = self.top_level_test()?;
        let levels = owned_levels(sample_levels);
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: levels.as_deref(),
                    model_frame: None,
                }),
        )?
        .into_results_pair()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run a DESeq2 `results(contrast=...)` request.
    pub fn fit_formula_with_results_contrast_request(
        &self,
        counts: &CountMatrix,
        formula: &str,
        contrast: &ResultsContrast,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: None,
                }),
        )?
        .into_results_pair()
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
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Results {
                    contrast,
                    sample_levels: None,
                    model_frame: Some(model_frame),
                }),
        )?
        .into_results_pair()
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_fit_only()
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_fit_only()
    }
}
