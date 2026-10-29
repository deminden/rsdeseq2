//! Top-level Wald contrast workflows.

use super::super::{
    ContrastSpec, CountMatrix, DeseqBuilder, DeseqError, DeseqFit, DeseqResults, DesignMatrix,
    FactorLevelContrast, FormulaModelFrame, ResultsContrast, TestType,
    factor_level_contrast_from_model_frame, factor_level_contrast_from_sample_levels,
};

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_with_results_contrast(counts, design, contrast));
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast(counts, design, reduced_design, contrast)
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_contrast(counts, design, contrast));
        }

        self.fit_with_results_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_with_results_contrast_spec(counts, design, contrast));
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_contrast_spec(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_spec(counts, design, reduced_design, contrast)
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
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)?
        {
            return self.fit_with_results_factor_level_contrast(counts, design, factor_contrast);
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_formula_with_results_contrast_request(counts, formula, contrast)
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_with_results_contrast_request_from_model_frame(
                    counts,
                    design,
                    contrast,
                    model_frame,
                )
            });
        }

        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) =
            factor_level_contrast_from_model_frame(contrast, model_frame)?
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_contrast_spec(counts, design, contrast));
        }

        self.fit_with_results_contrast_spec(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_with_results_factor_level_contrast(counts, design, contrast));
        }

        match self.test {
            TestType::Wald => self
                .with_automatic_dispersion_means()
                .fit_wald_glm_mu_factor_level_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_factor_level_contrast(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_factor_level_contrast(counts, design, contrast));
        }

        self.fit_with_results_factor_level_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }
}
