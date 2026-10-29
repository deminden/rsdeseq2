//! Top-level likelihood-ratio Cook's replacement workflows.

use super::super::{
    ContrastSpec, CooksReplacementLrtOutput, CooksReplacementOptions, CountMatrix, DeseqBuilder,
    DeseqError, DesignMatrix, FactorLevelContrast, FormulaModelFrame, ResultsContrast,
    default_results_coefficient, factor_level_contrast_from_model_frame,
    factor_level_contrast_from_sample_levels, resolve_coefficient_index,
};

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results_with_cooks_replacement(
                    counts,
                    full_formula,
                    reduced_formula,
                    replacement_options,
                )
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results_name_with_cooks_replacement(
                    counts,
                    full_formula,
                    reduced_formula,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

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

    /// Run the currently implemented top-level LRT workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    replacement_options,
                )
            });
        }

        let coefficient = default_results_coefficient(full_design)?;
        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_name_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    coefficient_name,
                    replacement_options,
                )
            });
        }

        let coefficient = resolve_coefficient_index(full_design, coefficient_name)?;
        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_contrast_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
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
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)?
        {
            return self
                .with_automatic_dispersion_means()
                .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    factor_contrast,
                    replacement_options,
                );
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
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                        counts,
                        full_design,
                        reduced_design,
                        contrast,
                        replacement_options,
                    )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_spec_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    model_frame,
                    replacement_options,
                )
            });
        }

        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) =
            factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder
                .with_automatic_dispersion_means()
                .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
                    counts,
                    full_formula,
                    reduced_formula,
                    contrast,
                    replacement_options,
                )
            });
        }

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
}
