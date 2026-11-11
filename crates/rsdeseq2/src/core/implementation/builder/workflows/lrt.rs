//! Top-level likelihood-ratio fit workflows.

use super::super::{
    ContrastSpec, CountMatrix, DeseqBuilder, DeseqError, DeseqFit, DeseqResults, DesignMatrix,
    FactorLevelContrast, FormulaModelFrame, ResultsContrast, default_results_coefficient,
    factor_level_contrast_from_model_frame, factor_level_contrast_from_sample_levels,
    resolve_coefficient_index,
};

impl DeseqBuilder {
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_lrt(counts, full_design, reduced_design));
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_with_results(counts, full_design, reduced_design)
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_lrt_formula(counts, full_formula, reduced_formula));
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_lrt_with_results(counts, full_design, reduced_design));
        }

        let coefficient = default_results_coefficient(full_design)?;
        self.with_automatic_dispersion_means().fit_lrt_glm_mu(
            counts,
            full_design,
            reduced_design,
            coefficient,
        )
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results(counts, full_formula, reduced_formula)
            });
        }

        let full_formula_design = self.formula_design_for_fitting(full_formula)?;
        let reduced_design = self.reduced_formula_design(&full_formula_design, reduced_formula)?;
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_name(
                    counts,
                    full_design,
                    reduced_design,
                    coefficient_name,
                )
            });
        }

        let coefficient = resolve_coefficient_index(full_design, coefficient_name)?;
        self.with_automatic_dispersion_means().fit_lrt_glm_mu(
            counts,
            full_design,
            reduced_design,
            coefficient,
        )
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results_name(
                    counts,
                    full_formula,
                    reduced_formula,
                    coefficient_name,
                )
            });
        }

        let full_formula_design = self.formula_design_for_fitting(full_formula)?;
        let reduced_design = self.reduced_formula_design(&full_formula_design, reduced_formula)?;
        self.with_formula_offsets(counts, &full_formula_design)?
            .fit_lrt_with_results_name(
                counts,
                &full_formula_design.design.standard_design,
                &reduced_design,
                coefficient_name,
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_name(counts, full_design, reduced_design, coefficient_name)
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_with_results_name(counts, full_design, reduced_design, coefficient_name)
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_name(counts, full_formula, reduced_formula, coefficient_name)
            });
        }

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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast(counts, full_design, reduced_design, contrast)
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    pub fn fit_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_lrt_contrast(counts, full_design, reduced_design, contrast));
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_with_results_contrast(counts, full_design, reduced_design, contrast)
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast_spec(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                )
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_contrast_spec(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    pub fn fit_lrt_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                )
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_glm_mu_factor_level_contrast(counts, full_design, reduced_design, contrast)
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
        if sample_levels.is_none()
            && let Some(factor_contrast) = self.model_frame_factor_level_contrast(contrast)?
        {
            return self
                .with_automatic_dispersion_means()
                .fit_lrt_with_results_factor_level_contrast(
                    counts,
                    full_design,
                    reduced_design,
                    factor_contrast,
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
                    .fit_lrt_with_results_factor_level_contrast(
                        counts,
                        full_design,
                        reduced_design,
                        contrast,
                    )
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                let contrast_spec = contrast.as_contrast_spec();
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_contrast_spec(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_with_results_contrast_request_from_model_frame(
                    counts,
                    full_design,
                    reduced_design,
                    contrast,
                    model_frame,
                )
            });
        }

        let builder = self.clone().try_model_frame(model_frame.clone())?;
        if let Some(factor_contrast) =
            factor_level_contrast_from_model_frame(contrast, model_frame)?
        {
            return builder
                .with_automatic_dispersion_means()
                .fit_lrt_with_results_factor_level_contrast(
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_formula_with_results_contrast_request(
                    counts,
                    full_formula,
                    reduced_formula,
                    contrast,
                )
            });
        }

        let full_formula_design = self.formula_design_for_fitting(full_formula)?;
        let reduced_design = self.reduced_formula_design(&full_formula_design, reduced_formula)?;
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_contrast_spec(counts, full_design, reduced_design, contrast)
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_with_results_contrast_spec(counts, full_design, reduced_design, contrast)
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.fit_lrt_factor_level_contrast(counts, full_design, reduced_design, contrast)
            });
        }

        self.with_automatic_dispersion_means()
            .fit_lrt_with_results_factor_level_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
            )
            .map(|(fit, _results)| fit)
    }

    pub(super) fn reduced_design_for_top_level_lrt(&self) -> Result<&DesignMatrix, DeseqError> {
        self.reduced_design
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "top-level LRT fit without a reduced design".to_string(),
            })
    }
}
