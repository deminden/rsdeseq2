//! Stored-fit inference and contrast policy behind stable public forwarding APIs.
use super::implementation::validate_fit_counts_shape;
use super::{CountMatrix, DeseqFit};
use crate::contrasts::{
    ContrastSpec, FactorLevelContrast, ResultsContrast, contrast_all_zero_factor_levels,
    contrast_all_zero_numeric, factor_level_contrast_from_model_frame, resolve_contrast,
};
use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::{WaldContrastOutput, WaldTestOptions};
use crate::results::DeseqResults;

impl DeseqFit {
    /// Assemble `results(contrast=...)` rows from this already-fitted object.
    ///
    /// The stored test state chooses Wald or LRT behavior. Character triplet
    /// contrasts use retained formula/model-frame metadata when available.
    pub fn results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        match (self.wald.is_some(), self.lrt.is_some()) {
            (true, false) => self.wald_results_contrast_request(counts, contrast, gene_names),
            (false, true) => self.lrt_results_contrast_request(counts, contrast, gene_names),
            (true, true) => Err(DeseqError::InvalidOptions {
                reason:
                    "fitted-object contrast results require either Wald or LRT output, not both"
                        .to_string(),
            }),
            (false, false) => Err(DeseqError::InvalidOptions {
                reason: "fitted-object contrast results require stored Wald or LRT output"
                    .to_string(),
            }),
        }
    }

    /// Assemble Wald `results(contrast=...)` rows from this already-fitted object.
    ///
    /// Character triplet contrasts use stored formula/model-frame metadata for
    /// factor levels and references. List and numeric contrasts use the design
    /// coefficient names carried by the fitted GLM state.
    pub fn wald_results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        self.wald_results_contrast_request_with_options(
            counts,
            contrast,
            gene_names,
            &WaldTestOptions::default(),
        )
    }

    /// Assemble stored-fit contrast results with explicit Wald options.
    pub fn wald_results_contrast_request_with_options(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
        options: &WaldTestOptions,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "Wald results contrast request")?;
        let fit = self.inference_view()?;
        let numeric_contrast = match contrast {
            ResultsContrast::Character { .. } => {
                let model_frame =
                    self.model_frame
                        .as_ref()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "character results contrast requires stored formula model-frame metadata; use wald_results_factor_level_contrast when supplying explicit sample levels".to_string(),
                        })?;
                let factor_contrast = factor_level_contrast_from_model_frame(contrast, model_frame)?
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "character results contrast does not match stored formula model-frame metadata".to_string(),
                    })?;
                return self.wald_results_factor_level_contrast_with_options(
                    counts,
                    factor_contrast,
                    gene_names,
                    options,
                );
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                contrast.as_contrast_spec()
            }
        };
        let numeric = resolve_contrast(fit.model_matrix, &numeric_contrast)?;
        let mut output = crate::glm::wald::wald_test_contrast_view(&fit, &numeric, options)?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, fit.model_matrix, &numeric)?;
        apply_contrast_all_zero_to_wald_contrast(&mut output, &contrast_all_zero, &self.all_zero)?;
        let mut results = crate::results::build_wald_contrast_results_view(
            &self.base_mean,
            &fit,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        results.apply_wald_test_options(options);
        self.attach_fit_diagnostics_to_results(&mut results)?;
        results.set_resolved_contrast_metadata(
            numeric_contrast.result_name(),
            numeric_contrast.comparison(),
            &numeric,
        );
        Ok(results)
    }

    /// Assemble Wald `results(contrast=...)` rows for an explicit factor-level contrast.
    ///
    /// This is the fitted-object route for callers that already know the
    /// per-sample factor levels. Formula-built fits can instead use
    /// [`DeseqFit::wald_results_contrast_request`] with stored model-frame
    /// metadata.
    pub fn wald_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        contrast: FactorLevelContrast<'_>,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        self.wald_results_factor_level_contrast_with_options(
            counts,
            contrast,
            gene_names,
            &WaldTestOptions::default(),
        )
    }

    /// Assemble stored-fit contrast results with explicit Wald options.
    pub fn wald_results_factor_level_contrast_with_options(
        &self,
        counts: &CountMatrix,
        contrast: FactorLevelContrast<'_>,
        gene_names: Option<&[String]>,
        options: &WaldTestOptions,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "Wald factor-level results contrast")?;
        let fit = self.inference_view()?;
        let contrast_spec = factor_level_contrast_spec(contrast);
        let numeric = resolve_contrast(fit.model_matrix, &contrast_spec)?;
        let contrast_all_zero = contrast_all_zero_factor_levels(
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let mut output = crate::glm::wald::wald_test_contrast_view(&fit, &numeric, options)?;
        apply_contrast_all_zero_to_wald_contrast(&mut output, &contrast_all_zero, &self.all_zero)?;
        let mut results = crate::results::build_wald_contrast_results_view(
            &self.base_mean,
            &fit,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        results.apply_wald_test_options(options);
        self.attach_fit_diagnostics_to_results(&mut results)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric);
        Ok(results)
    }

    /// Assemble LRT `results(contrast=...)` rows from this already-fitted object.
    ///
    /// The LRT statistic and p-value come from the stored full-vs-reduced test;
    /// the requested contrast only changes the reported effect-size columns.
    pub fn lrt_results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "LRT results contrast request")?;
        let fit = self.inference_view()?;
        let lrt = self
            .lrt
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "stored LRT output is required for LRT contrast results".to_string(),
            })?;
        let numeric_contrast = match contrast {
            ResultsContrast::Character { .. } => {
                let model_frame =
                    self.model_frame
                        .as_ref()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "character results contrast requires stored formula model-frame metadata; use lrt_results_factor_level_contrast when supplying explicit sample levels".to_string(),
                        })?;
                let factor_contrast = factor_level_contrast_from_model_frame(contrast, model_frame)?
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "character results contrast does not match stored formula model-frame metadata".to_string(),
                    })?;
                return self.lrt_results_factor_level_contrast(counts, factor_contrast, gene_names);
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                contrast.as_contrast_spec()
            }
        };
        let numeric = resolve_contrast(fit.model_matrix, &numeric_contrast)?;
        let output = fit.estimates()?.project_contrast(&numeric)?;
        let mut results = crate::results::build_lrt_projection_results_view(
            &self.base_mean,
            &fit,
            lrt,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        self.attach_fit_diagnostics_to_results(&mut results)?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, fit.model_matrix, &numeric)?;
        apply_contrast_all_zero_to_lrt_results(&mut results, &contrast_all_zero, &self.all_zero)?;
        results.set_resolved_contrast_metadata(
            numeric_contrast.result_name(),
            numeric_contrast.comparison(),
            &numeric,
        );
        Ok(results)
    }

    /// Assemble LRT `results(contrast=...)` rows for an explicit factor-level contrast.
    ///
    /// The stored full-vs-reduced LRT supplies the statistic and p-value; the
    /// explicit factor contrast supplies reported effect-size columns and
    /// character-contrast all-zero handling.
    pub fn lrt_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        contrast: FactorLevelContrast<'_>,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "LRT factor-level results contrast")?;
        let fit = self.inference_view()?;
        let lrt = self
            .lrt
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "stored LRT output is required for LRT contrast results".to_string(),
            })?;
        let contrast_spec = factor_level_contrast_spec(contrast);
        let numeric = resolve_contrast(fit.model_matrix, &contrast_spec)?;
        let contrast_all_zero = contrast_all_zero_factor_levels(
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let output = fit.estimates()?.project_contrast(&numeric)?;
        let mut results = crate::results::build_lrt_projection_results_view(
            &self.base_mean,
            &fit,
            lrt,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        self.attach_fit_diagnostics_to_results(&mut results)?;
        apply_contrast_all_zero_to_lrt_results(&mut results, &contrast_all_zero, &self.all_zero)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric);
        Ok(results)
    }

    fn inference_view(&self) -> Result<crate::glm::inference::GlmInferenceView<'_>, DeseqError> {
        let missing = |field: &str| DeseqError::InvalidOptions {
            reason: format!("stored {field} is required for fitted-object results"),
        };
        Ok(crate::glm::inference::GlmInferenceView {
            beta: self.beta.as_ref().ok_or_else(|| missing("beta matrix"))?,
            beta_se: self
                .beta_se
                .as_ref()
                .ok_or_else(|| missing("beta standard-error matrix"))?,
            beta_covariance: self.beta_covariance.as_ref(),
            beta_converged: self
                .beta_converged
                .as_deref()
                .ok_or_else(|| missing("beta convergence flags"))?,
            model_matrix: self.design.as_ref().ok_or_else(|| missing("design"))?,
        })
    }

    fn attach_fit_diagnostics_to_results(
        &self,
        results: &mut DeseqResults,
    ) -> Result<(), DeseqError> {
        if results.rows.len() != self.all_zero.len() {
            return Err(invalid_dimensions(
                "fitted-object result rows",
                self.all_zero.len(),
                results.rows.len(),
            ));
        }
        if let Some(max_cooks) = &self.max_cooks {
            if max_cooks.len() != results.rows.len() {
                return Err(invalid_dimensions(
                    "fitted-object maxCooks rows",
                    results.rows.len(),
                    max_cooks.len(),
                ));
            }
            for (row, max_cook) in results.rows.iter_mut().zip(max_cooks.iter().copied()) {
                row.max_cooks = max_cook;
            }
        }
        for (row, all_zero) in results.rows.iter_mut().zip(&self.all_zero) {
            if *all_zero {
                row.converged = None;
                row.max_cooks = None;
            }
        }
        Ok(())
    }
}

fn factor_level_contrast_spec(contrast: FactorLevelContrast<'_>) -> ContrastSpec {
    match contrast.reference {
        Some(reference) => ContrastSpec::factor_level_with_reference(
            contrast.factor,
            contrast.numerator,
            contrast.denominator,
            reference,
        ),
        None => {
            ContrastSpec::factor_level(contrast.factor, contrast.numerator, contrast.denominator)
        }
    }
}

pub(super) fn factor_level_result_metadata(contrast: FactorLevelContrast<'_>) -> (String, String) {
    let spec = factor_level_contrast_spec(contrast);
    (spec.result_name(), spec.comparison())
}

pub(super) fn apply_contrast_all_zero_to_wald_contrast(
    contrast: &mut WaldContrastOutput,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    crate::contrasts::ContrastZeroPolicy::new(contrast_all_zero, all_zero)?.apply_wald(
        &mut contrast.log2_fold_change,
        &mut contrast.wald.stat,
        &mut contrast.wald.pvalue,
    )
}

pub(super) fn apply_contrast_all_zero_to_lrt_results(
    results: &mut DeseqResults,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    crate::contrasts::ContrastZeroPolicy::new(contrast_all_zero, all_zero)?
        .apply_effects(results.rows.iter_mut().map(|row| &mut row.log2_fold_change))
}
