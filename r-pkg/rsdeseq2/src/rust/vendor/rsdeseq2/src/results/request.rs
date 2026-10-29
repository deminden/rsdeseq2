//! Validated result policy, independent of fit execution and presentation.
use super::{
    DeseqResults, apply_cooks_cutoff, apply_cooks_cutoff_with_low_count_heuristic,
    resolve_cooks_cutoff,
};
use crate::contrasts::ResultsContrast;
use crate::core::{CountMatrix, DeseqFit};
use crate::errors::DeseqError;
use crate::glm::{WaldTestOptions, wald::validate_wald_options};
use crate::independent_filtering::{
    IndependentFilteringOptions, apply_independent_filtering, validate_alpha, validate_theta,
};
use crate::options::CooksCutoff;

/// A contrast and validated statistical/presentation options for a stored fit.
///
/// Construction checks scalar options; evaluation checks model-dependent state.
/// Each evaluation recomputes its own Cook's masking and filtering/BH decisions.
#[derive(Clone, Debug)]
pub struct ResultRequest {
    contrast: ResultsContrast,
    wald: WaldTestOptions,
    cooks: CooksCutoff,
    filtering: IndependentFilteringOptions,
    low_count_heuristic: bool,
}

impl ResultRequest {
    /// Validate explicit options without changing or refitting a stored model.
    pub fn new(
        contrast: ResultsContrast,
        wald: WaldTestOptions,
        cooks: CooksCutoff,
        filtering: IndependentFilteringOptions,
    ) -> Result<Self, DeseqError> {
        validate_wald_options(&wald)?;
        validate_alpha(filtering.alpha)?;
        if let Some(theta) = &filtering.theta {
            validate_theta(theta)?;
        }
        if let CooksCutoff::Threshold(value) = cooks
            && (!value.is_finite() || value < 0.0)
        {
            return Err(DeseqError::InvalidOptions {
                reason: "Cook's threshold must be finite and non-negative".into(),
            });
        }
        Ok(Self {
            contrast,
            wald,
            cooks,
            filtering,
            low_count_heuristic: false,
        })
    }

    /// Enable the two-group heuristic after the caller establishes its design context.
    pub fn low_count_heuristic(mut self, enabled: bool) -> Self {
        self.low_count_heuristic = enabled;
        self
    }
}

impl DeseqFit {
    /// Evaluate a validated request using borrowed stored inference arrays.
    ///
    /// This method serves ordinary stored Wald/LRT fits. Replacement outputs
    /// retain their separate merge/context policy. LRT requests must keep default
    /// Wald options: their p-values represent the stored full-vs-reduced test.
    pub fn results_for_request(
        &self,
        counts: &CountMatrix,
        request: &ResultRequest,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        let design = self
            .design
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "stored design is required for a result request".into(),
            })?;
        let mut results = match (self.wald.is_some(), self.lrt.is_some()) {
            (true, false) => self.wald_results_contrast_request_with_options(
                counts,
                &request.contrast,
                gene_names,
                &request.wald,
            )?,
            (false, true) if request.wald == WaldTestOptions::default() => {
                self.lrt_results_contrast_request(counts, &request.contrast, gene_names)?
            }
            (false, true) => {
                return Err(DeseqError::InvalidOptions {
                    reason: "Wald threshold/t options cannot change a stored LRT test".into(),
                });
            }
            _ => {
                return Err(DeseqError::InvalidOptions {
                    reason: "result requests require exactly one stored Wald or LRT test".into(),
                });
            }
        };
        let cutoff =
            resolve_cooks_cutoff(request.cooks, design.n_samples(), design.n_coefficients())?;
        if cutoff.is_some() && self.max_cooks.is_none() {
            return Err(DeseqError::InvalidOptions {
                reason: "Cook's filtering requires retained Cook's context".into(),
            });
        }
        if request.low_count_heuristic && cutoff.is_some() {
            let cooks = self
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "low-count filtering requires retained Cook's distances".into(),
                })?;
            apply_cooks_cutoff_with_low_count_heuristic(&mut results, cutoff, counts, cooks)?;
        } else {
            apply_cooks_cutoff(&mut results, cutoff)?;
        }
        apply_independent_filtering(&mut results, &request.filtering)?;
        Ok(results)
    }
}
