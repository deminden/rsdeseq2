//! Validated result policy, independent of fit execution and presentation.
use super::{
    CooksPolicy, DeseqResults, LowCountHeuristic, mask_cooks_policy_rows, resolve_cooks_cutoff,
    two_level_design_matches,
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
    low_count_heuristic: Option<bool>,
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
            low_count_heuristic: None,
        })
    }

    /// Override the two-group heuristic. By default, stored formula factors and
    /// a matching two-column design establish eligibility automatically.
    pub fn low_count_heuristic(mut self, enabled: bool) -> Self {
        self.low_count_heuristic = Some(enabled);
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
        let mut low_count = request.low_count_heuristic.unwrap_or(false);
        if cutoff.is_some()
            && request.low_count_heuristic.is_none()
            && let Some(frame) = &self.model_frame
        {
            for factor in &frame.factors {
                let levels = factor.levels.as_ref().map_or_else(
                    || {
                        factor
                            .sample_levels
                            .iter()
                            .collect::<std::collections::BTreeSet<_>>()
                            .len()
                    },
                    Vec::len,
                );
                if two_level_design_matches(design.matrix().into(), &factor.sample_levels, levels)?
                {
                    low_count = true;
                    break;
                }
            }
        }
        let heuristic = if low_count && cutoff.is_some() {
            let cooks = self
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "low-count filtering requires retained Cook's distances".into(),
                })?;
            Some(LowCountHeuristic::new(
                crate::matrix::CountView::from(counts),
                cooks.into(),
            )?)
        } else {
            None
        };
        if cutoff.is_some() {
            mask_cooks_policy_rows(&mut results, CooksPolicy::new(cutoff, heuristic)?)?;
        }
        apply_independent_filtering(&mut results, &request.filtering)?;
        Ok(results)
    }
}
