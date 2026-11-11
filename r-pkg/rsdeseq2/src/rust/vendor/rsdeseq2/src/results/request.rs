//! Validated result policy, independent of fit execution and presentation.
use super::{
    CooksPolicy, DeseqResults, LowCountHeuristic, ReplacementRows, low_count_eligible_for_frame,
    mask_cooks_policy_rows, resolve_cooks_cutoff,
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

    /// Override the two-group heuristic. Automatic scope requires an original
    /// single-variable factor formula; numerical partition matching alone does
    /// not establish formula provenance.
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
        self.evaluate_result_request(counts, request, gene_names, None)
    }

    /// Evaluate a later request on a merged replacement fit with its retained plan.
    ///
    /// Pass the original counts and the plan returned by `into_merged_fit`.
    /// Replacement zero rows receive effect/SE/stat/p = 0/0/0/1 before the
    /// shared Cook's policy and this request's filtering/BH step. Original
    /// Cook's maxima remain diagnostics; refitted rows are not masked again.
    pub fn results_for_request_with_replacement(
        &self,
        counts: &CountMatrix,
        request: &ResultRequest,
        gene_names: Option<&[String]>,
        plan: &crate::cooks::CooksRefitPlan,
    ) -> Result<DeseqResults, DeseqError> {
        self.evaluate_result_request(counts, request, gene_names, Some(plan))
    }

    fn evaluate_result_request(
        &self,
        counts: &CountMatrix,
        request: &ResultRequest,
        gene_names: Option<&[String]>,
        plan: Option<&crate::cooks::CooksRefitPlan>,
    ) -> Result<DeseqResults, DeseqError> {
        let replacement = plan.map(ReplacementRows::from_plan).transpose()?;
        if let Some(rows) = replacement
            && rows.genes() != self.base_mean.len()
        {
            return Err(crate::errors::invalid_dimensions(
                "replacement result rows",
                self.base_mean.len(),
                rows.genes(),
            ));
        }
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
        if let Some(plan) = plan
            && plan.n_refit > 0
            && plan.requires_refit()
        {
            if plan.post_refit_max_cooks.len() != results.rows.len() {
                return Err(crate::errors::invalid_dimensions(
                    "post-refit maximum Cook's rows",
                    results.rows.len(),
                    plan.post_refit_max_cooks.len(),
                ));
            }
            for (row, &maximum) in results.rows.iter_mut().zip(&plan.post_refit_max_cooks) {
                row.max_cooks = maximum;
            }
        }
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
            low_count = low_count_eligible_for_frame(design, frame)?;
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
        let mut policy = CooksPolicy::new(cutoff, heuristic)?;
        if let Some(rows) = replacement {
            for (gene, row) in results.rows.iter_mut().enumerate() {
                if rows.is_zero(gene) {
                    ReplacementRows::set_zero_values(
                        &mut row.log2_fold_change,
                        &mut row.lfc_se,
                        &mut row.stat,
                        &mut row.pvalue,
                    );
                }
            }
            policy = policy.with_replacement(rows);
        }
        if cutoff.is_some() || replacement.is_some() {
            mask_cooks_policy_rows(&mut results, policy)?;
        }
        apply_independent_filtering(&mut results, &request.filtering)?;
        Ok(results)
    }
}
