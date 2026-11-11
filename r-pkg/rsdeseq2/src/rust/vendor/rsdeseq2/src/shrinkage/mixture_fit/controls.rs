//! Numerical policies for nondefault sequential-quadratic mixture fitting.
use super::{DeseqError, invalid};
#[derive(Clone, Copy, Debug)]
pub struct MixtureSqpOptions {
    pub active_tolerance: f64,
    /// None retains the established native direction policy.
    pub search_direction_threshold: Option<f64>,
    pub sufficient_decrease: f64,
    pub step_reduction: f64,
    pub minimum_step: f64,
    pub ridge_multiplier: f64,
    /// None (or zero) retains full-matrix fitting. Positive values request a
    /// deterministic truncated SVD, retaining at least two singular directions.
    pub svd_tolerance: Option<f64>,
    /// Resolve the active-set budget as min(20, components + 1).
    pub automatic_active_iterations: bool,
}
impl Default for MixtureSqpOptions {
    fn default() -> Self {
        Self {
            active_tolerance: 1e-10,
            search_direction_threshold: None,
            sufficient_decrease: 0.01,
            step_reduction: 0.75,
            minimum_step: 1e-8,
            ridge_multiplier: 10.0,
            svd_tolerance: None,
            automatic_active_iterations: false,
        }
    }
}
impl MixtureSqpOptions {
    pub(super) fn validate(self) -> Result<(), DeseqError> {
        if !self.active_tolerance.is_finite()
            || self.active_tolerance < 0.0
            || self
                .search_direction_threshold
                .is_some_and(|x| !x.is_finite() || x < 0.0)
            || !self.sufficient_decrease.is_finite()
            || self.sufficient_decrease <= 0.0
            || self.sufficient_decrease >= 1.0
            || !self.step_reduction.is_finite()
            || self.step_reduction <= 0.0
            || self.step_reduction >= 1.0
            || !self.minimum_step.is_finite()
            || self.minimum_step <= 0.0
            || self.minimum_step > 1.0
            || !self.ridge_multiplier.is_finite()
            || self.ridge_multiplier <= 0.0
            || self
                .svd_tolerance
                .is_some_and(|x| !x.is_finite() || x < 0.0)
        {
            return Err(invalid(
                "invalid SQP active-set, line-search, ridge or SVD controls",
            ));
        }
        Ok(())
    }
}
