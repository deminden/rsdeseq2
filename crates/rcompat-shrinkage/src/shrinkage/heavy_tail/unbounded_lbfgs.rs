//! Unbounded L-BFGS for posterior modes.
//!
//! Not available in this revision: the previous implementation was removed
//! and is being replaced by one written from
//! `docs/provenance/specs/unbounded_lbfgs.md`.
use rcompat_lbfgsb::{OptimCounts, OptimResult};

/// Minimise with six corrections, 300 iterations and tolerances 1e-8.
pub(super) fn minimize(f: impl FnMut(&[f64]) -> (f64, Vec<f64>), start: Vec<f64>) -> OptimResult {
    minimize_with_tolerances(f, start, 1e-8, 1e-8)
}

/// Minimise with the given function and gradient tolerances.
pub(super) fn minimize_with_tolerances(
    f: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    start: Vec<f64>,
    function_tolerance: f64,
    gradient_tolerance: f64,
) -> OptimResult {
    let _ = (f, function_tolerance, gradient_tolerance);
    OptimResult {
        par: start,
        value: f64::NAN,
        counts: OptimCounts {
            function: 0,
            gradient: 0,
        },
        convergence: 52,
        message: "unbounded L-BFGS is being rewritten".to_string(),
    }
}
