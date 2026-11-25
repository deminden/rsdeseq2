//! ashr-compatible policy of the mixture SQP solver.
//!
//! Not available in this revision: the previous implementation was removed
//! and is being replaced by one written from
//! `docs/provenance/specs/mixture_sqp_compatibility.md`.
use super::{DMatrix, DVector, ERROR_SCOPE, Error, Problem};
use crate::errors::scoped_invalid_options;
use crate::shrinkage::MixtureSqpOptions;

impl Problem {
    pub(super) fn ashr_hessian_into(&self, weights: &[f64], hessian: &mut DMatrix<f64>) {
        let _ = weights;
        hessian.fill(f64::NAN);
    }

    pub(super) fn ashr_expectation_step(&self, weights: &[f64]) -> Result<Vec<f64>, Error> {
        let _ = weights;
        Err(unavailable())
    }
}

pub(super) fn quadratic(
    h: &DMatrix<f64>,
    c: &DVector<f64>,
    start: &DVector<f64>,
    max_iterations: usize,
    options: MixtureSqpOptions,
) -> Result<DVector<f64>, Error> {
    let _ = (h, c, start, max_iterations, options);
    Err(unavailable())
}

fn unavailable() -> Error {
    scoped_invalid_options(ERROR_SCOPE, "ashr-compatible SQP policy is being rewritten")
}
