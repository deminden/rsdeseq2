//! Deferred R-style controls for the compiled negative-binomial refinement step.
use super::{
    negative_binomial::Row,
    prior::invalid,
    regression::{RegressionMethod, RegressionOptimizer, RegressionPosteriorOptions},
};
use crate::DeseqError;

#[derive(Clone, Debug, Default)]
pub enum NegativeBinomialRefinementMethod {
    #[default]
    Bfgs,
    Lbfgsb,
    /// An unrecognized method is irrelevant when the prefit is already valid.
    Unsupported(String),
}
#[derive(Clone, Debug)]
pub struct NegativeBinomialRefinement {
    pub method: NegativeBinomialRefinementMethod,
    /// Malformed bounds are rejected only if bounded refinement runs.
    /// Missing vector entries and non-finite endpoints denote absent bounds.
    /// R uses the first two values and ignores additional entries.
    pub bounds: Option<Vec<f64>>,
}
impl Default for NegativeBinomialRefinement {
    fn default() -> Self {
        Self {
            method: Default::default(),
            bounds: Some(vec![f64::NEG_INFINITY, f64::INFINITY]),
        }
    }
}
pub(super) fn optimize(
    row: &Row<'_>,
    start: Vec<f64>,
    control: &NegativeBinomialRefinement,
) -> Result<rcompat_lbfgsb::OptimResult, DeseqError> {
    let (optimizer, bounds) = match &control.method {
        NegativeBinomialRefinementMethod::Bfgs => (
            RegressionOptimizer::Bfgs,
            (f64::NEG_INFINITY, f64::INFINITY),
        ),
        NegativeBinomialRefinementMethod::Lbfgsb => {
            let b = control
                .bounds
                .as_ref()
                .ok_or_else(|| invalid("bounded refinement requires numeric bounds"))?;
            // R's optimizer treats non-finite endpoints as absent bounds.
            let lower = b
                .first()
                .copied()
                .filter(|v| v.is_finite())
                .unwrap_or(f64::NEG_INFINITY);
            let upper = b
                .get(1)
                .copied()
                .filter(|v| v.is_finite())
                .unwrap_or(f64::INFINITY);
            (RegressionOptimizer::Lbfgsb, (lower, upper))
        }
        NegativeBinomialRefinementMethod::Unsupported(name) => {
            return Err(invalid(&format!(
                "unsupported refinement optimizer '{name}'"
            )));
        }
    };
    super::regression::optimize_nb_r(
        row,
        &RegressionPosteriorOptions {
            optimizer,
            bounds,
            method: RegressionMethod::NegativeBinomialR,
            ..Default::default()
        },
        start,
    )
}
