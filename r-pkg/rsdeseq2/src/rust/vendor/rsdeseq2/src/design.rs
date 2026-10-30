//! Design matrix construction helpers.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::design::*` paths remain source-compatible.

mod formula_origin;
mod implementation;
mod model_matrix;
mod origin;

pub use implementation::*;
pub use model_matrix::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
pub use origin::DesignOrigin;
