//! Design matrix construction helpers.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::design::*` paths remain source-compatible.

mod formula_origin;
mod implementation;
mod origin;

pub use implementation::*;
pub use origin::DesignOrigin;
