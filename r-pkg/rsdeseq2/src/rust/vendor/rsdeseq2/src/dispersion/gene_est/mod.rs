//! Gene-wise dispersion estimation primitives.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::dispersion::gene_est::*` paths remain source-compatible.

mod implementation;

pub use implementation::*;
