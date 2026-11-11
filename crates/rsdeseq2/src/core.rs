//! Rust core for `rsdeseq2`.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::core::*` paths remain source-compatible.

mod execution;
mod fitted_results;
mod implementation;
mod replacement_fit;

pub use implementation::*;
