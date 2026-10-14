//! Rust core for `rsdeseq2`.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::core::*` paths remain source-compatible.

mod implementation;

pub use implementation::*;
