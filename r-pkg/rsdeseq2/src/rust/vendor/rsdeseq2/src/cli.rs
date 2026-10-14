//! Command-line interface for the `rsdeseq2` binary.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::cli::*` paths remain source-compatible.

mod implementation;

pub use implementation::*;
