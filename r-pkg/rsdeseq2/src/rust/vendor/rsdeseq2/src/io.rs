//! TSV readers, writers, and label-alignment helpers.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::io::*` paths remain source-compatible.

mod implementation;

pub use implementation::*;
