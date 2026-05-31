//! DESeq-style result table assembly.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::results::*` paths remain source-compatible.

mod implementation;

pub use implementation::*;
