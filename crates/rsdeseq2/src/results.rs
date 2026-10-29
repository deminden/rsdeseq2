//! DESeq-style result table assembly.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::results::*` paths remain source-compatible.

mod implementation;
mod request;

pub use implementation::*;
pub use request::ResultRequest;

pub(crate) use implementation::{
    build_lrt_projection_results_view, build_wald_contrast_results_view,
};
