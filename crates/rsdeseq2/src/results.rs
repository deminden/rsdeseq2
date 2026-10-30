//! DESeq-style result table assembly.
//!
//! Public items are re-exported from a split private implementation so
//! existing `rsdeseq2::results::*` paths remain source-compatible.

mod cooks_policy;
mod implementation;
mod replacement_policy;
mod request;

pub use cooks_policy::{
    CooksMask, CooksPolicy, LowCountHeuristic, low_count_eligible, two_level_design_matches,
};
pub use implementation::*;
pub use replacement_policy::{ReplacementRows, RowFlags};
pub use request::ResultRequest;

pub(crate) use implementation::{
    build_lrt_projection_results_view, build_wald_contrast_results_view, mask_cooks_policy_rows,
};
