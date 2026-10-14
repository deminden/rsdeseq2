//! Rust core for `rsdeseq2`.
//!
//! This crate currently implements early DESeq2-compatible stages:
//! count/design matrix validation, size-factor estimation, normalized counts,
//! gene/sample normalization factors for implemented fixed-dispersion paths,
//! base means, optionally weighted fixed-dispersion IRLS with scalar or
//! per-coefficient ridge values, normal-equation or QR updates, and bounded
//! optim fallback refits for routed rows,
//! observation-weight preprocessing, primitive Wald linear contrasts, Wald/LRT
//! statistics, Cook's diagnostics and low-count heuristic helper, primitive
//! outlier count replacement, limited GLM-mu native Wald/LRT replacement refit
//! metadata, BH adjusted p-values, independent filtering, and an inspectable
//! fit-state skeleton, plus implemented normTransform, VST helpers, and a
//! fit-state and builder-level rlog sample-effect ridge-GLM path.

pub(crate) mod all_zero;
pub mod bindings;
pub mod cli;
pub mod contrasts;
pub mod cooks;
pub mod core;
pub mod design;
pub mod diagnostics;
pub mod dispersion;
pub mod errors;
pub mod glm;
pub mod independent_filtering;
pub mod io;
pub mod math;
pub mod matrix;
pub mod multiple_testing;
pub mod normalization;
pub mod options;
pub mod prelude;
pub mod results;
pub mod transform;

pub use crate::errors::DeseqError;
