//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

use super::execution::ExecutionContext;
use super::fitted_results::{
    apply_contrast_all_zero_to_lrt_results, apply_contrast_all_zero_to_wald_contrast,
    factor_level_result_metadata,
};
use crate::numeric::checked;

// Count matrix storage, row access, and count summaries.
#[path = "implementation/counts.rs"]
mod counts;
pub use counts::*;

// Inspectable fit state, builder options, and transform methods on fit state.
#[path = "implementation/fit_state.rs"]
mod fit_state;
pub use fit_state::*;

// Private pipeline input/output bundles and public transform/replacement outputs.
#[path = "implementation/outputs.rs"]
mod outputs;
pub use outputs::*;

// Request and output types for `DeseqBuilder::run`.
#[path = "implementation/analysis.rs"]
mod analysis;
pub use analysis::*;

// Builder configuration setters and accessors.
#[path = "implementation/builder/config.rs"]
mod config;
use config::*;

// Builder stage boundaries retain the existing public forwarding methods.
#[path = "implementation/builder/stages/gene_wise.rs"]
mod builder_gene_wise;
#[path = "implementation/builder/stages/map.rs"]
mod builder_map;
#[path = "implementation/builder/stages/normalization.rs"]
mod builder_normalization;
#[path = "implementation/builder/stages/rlog.rs"]
mod builder_rlog;
#[path = "implementation/builder/stages/trends.rs"]
mod builder_trends;
#[path = "implementation/builder/stages/vst.rs"]
mod builder_vst;

// Builder Wald and LRT orchestration before replacement refits.
#[path = "implementation/builder/wald_lrt.rs"]
mod wald_lrt;

// Cook's replacement-refit builder workflows.
#[path = "implementation/builder/cooks_replacement.rs"]
mod cooks_replacement;

// Native Wald/LRT fit-state attachment helpers.
#[path = "implementation/builder/native_attachments.rs"]
mod native_attachments;

// Fixed-dispersion Wald builder workflows.
#[path = "implementation/builder/fixed_wald.rs"]
mod fixed_wald;

// Fixed-dispersion LRT builder workflows.
#[path = "implementation/builder/fixed_lrt.rs"]
mod fixed_lrt;

// Shared normalization and fixed-dispersion pipeline components.
#[path = "implementation/builder/normalization_components.rs"]
mod normalization_components;

// DESeq2's beta prior (`DESeq(betaPrior = TRUE)`) for `run`.
#[path = "implementation/builder/beta_prior.rs"]
mod beta_prior;

// Explicit-path `fit*` methods, each a `run` request.
#[path = "implementation/builder/leaf_forwarders.rs"]
mod leaf_forwarders;

// Public top-level builder workflows.
#[path = "implementation/builder/workflows/mod.rs"]
mod workflows;

// Shared validation, compaction, expansion, and metadata helpers.
#[path = "implementation/helpers.rs"]
mod helpers;
use helpers::*;

// Unit tests for the core module implementation.
#[cfg(test)]
#[path = "implementation/tests.rs"]
mod tests;
