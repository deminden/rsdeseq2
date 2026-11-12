//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

use super::execution::ExecutionContext;
use super::fitted_results::{
    apply_contrast_all_zero_to_lrt_results, apply_contrast_all_zero_to_wald_contrast,
    factor_level_result_metadata,
};

// Count matrix storage, row access, and count summaries.
include!("implementation/counts.rs");

// Inspectable fit state, builder options, and transform methods on fit state.
include!("implementation/fit_state.rs");

// Private pipeline input/output bundles and public transform/replacement outputs.
include!("implementation/outputs.rs");

// Builder configuration setters and accessors.
include!("implementation/builder/config.rs");

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
include!("implementation/builder/wald_lrt.rs");

// Cook's replacement-refit builder workflows.
include!("implementation/builder/cooks_replacement.rs");

// Native Wald/LRT fit-state attachment helpers.
include!("implementation/builder/native_attachments.rs");

// Fixed-dispersion Wald builder workflows.
include!("implementation/builder/fixed_wald.rs");

// Fixed-dispersion LRT builder workflows.
include!("implementation/builder/fixed_lrt.rs");

// Shared normalization and fixed-dispersion pipeline components.
include!("implementation/builder/normalization_components.rs");

// Public top-level builder workflows.
#[path = "implementation/builder/workflows/mod.rs"]
mod workflows;

// Shared validation, compaction, expansion, and metadata helpers.
include!("implementation/helpers.rs");

// Unit tests for the core module implementation.
include!("implementation/tests.rs");
