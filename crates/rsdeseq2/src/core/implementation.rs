//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Count matrix storage, row access, and count summaries.
include!("implementation/counts.rs");

// Inspectable fit state, builder options, and transform methods on fit state.
include!("implementation/fit_state.rs");

// Private pipeline input/output bundles and public transform/replacement outputs.
include!("implementation/outputs.rs");

// Builder configuration setters and accessors.
include!("implementation/builder/config.rs");

// Builder normalization, dispersion, VST, and rlog stages.
include!("implementation/builder/dispersion_transforms.rs");

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
include!("implementation/builder/top_level.rs");

// Shared validation, compaction, expansion, and metadata helpers.
include!("implementation/helpers.rs");

// Unit tests for the core module implementation.
include!("implementation/tests.rs");
