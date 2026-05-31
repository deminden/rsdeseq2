//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Result row, table, metadata, and beta-prior workflow types.
include!("implementation/types.rs");

// Core Wald builders and caller-supplied expanded beta-prior workflows.
include!("implementation/beta_prior_core.rs");

// One-factor expanded beta-prior result workflows.
include!("implementation/beta_prior_factor.rs");

// Additive expanded beta-prior result workflows.
include!("implementation/beta_prior_additive.rs");

// Formula-driven expanded beta-prior result workflows.
include!("implementation/beta_prior_formula.rs");

// Cook's replacement merge helpers for beta-prior result workflows.
include!("implementation/cooks_replacement.rs");

// Wald/LRT result builders and Cook filtering.
include!("implementation/wald_lrt.rs");

// Validation, column assembly, and result metadata descriptions.
include!("implementation/metadata_validation.rs");
