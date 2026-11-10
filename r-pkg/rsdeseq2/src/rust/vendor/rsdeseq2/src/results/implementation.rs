//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Result row, table, metadata, and beta-prior workflow types.
include!("implementation/types.rs");

// Shared adjusted/public and deferred/private row construction.
include!("implementation/row_builders.rs");

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

// Cook filtering and explicit adjustment updates.
include!("implementation/wald_lrt.rs");

// Validation, column assembly, and result metadata descriptions.
include!("implementation/metadata_validation.rs");
