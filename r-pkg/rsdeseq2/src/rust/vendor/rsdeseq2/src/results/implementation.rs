//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Result row, table, metadata, and beta-prior workflow types.
#[path = "implementation/types.rs"]
mod types;
pub use types::*;

// Shared adjusted/public and deferred/private row construction.
#[path = "implementation/row_builders.rs"]
mod row_builders;
pub use row_builders::*;

// Core Wald builders and caller-supplied expanded beta-prior workflows.
#[path = "implementation/beta_prior_core.rs"]
mod beta_prior_core;
pub use beta_prior_core::*;

// One-factor expanded beta-prior result workflows.
#[path = "implementation/beta_prior_factor.rs"]
mod beta_prior_factor;
pub use beta_prior_factor::*;

// Additive expanded beta-prior result workflows.
#[path = "implementation/beta_prior_additive.rs"]
mod beta_prior_additive;
pub use beta_prior_additive::*;

// Formula parsing adapters and shared beta-prior execution policies.
#[path = "implementation/beta_prior_formula/mod.rs"]
mod beta_prior_formula;
pub use beta_prior_formula::*;

// Cook's replacement merge helpers for beta-prior result workflows.
#[path = "implementation/cooks_replacement.rs"]
mod cooks_replacement;
use cooks_replacement::*;

// Cook filtering and explicit adjustment updates.
#[path = "implementation/wald_lrt.rs"]
mod wald_lrt;
pub use wald_lrt::*;

// Validation, column assembly, and result metadata descriptions.
#[path = "implementation/metadata_validation.rs"]
mod metadata_validation;
use metadata_validation::*;
