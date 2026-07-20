//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

use rayon::prelude::*;

// Public options, inputs, outputs, and private diagnostics types.
include!("implementation/types.rs");

// Linear-mu and GLM-mu gene-wise estimation entry points.
include!("implementation/estimate.rs");

// Initial estimates plus line-search and grid dispersion fitting.
include!("implementation/fit.rs");

// Negative-binomial likelihood kernels and derivatives.
include!("implementation/likelihood.rs");

// Cox-Reid adjustment and derivative helpers.
include!("implementation/cox_reid.rs");

// Dispersion prior and posterior objective helpers.
include!("implementation/posterior.rs");

// Validation, compaction, and checked numeric helpers.
include!("implementation/validation.rs");

// Unit tests for gene-wise dispersion numeric guards.
include!("implementation/tests.rs");
