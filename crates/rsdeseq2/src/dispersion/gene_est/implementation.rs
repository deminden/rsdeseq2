//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

use crate::normalization::{
    validate_normalization_factor_matrix, validate_normalization_factor_value,
    validate_size_factor_values,
};
use crate::numeric::checked::{self, CheckedExt};
use rayon::prelude::*;

// Public options, inputs, outputs, and private diagnostics types.
#[path = "implementation/types.rs"]
mod types;
pub use types::*;

// Linear-mu and GLM-mu gene-wise estimation entry points.
#[path = "implementation/estimate.rs"]
mod estimate;
pub use estimate::*;

// Initial estimates plus line-search and grid dispersion fitting.
#[path = "implementation/fit.rs"]
mod fit;
pub use fit::*;

// Negative-binomial likelihood kernels and derivatives.
#[path = "implementation/likelihood.rs"]
mod likelihood;
pub use likelihood::*;

// Cox-Reid adjustment and derivative helpers.
#[path = "implementation/cox_reid.rs"]
mod cox_reid;
pub use cox_reid::*;

// Dispersion prior and posterior objective helpers.
#[path = "implementation/posterior.rs"]
mod posterior;
pub use posterior::*;

// Validation, compaction, and checked numeric helpers.
#[path = "implementation/validation.rs"]
mod validation;
use validation::*;

// Unit tests for gene-wise dispersion numeric checks.
#[cfg(test)]
#[path = "implementation/tests.rs"]
mod tests;
