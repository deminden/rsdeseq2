//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Design matrix types, validation, and basic matrix access.
#[path = "implementation/basics.rs"]
mod basics;
pub use basics::*;

// Expanded factor and additive design builders.
#[path = "implementation/expanded.rs"]
mod expanded;
pub use expanded::*;

// Public formula design entry points.
#[path = "implementation/formula_entry.rs"]
mod formula_entry;
pub use formula_entry::*;

// Model-frame construction, validation and alias resolution.
#[path = "implementation/formula_model_frame.rs"]
mod formula_model_frame;
pub use formula_model_frame::*;

// Private state and derived-covariate types shared by the formula implementation.
#[path = "implementation/formula_types.rs"]
mod formula_types;
use formula_types::*;
#[path = "implementation/formula_coding.rs"]
mod formula_coding;
use formula_coding::*;

// Formula right-hand sides and offset extraction.
#[path = "implementation/formula_offsets.rs"]
mod formula_offsets;
use formula_offsets::*;
#[path = "implementation/formula_offset_keys.rs"]
mod formula_offset_keys;
use formula_offset_keys::*;

// Numeric transform discovery, parsing and covariate resolution.
#[path = "implementation/formula_numeric_transforms.rs"]
mod formula_numeric_transforms;
use formula_numeric_transforms::*;

// Translation of validated parser arguments into the shared evaluator.
#[path = "implementation/formula_evaluation.rs"]
mod formula_evaluation;
use formula_evaluation::*;

// Numeric scaling transforms and their options.
#[path = "implementation/formula_scaling.rs"]
mod formula_scaling;
use formula_scaling::*;

// Polynomial transform argument handling.
#[path = "implementation/formula_polynomials.rs"]
mod formula_polynomials;
use formula_polynomials::*;

// Shared transform arguments, options and string parsing.
#[path = "implementation/formula_arguments.rs"]
mod formula_arguments;
use formula_arguments::*;

// Raw and orthogonal polynomial column construction.
#[path = "implementation/formula_polynomial_math.rs"]
mod formula_polynomial_math;
use formula_polynomial_math::*;

// Identity, scalar, binary and power arithmetic transforms.
#[path = "implementation/formula_numeric_arithmetic.rs"]
mod formula_numeric_arithmetic;
use formula_numeric_arithmetic::*;

// Formula factor transforms such as factor(), relevel(), and droplevels().
#[path = "implementation/formula_factors.rs"]
mod formula_factors;
use formula_factors::*;
#[path = "implementation/formula_factor_evaluation.rs"]
mod formula_factor_evaluation;
use formula_factor_evaluation::*;

// Formula term splitting and parenthesized-term expansion.
#[path = "implementation/formula_terms.rs"]
mod formula_terms;
use formula_terms::*;

// Formula state mutation for add/remove term handling.
#[path = "implementation/formula_state.rs"]
mod formula_state;
use formula_state::*;

// Formula design matrix assembly.
#[path = "implementation/formula_assembly.rs"]
mod formula_assembly;
use formula_assembly::*;

// Compatibility labels/order translated to the shared matrix product plan.
#[path = "implementation/formula_matrix_plan.rs"]
mod formula_matrix_plan;
use formula_matrix_plan::*;
#[path = "implementation/formula_ordered_plan.rs"]
mod formula_ordered_plan;
use formula_ordered_plan::*;
#[path = "implementation/formula_fitting.rs"]
mod formula_fitting;
pub use formula_fitting::*;

// Formula-specific validation helpers.
#[path = "implementation/formula_validation.rs"]
mod formula_validation;
use formula_validation::*;

// Shared design validation and coefficient-name helpers.
#[path = "implementation/validation.rs"]
mod validation;
pub(crate) use validation::*;
