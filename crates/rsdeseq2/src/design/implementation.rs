//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Design matrix types, validation, and basic matrix access.
include!("implementation/basics.rs");

// Expanded factor and additive design builders.
include!("implementation/expanded.rs");

// Public formula design entry points.
include!("implementation/formula_entry.rs");

// Model-frame construction, validation and alias resolution.
include!("implementation/formula_model_frame.rs");

// Private state and derived-covariate types shared by the formula implementation.
include!("implementation/formula_types.rs");

// Formula right-hand sides and offset extraction.
include!("implementation/formula_offsets.rs");

// Numeric transform discovery, parsing and covariate resolution.
include!("implementation/formula_numeric_transforms.rs");

// Numeric scaling transforms and their options.
include!("implementation/formula_scaling.rs");

// Polynomial transform argument handling.
include!("implementation/formula_polynomials.rs");

// Shared transform arguments, options and string parsing.
include!("implementation/formula_arguments.rs");

// Raw and orthogonal polynomial column construction.
include!("implementation/formula_polynomial_math.rs");

// Identity, scalar, binary and power arithmetic transforms.
include!("implementation/formula_numeric_arithmetic.rs");

// Formula factor transforms such as factor(), relevel(), and droplevels().
include!("implementation/formula_factors.rs");

// Formula term splitting and parenthesized-term expansion.
include!("implementation/formula_terms.rs");

// Formula state mutation for add/remove term handling.
include!("implementation/formula_state.rs");

// Formula design matrix assembly.
include!("implementation/formula_assembly.rs");

// Formula-specific validation helpers.
include!("implementation/formula_validation.rs");

// Shared design validation and coefficient-name helpers.
include!("implementation/validation.rs");
