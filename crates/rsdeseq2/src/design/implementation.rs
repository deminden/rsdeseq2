//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Design matrix types, validation, and basic matrix access.
include!("implementation/basics.rs");

// Expanded factor and additive design builders.
include!("implementation/expanded.rs");

// Formula entry points, offsets, and numeric transforms.
include!("implementation/formula_entry.rs");

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
