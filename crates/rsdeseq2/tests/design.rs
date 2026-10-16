//! Integration coverage grouped by workflow; shared fixtures remain in this test crate.

use rsdeseq2::design::formula_has_offset_terms;
use rsdeseq2::prelude::*;
use std::assert_matches;

#[path = "design/matrix_and_expansion.rs"]
mod matrix_and_expansion;

#[path = "design/interactions.rs"]
mod interactions;

#[path = "design/factor_transforms.rs"]
mod factor_transforms;

#[path = "design/model_frame.rs"]
mod model_frame;

#[path = "design/aliases.rs"]
mod aliases;

#[path = "design/formula_operators.rs"]
mod formula_operators;

#[path = "design/grouped_terms.rs"]
mod grouped_terms;

#[path = "design/numeric_arithmetic.rs"]
mod numeric_arithmetic;

#[path = "design/numeric_functions.rs"]
mod numeric_functions;

#[path = "design/polynomials.rs"]
mod polynomials;

#[path = "design/offsets_and_validation.rs"]
mod offsets_and_validation;
