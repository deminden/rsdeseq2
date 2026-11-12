//! Formula adapters share one execution path for each normalization policy.
use super::*;

mod input;
use input::{FormulaPriorInput, FormulaPriorNormalizedInput};

mod size_factors;
use size_factors::*;
#[path = "normalized.rs"]
mod normalization_factors;
use normalization_factors::*;

mod primitive;
pub use primitive::*;
mod model_frame;
pub use model_frame::*;
