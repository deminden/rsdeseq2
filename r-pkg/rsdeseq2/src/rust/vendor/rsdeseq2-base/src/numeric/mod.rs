//! Shared numerical building blocks for the statistical stages.

pub mod checked;
pub mod math;
pub mod profile;
pub mod special;
pub mod validate;

pub use profile::{NumericProfile, numeric_profile, set_numeric_profile};
