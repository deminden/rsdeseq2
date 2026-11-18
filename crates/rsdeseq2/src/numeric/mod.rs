//! Shared numerical building blocks for the statistical stages.

pub(crate) mod checked;
pub(crate) mod math;
pub mod profile;
pub(crate) mod special;
pub(crate) mod validate;

pub use profile::{NumericProfile, numeric_profile, set_numeric_profile};
