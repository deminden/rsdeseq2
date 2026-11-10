//! Standalone posterior primitives for the ongoing shrinkage-engine migration.
//!
//! These primitives do not fit a prior or replace the public R shrinkage methods.
mod normal_mixture;
mod svalue;

pub use normal_mixture::{NormalMixture, NormalMixturePosterior};
pub use svalue::svalues;
