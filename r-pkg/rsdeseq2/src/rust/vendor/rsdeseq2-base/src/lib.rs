//! Numeric foundation shared by the rsdeseq2 crates: errors, numeric
//! profiles, matrices, special functions and distributions.
//!
//! This crate is an implementation detail of `rsdeseq2`; its public surface
//! has no stability guarantee of its own. Items marked `#[doc(hidden)]` are
//! shared with the sibling crates only. Use the `rsdeseq2` crate, which
//! re-exports the supported paths.

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub mod errors;
pub mod math;
pub mod matrix;
pub mod nb;
#[doc(hidden)]
pub mod numeric;
#[doc(hidden)]
pub mod probability;

pub use errors::DeseqError;
pub use numeric::{NumericProfile, numeric_profile, set_numeric_profile};
