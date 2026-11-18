//! R-compatible empirical-Bayes shrinkage: ashr-style mixture posteriors
//! and apeglm-style heavy-tailed posteriors.
//!
//! Part of the rcompat series; used by `rsdeseq2` as `rsdeseq2::shrinkage`.
//! Errors are [`rcompat_numeric::Error`].

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub mod shrinkage;

pub use rcompat_numeric::Error;

// The shrinkage sources use `crate::<module>` paths for these modules.
pub(crate) use rcompat_formula::design;
pub(crate) use rcompat_numeric::{errors, math, matrix, numeric, probability};
