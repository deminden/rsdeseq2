//! Empirical-Bayes shrinkage engines for rsdeseq2: ashr-style mixture
//! posteriors and apeglm-style heavy-tailed posteriors.
//!
//! This crate is an implementation detail of `rsdeseq2`, which re-exports it
//! as `rsdeseq2::shrinkage`; its own surface has no stability guarantee.
//! Items marked `#[doc(hidden)]` are shared with the sibling crates only.

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub mod shrinkage;

// The shrinkage sources use `crate::<module>` paths for the foundation modules.
pub(crate) use rsdeseq2_base::DeseqError;
pub(crate) use rsdeseq2_base::{errors, math, matrix, numeric, probability};
pub(crate) use rsdeseq2_formula::design;
