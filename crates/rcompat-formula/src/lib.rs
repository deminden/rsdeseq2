//! R formula evaluation and `model.matrix` construction.
//!
//! Part of the rcompat series: parses and evaluates R model formulas
//! (terms, interactions, `offset()`, `I()`, `poly()`, `log()`, factor
//! codings and contrasts, R's collation) against a model frame and builds the
//! design matrices R builds, including DESeq2-style expanded designs. Used by
//! `rsdeseq2` as `rsdeseq2::design`. Errors are [`rcompat_numeric::Error`].

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub mod design;

pub use rcompat_numeric::Error;

// The design sources use `crate::<module>` paths for the numeric modules.
pub(crate) use rcompat_numeric::{errors, math, matrix, numeric};
