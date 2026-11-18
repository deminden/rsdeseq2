//! DESeq2-compatible RNA-seq differential expression in Rust.
//!
//! The library provides normalization, dispersion estimation, negative-binomial
//! GLMs, Wald and likelihood-ratio tests, contrasts, Cook's diagnostics and
//! replacement/refitting, independent filtering, transformations and supported
//! normal, ashr and apeglm shrinkage. The CLI and optional R interface share
//! this numerical core; standalone Rust execution requires no R runtime.
//!
//! Start with [`prelude`] and [`core::DeseqBuilder`]. Individual statistical
//! stages are available through their public modules. Compatibility is bounded
//! by the documented workflow contracts; the complete Bioconductor object API
//! is outside scope. See the repository's compatibility and testing guides for
//! supported settings and independent reference validation.

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub(crate) mod all_zero;
#[cfg(feature = "cli")]
pub mod cli;
pub mod contrasts;
pub mod cooks;
pub mod core;
pub mod design;
pub mod diagnostics;
pub mod dispersion;
pub use rsdeseq2_base::errors;
pub mod glm;
pub mod independent_filtering;
pub mod io;
pub use rsdeseq2_base::math;
pub use rsdeseq2_base::matrix;
pub mod multiple_testing;
pub mod normalization;
pub use numeric::{NumericProfile, numeric_profile, set_numeric_profile};
pub(crate) use rsdeseq2_base::numeric;
pub mod options;
pub mod prelude;
pub(crate) use rsdeseq2_base::probability;
pub mod results;
#[cfg(feature = "shrinkage")]
pub mod shrinkage;
pub mod transform;

pub use crate::errors::DeseqError;
