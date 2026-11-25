//! R-compatible numeric primitives.
//!
//! Part of the rcompat series (with `rcompat-locfit`, `rcompat-lbfgsb`,
//! `rcompat-formula` and `rcompat-shrinkage`): ports of R behaviour that
//! reproduce R's results, used by `rsdeseq2`.
//!
//! - [`numeric`]: the process-wide [`NumericProfile`] (`portable` pure-Rust
//!   maths, identical on every OS, or the host R build's libm and reductions),
//!   profile-aware elementary functions, special functions and the gamma
//!   distribution, and checked arithmetic.
//! - [`probability`]: R's `pnorm`/`pt`/`pgamma`-style tails and quantiles.
//! - [`math`]: medians, rank tests, a bounded optimizer, trigamma and small
//!   dense linear algebra ([`math::dense`]).
//! - [`matrix`]: a dense row-major matrix and borrowed views.
//! - [`errors`]: the shared [`Error`] type and its coarse [`ErrorKind`].

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code, clippy::undocumented_unsafe_blocks)]

pub mod errors;
pub mod math;
pub mod matrix;
pub mod numeric;
pub mod probability;

pub use errors::{Error, ErrorKind, Requirement};
pub use numeric::{NumericProfile, numeric_profile, set_numeric_profile};
