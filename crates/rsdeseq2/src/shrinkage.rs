//! Standalone posterior primitives for the ongoing shrinkage-engine migration.
//!
//! These kernels do not yet replace the public R shrinkage methods.
mod mixture_fit;
mod normal_mixture;
mod svalue;

pub use mixture_fit::{MixtureFit, MixtureFitOptions, fit_mixture_proportions};
pub use normal_mixture::{NormalMixture, NormalMixturePosterior};
pub use svalue::svalues;
