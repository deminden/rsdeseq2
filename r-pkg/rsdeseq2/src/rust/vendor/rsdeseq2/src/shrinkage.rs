//! Standalone posterior primitives for the ongoing shrinkage-engine migration.
//!
//! These kernels do not yet replace the public R shrinkage methods.
mod adaptive_grid;
mod fixed_mixture;
mod mixture_fit;
mod normal_mixture;
mod svalue;
mod truncated_normal_mixture;
mod truncated_standard;
mod uniform_mixture;

pub use mixture_fit::{MixtureFit, MixtureFitOptions, fit_mixture_proportions};
pub use normal_mixture::{NormalMixture, NormalMixturePosterior};
pub use svalue::svalues;

pub use uniform_mixture::{UniformMixture, UniformMixturePosterior};

pub use truncated_normal_mixture::{TruncatedNormalMixture, TruncatedNormalMixturePosterior};

pub use fixed_mixture::{FixedMixture, FixedMixturePosterior};

pub use adaptive_grid::{
    MixtureComponents, MixtureFamily, MixtureGrid, MixtureGridOptions, MixturePriorPenalty,
    prepare_mixture_grid,
};
