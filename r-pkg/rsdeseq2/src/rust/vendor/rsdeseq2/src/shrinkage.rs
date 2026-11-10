//! Standalone posterior primitives for the ongoing shrinkage-engine migration.
//!
//! These kernels do not yet replace the public R shrinkage methods.
mod adaptive_grid;
mod adaptive_result;
mod beta_normalization;
mod component_summary;
mod expression_likelihood;
mod expression_quadrature;
mod fixed_mixture;
mod gamma_normalization;
mod likelihood;
mod likelihood_diagnostics;
mod location_batch;
mod location_missing;
mod location_query;
mod mixture_fit;
mod normal_mixture;
mod student_t;
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
    AdaptiveMethod, MixtureComponents, MixtureFamily, MixtureGrid, MixtureGridOptions,
    MixturePriorPenalty, prepare_mixture_grid,
};

pub use location_batch::{LocationShrinkage, LocationShrinkageOptions, shrink_location_batch};

pub use expression_likelihood::{RestrictedFunction, RestrictedLikelihood};

pub use likelihood::LocationLikelihood;

pub use location_query::{LocationCdfOptions, location_cdf_batch, threshold_svalues};

pub use adaptive_result::{AdaptiveShrinkageQuery, AdaptiveShrinkageResult, adaptive_shrinkage};
