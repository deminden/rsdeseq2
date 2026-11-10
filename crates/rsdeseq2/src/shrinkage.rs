//! Standalone posterior primitives for the ongoing shrinkage-engine migration.
//!
//! Default R shrinkage routes share these implementations. Advanced optional
//! methods remain under development; see the documented migration scope.
mod adaptive_grid;
mod adaptive_result;
mod beta_normalization;
mod component_summary;
mod expression_likelihood;
mod expression_quadrature;
mod fixed_mixture;
mod gamma_normalization;
mod heavy_tail;
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
pub use heavy_tail::RegressionMethod;
pub use heavy_tail::RestrictedRegressionLikelihood;
pub use heavy_tail::{CauchyCoefficientPrior, estimate_cauchy_prior_variance};
pub use heavy_tail::{
    CoefficientMle, CoefficientScale, NegativeBinomialPosterior, NegativeBinomialPosteriorInput,
    NegativeBinomialPosteriorOptions, PosteriorMapDiagnostic, fit_negative_binomial_posterior,
};
pub use heavy_tail::{PosteriorCoefficientSummary, PosteriorOffsets};
pub use heavy_tail::{
    RegressionLikelihood, RegressionOptimizer, RegressionPosteriorInput,
    RegressionPosteriorOptions, fit_regression_posterior,
};
