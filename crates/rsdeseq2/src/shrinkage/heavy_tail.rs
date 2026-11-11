//! Standalone heavy-tailed negative-binomial posterior fitting.
mod analysis;
mod beta_binomial;
mod bfgs;
mod compiled_refinement;
mod density;
mod fit;
mod grid;
mod negative_binomial;
mod prior;
mod random;
mod regression;
mod regression_expression;
mod unbounded_lbfgs;
pub use fit::{
    CoefficientMle, CoefficientScale, NegativeBinomialPosterior, NegativeBinomialPosteriorInput,
    NegativeBinomialPosteriorOptions, PosteriorMapDiagnostic, fit_negative_binomial_posterior,
};
pub use fit::{PosteriorCoefficientSummary, PosteriorOffsets};
pub use prior::{CauchyCoefficientPrior, estimate_cauchy_prior_variance};
pub use regression::RegressionMethod;
pub use regression::{
    RegressionLikelihood, RegressionOptimizer, RegressionPosteriorInput,
    RegressionPosteriorOptions, fit_regression_posterior,
};
pub use regression_expression::RestrictedRegressionLikelihood;

pub use random::{PosteriorRng, PosteriorUniformKind};

pub use fit::fit_negative_binomial_posterior_with_rng;
pub use regression::fit_regression_posterior_with_rng;

pub use analysis::{
    PosteriorAnalysis, PosteriorAnalysisOptions, PosteriorGridDiagnostic, PosteriorIntervalKind,
    PosteriorIntervalOptions,
};

pub use compiled_refinement::{NegativeBinomialRefinement, NegativeBinomialRefinementMethod};
