//! Standalone heavy-tailed negative-binomial posterior fitting.
mod beta_binomial;
mod bfgs;
mod density;
mod fit;
mod negative_binomial;
mod prior;
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
