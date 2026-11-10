//! Standalone heavy-tailed negative-binomial posterior fitting.
mod fit;
mod negative_binomial;
mod prior;
mod unbounded_lbfgs;
pub use fit::{
    CoefficientMle, CoefficientScale, NegativeBinomialPosterior, NegativeBinomialPosteriorInput,
    NegativeBinomialPosteriorOptions, PosteriorMapDiagnostic, fit_negative_binomial_posterior,
};
pub use fit::{PosteriorCoefficientSummary, PosteriorOffsets};
pub use prior::{CauchyCoefficientPrior, estimate_cauchy_prior_variance};
