//! Standalone heavy-tailed negative-binomial posterior fitting.
use super::try_dense_matrix;
use crate::{Error, matrix::try_alloc_matrix};
use rcompat_numeric::math::dense::DenseMatrix;

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

mod work;
pub use work::{PosteriorCancellation, PosteriorWorkControl};

fn finite_difference_hessian(
    beta: &[f64],
    step: f64,
    gradient: impl Fn(&[f64]) -> Vec<f64>,
) -> Result<DenseMatrix, Error> {
    let p = beta.len();
    let mut hessian = try_dense_matrix(p, p, "posterior Hessian", |_, _| 0.0)?;
    let mut perturbed = try_alloc_matrix(p, 1, "posterior Hessian coefficients")?;
    perturbed.extend_from_slice(beta);
    for j in 0..p {
        perturbed[j] = beta[j] + step;
        let plus = gradient(&perturbed);
        perturbed[j] = beta[j] - step;
        let minus = gradient(&perturbed);
        perturbed[j] = beta[j];
        for i in 0..p {
            hessian[(i, j)] = (plus[i] - minus[i]) / (2.0 * step);
        }
    }
    for j in 0..p {
        for i in 0..j {
            let mean = (hessian[(i, j)] + hessian[(j, i)]) / 2.0;
            hessian[(i, j)] = mean;
            hessian[(j, i)] = mean;
        }
    }
    Ok(hessian)
}
