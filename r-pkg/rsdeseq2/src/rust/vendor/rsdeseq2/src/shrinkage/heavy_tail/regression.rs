//! General per-observation regression likelihoods with shared posterior policy.
use super::{
    density,
    fit::{
        CoefficientMle, NegativeBinomialPosterior, NegativeBinomialPosteriorOptions,
        PosteriorOffsets,
    },
    prior::{CauchyCoefficientPrior, invalid},
    regression_expression::RestrictedRegressionLikelihood,
};
use crate::{DeseqError, matrix::MatrixView};
use nalgebra::DMatrix;
use rcompat_lbfgsb::{Bounds, OptimControl, optim_lbfgsb_with_gradient};
use std::{cell::RefCell, sync::Arc};
#[derive(Clone, Debug)]
pub enum RegressionLikelihood {
    NegativeBinomial,
    Restricted(Arc<RestrictedRegressionLikelihood>),
}
#[derive(Clone, Copy, Debug, Default)]
pub enum RegressionOptimizer {
    #[default]
    Bfgs,
    Lbfgsb,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RegressionMethod {
    #[default]
    General,
    NegativeBinomialR,
}
#[derive(Clone, Copy, Debug)]
pub struct RegressionPosteriorInput<'a> {
    pub observations: MatrixView<'a, f64>,
    pub design: MatrixView<'a, f64>,
    pub parameters: Option<MatrixView<'a, f64>>,
    pub offsets: PosteriorOffsets<'a>,
    pub weights: Option<MatrixView<'a, f64>>,
    pub mle: Option<CoefficientMle<'a>>,
}
#[derive(Clone, Debug)]
pub struct RegressionPosteriorOptions {
    /// Shared scale, coefficient, threshold, adaptation and Laplace controls.
    pub posterior: NegativeBinomialPosteriorOptions,
    pub likelihood: RegressionLikelihood,
    /// The shrunken prior is dt(beta/scale,df,ncp), as in the reference method.
    pub prior_df: f64,
    pub prior_noncentrality: f64,
    pub nuisance_mean: f64,
    pub optimizer: RegressionOptimizer,
    pub method: RegressionMethod,
    pub log_link: bool,
    pub bounds: (f64, f64),
}
impl Default for RegressionPosteriorOptions {
    fn default() -> Self {
        Self {
            posterior: Default::default(),
            likelihood: RegressionLikelihood::NegativeBinomial,
            prior_df: 1.0,
            prior_noncentrality: 0.0,
            nuisance_mean: 0.0,
            optimizer: Default::default(),
            method: Default::default(),
            log_link: true,
            bounds: (f64::NEG_INFINITY, f64::INFINITY),
        }
    }
}
pub fn fit_regression_posterior(
    input: RegressionPosteriorInput<'_>,
    options: &RegressionPosteriorOptions,
) -> Result<NegativeBinomialPosterior, DeseqError> {
    if input.parameters.is_none() {
        return Err(invalid("regression likelihood parameters must be supplied"));
    }
    if options.prior_df.is_nan()
        || options.prior_df <= 0.0
        || !options.prior_noncentrality.is_finite()
        || !options.nuisance_mean.is_finite()
    {
        return Err(invalid("invalid general prior parameters"));
    }
    if let Some(parameters) = input.parameters
        && parameters.n_rows() != input.observations.n_rows()
    {
        return Err(invalid("likelihood parameters must match observation rows"));
    }
    if matches!(options.likelihood, RegressionLikelihood::NegativeBinomial)
        && input.parameters.is_none_or(|p| p.n_cols() != 1)
    {
        return Err(invalid(
            "negative-binomial likelihood requires one dispersion per row",
        ));
    }
    if options.method == RegressionMethod::NegativeBinomialR
        && (!matches!(options.likelihood, RegressionLikelihood::NegativeBinomial)
            || options.prior_df != 1.0
            || options.prior_noncentrality != 0.0
            || options.nuisance_mean != 0.0)
    {
        return Err(invalid(
            "nbinomR requires the zero-centered Cauchy prior and NB likelihood",
        ));
    }
    let dispersions = if options.method == RegressionMethod::NegativeBinomialR {
        let parameters = input
            .parameters
            .ok_or_else(|| invalid("missing NB dispersions"))?;
        (0..parameters.n_rows())
            .map(|i| *parameters.at(i, 0))
            .collect()
    } else {
        vec![1.0; input.observations.n_rows()]
    };
    super::fit::fit_batch(
        super::fit::NegativeBinomialPosteriorInput {
            counts: input.observations,
            design: input.design,
            dispersions: &dispersions,
            offsets: input.offsets,
            weights: input.weights,
            mle: input.mle,
        },
        &options.posterior,
        Some((input.parameters, options)),
    )
}
pub(super) struct GeneralRow<'a> {
    pub y: &'a [f64],
    pub design: MatrixView<'a, f64>,
    pub parameters: Option<Vec<f64>>,
    pub weights: &'a [f64],
    pub offsets: &'a [f64],
    pub prior: &'a CauchyCoefficientPrior,
    pub options: &'a RegressionPosteriorOptions,
    pub error: RefCell<Option<DeseqError>>,
}
impl GeneralRow<'_> {
    pub fn take_error(&self) -> Option<DeseqError> {
        self.error.borrow_mut().take()
    }
    fn value_checked(&self, beta: &[f64]) -> Result<f64, DeseqError> {
        let ll = match &self.options.likelihood {
            RegressionLikelihood::Restricted(f) => f.evaluate(
                self.y,
                self.design,
                beta,
                self.parameters.as_deref(),
                self.offsets,
            )?,
            RegressionLikelihood::NegativeBinomial => {
                let dispersion = self
                    .parameters
                    .as_ref()
                    .and_then(|v| v.first())
                    .copied()
                    .ok_or_else(|| invalid("missing dispersion"))?;
                self.y
                    .iter()
                    .enumerate()
                    .map(|(i, &y)| {
                        let linear = beta
                            .iter()
                            .enumerate()
                            .map(|(j, b)| self.design.at(i, j) * b)
                            .sum::<f64>();
                        density::negative_binomial(
                            y,
                            1.0 / dispersion,
                            (linear + self.offsets[i]).exp(),
                        )
                    })
                    .collect()
            }
        };
        let mut value = -ll.iter().zip(self.weights).map(|(l, w)| l * w).sum::<f64>();
        for &j in &self.prior.nuisance {
            value -= density::normal(
                beta[j],
                self.options.nuisance_mean,
                self.prior.nuisance_scale,
            );
        }
        // The source negative-empty-index convention omits its shrunken sum.
        if !self.prior.nuisance.is_empty() {
            for (j, &b) in beta.iter().enumerate() {
                if !self.prior.nuisance.contains(&j) {
                    value -= density::student(
                        b / self.prior.scale,
                        self.options.prior_df,
                        self.options.prior_noncentrality,
                    );
                }
            }
        }
        Ok(if value.is_finite() { value } else { 1e100 })
    }
    fn value(&self, beta: &[f64]) -> f64 {
        match self.value_checked(beta) {
            Ok(v) => v,
            Err(e) => {
                self.error.borrow_mut().get_or_insert(e);
                1e100
            }
        }
    }
    pub fn objective_gradient(&self, beta: &[f64]) -> (f64, Vec<f64>) {
        let h = 1e-3;
        let gradient = (0..beta.len())
            .map(|j| {
                let mut plus = beta.to_vec();
                let mut minus = beta.to_vec();
                plus[j] += h;
                minus[j] -= h;
                (self.value(&plus) - self.value(&minus)) / (2.0 * h)
            })
            .collect();
        (self.value(beta), gradient)
    }
    pub fn hessian(&self, beta: &[f64], step: f64) -> DMatrix<f64> {
        let p = beta.len();
        let mut h = DMatrix::zeros(p, p);
        for j in 0..p {
            let mut plus = beta.to_vec();
            let mut minus = beta.to_vec();
            plus[j] += step;
            minus[j] -= step;
            let gp = self.objective_gradient(&plus).1;
            let gm = self.objective_gradient(&minus).1;
            for i in 0..p {
                h[(i, j)] = (gp[i] - gm[i]) / (2.0 * step);
            }
        }
        for j in 0..p {
            for i in 0..j {
                let mean = (h[(i, j)] + h[(j, i)]) / 2.0;
                h[(i, j)] = mean;
                h[(j, i)] = mean;
            }
        }
        h
    }
    pub fn optimize(&self, start: Vec<f64>) -> Result<rcompat_lbfgsb::OptimResult, DeseqError> {
        let result = match self.options.optimizer {
            RegressionOptimizer::Bfgs => {
                super::bfgs::minimize(|x| self.objective_gradient(x), start)
            }
            RegressionOptimizer::Lbfgsb => {
                let p = start.len();
                let bounds = Bounds::new(
                    vec![self.options.bounds.0; p],
                    vec![self.options.bounds.1; p],
                )
                .map_err(|e| invalid(&e.to_string()))?;
                optim_lbfgsb_with_gradient(
                    start,
                    bounds,
                    |x| self.value(x),
                    |x| self.objective_gradient(x).1,
                    OptimControl::default_for_dimension(p),
                )
                .map_err(|e| invalid(&e.to_string()))?
            }
        };
        if let Some(e) = self.error.borrow_mut().take() {
            return Err(e);
        }
        Ok(result)
    }
}

pub(super) fn optimize_nb_r(
    row: &super::negative_binomial::Row<'_>,
    options: &RegressionPosteriorOptions,
    start: Vec<f64>,
) -> Result<rcompat_lbfgsb::OptimResult, DeseqError> {
    let constant = -row.objective_gradient(&start).0 - 1.0;
    let f = |b: &[f64]| {
        let (v, g) = row.objective_gradient(b);
        (v + constant, g)
    };
    match options.optimizer {
        RegressionOptimizer::Bfgs => Ok(super::bfgs::minimize(f, start)),
        RegressionOptimizer::Lbfgsb => {
            let p = start.len();
            let bounds = Bounds::new(vec![options.bounds.0; p], vec![options.bounds.1; p])
                .map_err(|e| invalid(&e.to_string()))?;
            optim_lbfgsb_with_gradient(
                start,
                bounds,
                |b| f(b).0,
                |b| f(b).1,
                OptimControl::default_for_dimension(p),
            )
            .map_err(|e| invalid(&e.to_string()))
        }
    }
}
