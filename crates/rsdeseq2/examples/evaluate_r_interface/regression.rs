//! Data-only transport for the shared general-regression posterior engine.
use super::{
    bits, heavy_tail::Matrix, parse_numeric_input as number, parse_numeric_inputs as numbers,
};
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{
        CauchyCoefficientPrior, CoefficientMle, CoefficientScale, NegativeBinomialPosterior,
        NegativeBinomialPosteriorOptions, PosteriorOffsets, RegressionLikelihood, RegressionMethod,
        RegressionOptimizer, RegressionPosteriorInput, RegressionPosteriorOptions,
        RestrictedRegressionLikelihood, fit_regression_posterior,
    },
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{error::Error, sync::Arc};
#[derive(Deserialize)]
pub struct Input {
    pub regression_cases: Vec<Case>,
    #[serde(default)]
    pub density_cases: Vec<Leaf>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub observations: Matrix,
    pub design: Matrix,
    pub parameters: Option<Matrix>,
    pub offsets: Option<Matrix>,
    pub weights: Option<Matrix>,
    pub mle: Option<Matrix>,
    pub coefficient: usize,
    pub nuisance: Vec<usize>,
    pub scale: String,
    pub nuisance_scale: String,
    pub prior_df: String,
    pub prior_mean: String,
    pub nuisance_mean: String,
    pub method: String,
    pub optimizer: String,
    pub bounds: [String; 2],
    pub log_link: bool,
    pub threshold: String,
    pub likelihood: Option<Value>,
}
#[derive(Deserialize)]
pub struct Leaf {
    pub name: String,
    pub likelihood: Value,
    pub observations: Vec<String>,
    pub parameters: Vec<String>,
}
pub fn likelihood(v: &Value) -> Result<RestrictedRegressionLikelihood, Box<dyn Error>> {
    let (arguments, expression, constants) = super::likelihood::parts(v)?;
    Ok(RestrictedRegressionLikelihood::new(
        arguments, expression, constants,
    )?)
}
pub fn leaf(c: &Leaf) -> Result<Vec<f64>, Box<dyn Error>> {
    let y = numbers(&c.observations)?;
    let parameters = numbers(&c.parameters)?;
    let n = y.len();
    let design = vec![1.0; n];
    Ok(likelihood(&c.likelihood)?.evaluate(
        &y,
        MatrixView::new(n, 1, &design, MatrixOrder::RowMajor)?,
        &[0.0],
        Some(&parameters),
        &vec![0.0; n],
    )?)
}
pub fn run(c: &Case) -> Result<NegativeBinomialPosterior, Box<dyn Error>> {
    let y = numbers(&c.observations.values)?;
    let x = numbers(&c.design.values)?;
    let parameters = c
        .parameters
        .as_ref()
        .map(|m| numbers(&m.values))
        .transpose()?;
    let offsets = c.offsets.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let weights = c.weights.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let mle = c.mle.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let input = RegressionPosteriorInput {
        observations: c.observations.view(&y)?,
        design: c.design.view(&x)?,
        parameters: c
            .parameters
            .as_ref()
            .zip(parameters.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?,
        offsets: c
            .offsets
            .as_ref()
            .zip(offsets.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?
            .map(PosteriorOffsets::Log)
            .unwrap_or_default(),
        weights: c
            .weights
            .as_ref()
            .zip(weights.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?,
        mle: c
            .mle
            .as_ref()
            .zip(mle.as_ref())
            .map(|(m, v)| {
                if m.columns != 2 || m.rows != c.observations.rows {
                    return Err("invalid MLE fixture shape");
                }
                Ok(CoefficientMle {
                    estimates: &v[..m.rows],
                    standard_errors: &v[m.rows..],
                    scale: CoefficientScale::Natural,
                })
            })
            .transpose()?,
    };
    let options = RegressionPosteriorOptions {
        posterior: NegativeBinomialPosteriorOptions {
            prior: CauchyCoefficientPrior {
                nuisance: c.nuisance.clone(),
                scale: number(&c.scale)?,
                nuisance_scale: number(&c.nuisance_scale)?,
                adapted_variance: None,
            },
            coefficient: Some(c.coefficient),
            threshold: Some(number(&c.threshold)?),
            ..Default::default()
        },
        likelihood: if c.method.starts_with("betabin") {
            RegressionLikelihood::BetaBinomial
        } else {
            match &c.likelihood {
                None => RegressionLikelihood::NegativeBinomial,
                Some(v) => RegressionLikelihood::Restricted(Arc::new(likelihood(v)?)),
            }
        },
        method: match c.method.as_str() {
            "general" => RegressionMethod::General,
            "nbinomR" => RegressionMethod::NegativeBinomialR,
            "betabinR" => RegressionMethod::BetaBinomialR,
            "betabinCR" => RegressionMethod::BetaBinomialCR,
            "betabinC" => RegressionMethod::BetaBinomialC,
            _ => return Err("invalid fixture method".into()),
        },
        optimizer: match c.optimizer.as_str() {
            "BFGS" => RegressionOptimizer::Bfgs,
            "L-BFGS-B" => RegressionOptimizer::Lbfgsb,
            _ => return Err("invalid fixture optimizer".into()),
        },
        prior_df: number(&c.prior_df)?,
        prior_noncentrality: number(&c.prior_mean)?,
        nuisance_mean: number(&c.nuisance_mean)?,
        log_link: c.log_link,
        bounds: (number(&c.bounds[0])?, number(&c.bounds[1])?),
    };
    Ok(fit_regression_posterior(input, &options)?)
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let cases = input
        .regression_cases
        .iter()
        .map(|c| match run(c) {
            Ok(r) => json!({"name": c.name, "value": super::heavy_tail::encode(r, c.coefficient)}),
            Err(e) => json!({"name": c.name, "error": e.to_string()}),
        })
        .collect::<Vec<_>>();
    let leaves = input.density_cases.iter().map(|c| {
        Ok(json!({"name": c.name, "values": leaf(c)?.into_iter().map(bits).collect::<Vec<_>>()}))
    }).collect::<Result<Vec<Value>, Box<dyn Error>>>()?;
    Ok(json!({"cases": cases, "leaves": leaves}))
}
