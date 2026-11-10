//! Primitive conformance transport for standalone Cauchy/NB posterior batches.
use super::{bits, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::{
    DeseqError,
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{
        CauchyCoefficientPrior, CoefficientMle, CoefficientScale, NegativeBinomialPosterior,
        NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions, PosteriorOffsets,
        estimate_cauchy_prior_variance, fit_negative_binomial_posterior,
    },
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub heavy_tail_cases: Vec<Case>,
    #[serde(default)]
    pub prior_variance_cases: Vec<PriorCase>,
}
#[derive(Deserialize)]
pub struct PriorCase {
    pub name: String,
    pub estimates: Vec<String>,
    pub standard_errors: Vec<String>,
}
#[derive(Deserialize)]
pub struct Matrix {
    pub rows: usize,
    pub columns: usize,
    pub values: Vec<String>,
}
impl Matrix {
    pub fn view<'a>(&self, values: &'a [f64]) -> Result<MatrixView<'a, f64>, DeseqError> {
        MatrixView::new(self.rows, self.columns, values, MatrixOrder::ColumnMajor)
    }
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub counts: Matrix,
    pub design: Matrix,
    pub dispersions: Vec<String>,
    pub offsets: Option<Matrix>,
    pub size_factors: Option<Vec<String>>,
    pub normalization_factors: Option<Matrix>,
    pub weights: Option<Matrix>,
    pub mle: Option<Matrix>,
    pub mle_scale: Option<String>,
    pub threshold_scale: Option<String>,
    pub coefficient: usize,
    pub scale: String,
    pub nuisance: Vec<usize>,
    pub nuisance_scale: String,
    pub threshold: String,
    pub flip_sign: bool,
    pub map_only: bool,
}
pub fn run(c: &Case) -> Result<NegativeBinomialPosterior, Box<dyn Error>> {
    let y = numbers(&c.counts.values)?;
    let x = numbers(&c.design.values)?;
    let d = numbers(&c.dispersions)?;
    let offsets = c.offsets.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let size_factors = c.size_factors.as_ref().map(|v| numbers(v)).transpose()?;
    let factors = c
        .normalization_factors
        .as_ref()
        .map(|m| numbers(&m.values))
        .transpose()?;
    let weights = c.weights.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let mle = c.mle.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let options = NegativeBinomialPosteriorOptions {
        prior: CauchyCoefficientPrior {
            nuisance: c.nuisance.clone(),
            scale: number(&c.scale)?,
            nuisance_scale: number(&c.nuisance_scale)?,
            adapted_variance: None,
        },
        coefficient: Some(c.coefficient),
        threshold: Some(if c.threshold_scale.as_deref() == Some("log2") {
            CoefficientScale::Log2.to_natural(number(&c.threshold)?)
        } else {
            number(&c.threshold)?
        }),
        flip_sign: c.flip_sign,
        map_only: c.map_only,
        check_opposite_start: !c.map_only,
        ..Default::default()
    };
    let input = NegativeBinomialPosteriorInput {
        counts: c.counts.view(&y)?,
        design: c.design.view(&x)?,
        dispersions: &d,
        offsets: if let Some(s) = &size_factors {
            PosteriorOffsets::SizeFactors(s)
        } else if let Some((m, v)) = c.normalization_factors.as_ref().zip(factors.as_ref()) {
            PosteriorOffsets::Factors(m.view(v)?)
        } else {
            c.offsets
                .as_ref()
                .zip(offsets.as_ref())
                .map(|(m, v)| m.view(v))
                .transpose()?
                .map(PosteriorOffsets::Log)
                .unwrap_or_default()
        },
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
                if m.columns != 2 || m.rows != c.counts.rows {
                    return Err("invalid MLE fixture shape");
                }
                Ok(CoefficientMle {
                    estimates: &v[..m.rows],
                    standard_errors: &v[m.rows..],
                    scale: if c.mle_scale.as_deref() == Some("log2") {
                        CoefficientScale::Log2
                    } else {
                        CoefficientScale::Natural
                    },
                })
            })
            .transpose()?,
    };
    Ok(fit_negative_binomial_posterior(input, &options)?)
}
pub fn encode(r: NegativeBinomialPosterior, coefficient: usize) -> Value {
    let array = |v: Vec<f64>| v.into_iter().map(bits).collect::<Vec<_>>();
    let optional = |v: Vec<Option<f64>>| v.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>();
    let summary = r
        .coefficient_summary(coefficient, CoefficientScale::Log2)
        .unwrap();
    json!({"mode":array(r.mode.into_values()),"sd":r.standard_deviation.map(|v|array(v.into_values())),
        "fsr":r.false_sign_rate.map(optional),"sign_svalues":r.sign_svalues.map(optional),
        "fsos":r.false_sign_or_small_rate.map(optional),"threshold_svalues":r.threshold_svalues.map(optional),
        "prior_scale":array(vec![r.prior.scale]),"prior_variance":r.prior.adapted_variance.map(|v|array(vec![v])),
        "nuisance":r.prior.nuisance,"warnings":r.warnings,
        "summary_mode":array(summary.mode),"summary_sd":optional(summary.standard_deviation),
        "diagnostics":r.diagnostics.into_iter().map(|d|json!({"convergence":d.convergence,
            "function_evaluations":d.function_evaluations,"opposite_start_disagrees":d.opposite_start_disagrees,
            "laplace_valid":d.laplace_valid})).collect::<Vec<_>>()})
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let priors=input.prior_variance_cases.into_iter().map(|c|{
        Ok(json!({"name":c.name,"variance":bits(estimate_cauchy_prior_variance(&numbers(&c.estimates)?,&numbers(&c.standard_errors)?)?)}))
    }).collect::<Result<Vec<Value>,Box<dyn Error>>>()?;
    let cases = input
        .heavy_tail_cases
        .into_iter()
        .map(|c| Ok(json!({"name":c.name,"value":encode(run(&c)?,c.coefficient)})))
        .collect::<Result<Vec<Value>, Box<dyn Error>>>()?;
    Ok(json!({"cases":cases,"priors":priors}))
}
