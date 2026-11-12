//! Conformance transport for original-result queries of a location fit.
use super::{bits, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::shrinkage::{
    AdaptiveShrinkageQuery, AdaptiveShrinkageResult, LocationLikelihood, LocationShrinkageOptions,
    MixtureComponents, MixtureGrid, MixturePriorPenalty, adaptive_shrinkage,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub adaptive_cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Prior {
    component_class: String,
    weights: Vec<String>,
    means: Option<Vec<String>>,
    sds: Option<Vec<String>>,
    lower: Option<Vec<String>>,
    upper: Option<Vec<String>>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    query_estimates: Vec<String>,
    query_standard_errors: Vec<String>,
    prior: Prior,
    alpha: String,
    df: Option<String>,
    threshold: String,
    svalue: bool,
}
pub fn run(c: &Case) -> Result<AdaptiveShrinkageResult, Box<dyn Error>> {
    let required = |v: &Option<Vec<String>>| {
        v.as_ref()
            .ok_or_else(|| -> Box<dyn Error> { "missing prior field".into() })
            .and_then(|v| numbers(v))
    };
    let components = match c.prior.component_class.as_str() {
        "normalmix" => MixtureComponents::Normal {
            means: required(&c.prior.means)?,
            sds: required(&c.prior.sds)?,
        },
        "unimix" => MixtureComponents::Uniform {
            lower: required(&c.prior.lower)?,
            upper: required(&c.prior.upper)?,
        },
        _ => return Err("unrecognized prior class".into()),
    };
    let options = LocationShrinkageOptions {
        supplied_grid: Some(MixtureGrid::from_supplied(
            components,
            numbers(&c.prior.weights)?,
            &MixturePriorPenalty::Uniform,
        )?),
        fixed_prior: true,
        alpha: number(&c.alpha)?,
        likelihood: c
            .df
            .as_ref()
            .map(|s| number(s).map(LocationLikelihood::StudentT))
            .transpose()?
            .unwrap_or(LocationLikelihood::Normal),
        ..Default::default()
    };
    Ok(adaptive_shrinkage(
        &numbers(&c.estimates)?,
        &numbers(&c.standard_errors)?,
        &options,
        AdaptiveShrinkageQuery {
            estimates: &numbers(&c.query_estimates)?,
            standard_errors: &numbers(&c.query_standard_errors)?,
            lfc_threshold: number(&c.threshold)?,
            svalues: c.svalue,
        },
    )?)
}
pub fn encode(r: AdaptiveShrinkageResult) -> Value {
    let array = |v: Vec<f64>| v.into_iter().map(bits).collect::<Vec<_>>();
    json!({"posterior_mean":array(r.posterior_mean), "posterior_sd":array(r.posterior_sd),
        "svalue":r.svalues.map(array), "warnings":r.fit.warnings})
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let mut results = Vec::new();
    for c in input.adaptive_cases {
        results.push(json!({"name":c.name, "value":encode(run(&c)?)}));
    }
    Ok(json!({"cases":results}))
}
