//! Independent transport for selected adaptive shrinkage results.
use super::{bits, formula, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::{design::expression::Data, shrinkage::*};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub ashr_output_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    output: Value,
    threshold: String,
    query_svalues: bool,
    location: Option<super::location_shrinkage::Case>,
}
pub fn run(c: &Case) -> Result<SelectedAdaptiveShrinkageResult, Box<dyn Error>> {
    let output = match formula::value(&c.output).data {
        Data::Number { values, .. } => match values.as_slice() {
            [v] => AdaptiveOutputSelection::level(*v)?,
            _ => return Err("outputlevel requires one numeric level or a text subset".into()),
        },
        Data::Text(names) => {
            AdaptiveOutputSelection::names(names.iter().flatten().map(String::as_str))
        }
        Data::Null | Data::Logical(_) => AdaptiveOutputSelection::default(),
        _ => return Err("outputlevel requires one numeric level or a text subset".into()),
    };
    let x = numbers(&c.estimates)?;
    let s = numbers(&c.standard_errors)?;
    Ok(adaptive_shrinkage_with_output(
        &x,
        &s,
        &if let Some(location) = &c.location {
            super::location_shrinkage::options(location)?
        } else {
            LocationShrinkageOptions {
                supplied_grid: Some(MixtureGrid::from_supplied(
                    MixtureComponents::Normal {
                        means: vec![0.0; 3],
                        sds: vec![0.0, 0.5, 2.0],
                    },
                    vec![0.5, 0.3, 0.2],
                    &MixturePriorPenalty::Uniform,
                )?),
                fixed_prior: true,
                ..Default::default()
            }
        },
        AdaptiveShrinkageQuery {
            estimates: &x,
            standard_errors: &s,
            lfc_threshold: number(&c.threshold)?,
            svalues: c.query_svalues,
        },
        output,
    )?)
}
pub fn encode(r: SelectedAdaptiveShrinkageResult) -> Value {
    let optional = |v: Option<Vec<f64>>| v.map(|v| v.into_iter().map(bits).collect::<Vec<_>>());
    json!({"posterior_mean":optional(r.posterior_mean),"posterior_sd":optional(r.posterior_sd),
        "svalue":optional(r.svalues),"warnings":r.fit.warnings,
        "fitted_prior":r.fitted_prior.then(||super::location_shrinkage::mixture(r.fit.fitted_grid.components,r.fit.fitted_weights))})
}
pub fn evaluate(input: Input) -> Value {
    json!({"cases":input.ashr_output_cases.iter().map(|c| match run(c) {
        Ok(r)=>json!({"name":c.name,"ok":true,"value":encode(r)}),
        Err(_)=>json!({"name":c.name,"ok":false}),
    }).collect::<Vec<_>>()})
}
