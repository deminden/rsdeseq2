//! Estimated-mode conformance transport; both stages use standalone Rust.
use super::{formula, location_shrinkage, parse_numeric_inputs as numbers};
use rsdeseq2::{design::expression::Data, shrinkage::*};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub ashr_mode_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    #[serde(flatten)]
    pub data: location_shrinkage::Case,
    mode_input: Value,
}
pub fn search(c: &Case) -> Result<PriorModeSearch, Box<dyn Error>> {
    Ok(match formula::value(&c.mode_input).data {
        Data::Text(values)
            if values.len() == 2 || values.iter().flatten().any(|v| v == "estimate") =>
        {
            PriorModeSearch::ObservationRange
        }
        Data::Number { values, .. } if values.len() == 2 => {
            if values[0].is_nan() {
                PriorModeSearch::ObservationRange
            } else {
                PriorModeSearch::ExpandedObservationRange(values[0])
            }
        }
        _ => return Err("invalid estimated-mode request".into()),
    })
}
pub fn run(c: &Case) -> Result<(LocationModeEstimate, LocationShrinkage), Box<dyn Error>> {
    let mut options = location_shrinkage::options(&c.data)?;
    let x = numbers(&c.data.estimates)?;
    let se = numbers(&c.data.standard_errors)?;
    let mode = estimate_location_mode(&x, &se, &options, search(c)?)?;
    options.grid.mode = mode.mode;
    let mut result = shrink_location_batch(&x, &se, &options)?;
    let mut warnings = mode.warnings.clone();
    warnings.append(&mut result.warnings);
    result.warnings = warnings;
    Ok((mode, result))
}
pub fn encode(result: LocationShrinkage) -> Value {
    let mut v = location_shrinkage::encode(result);
    let object = v.as_object_mut().unwrap();
    object.retain(|name, _| {
        [
            "posterior_mean",
            "posterior_sd",
            "svalues",
            "fitted_prior",
            "warnings",
        ]
        .contains(&name.as_str())
    });
    let svalue = object.remove("svalues").unwrap();
    object.insert("svalue".into(), svalue);
    v
}
pub fn evaluate(input: Input) -> Value {
    json!({"cases":input.ashr_mode_cases.iter().map(|c| match run(c){
        Ok((_,r))=>json!({"name":c.data.name,"ok":true,"value":encode(r)}),
        Err(_)=>json!({"name":c.data.name,"ok":false}),
    }).collect::<Vec<_>>()})
}
