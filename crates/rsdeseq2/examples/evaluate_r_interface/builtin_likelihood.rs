//! Primitive built-in parameter transport, independent of the R bridge.
use super::{bits, formula, parse_numeric_inputs as numbers};
use rsdeseq2::{
    design::expression::Data,
    matrix::{MatrixOrder, MatrixView},
    shrinkage::BuiltinLocationLikelihood,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
pub fn model(v: &Value) -> Result<BuiltinLocationLikelihood, Box<dyn Error>> {
    match v["type"].as_str().ok_or("missing built-in type")? {
        "normalmix" => {
            let weights = numbers(&serde_json::from_value::<Vec<String>>(
                v["weights"].clone(),
            )?)?;
            let scales = formula::value(&v["scales"]);
            let (rows, cols) = scales
                .dimensions
                .map(|d| (d[0], d[1]))
                .unwrap_or((1, weights.len()));
            let Data::Number { values, .. } = scales.data else {
                return Err("noise scales require numeric storage".into());
            };
            Ok(BuiltinLocationLikelihood::normal_mixture_by_observation(
                weights,
                MatrixView::new(rows, cols, &values, MatrixOrder::RowMajor)?,
            )?)
        }
        "logF" => Ok(BuiltinLocationLikelihood::log_f(
            numbers(&serde_json::from_value::<Vec<String>>(v["df1"].clone())?)?,
            numbers(&serde_json::from_value::<Vec<String>>(v["df2"].clone())?)?,
        )?),
        _ => Err("unknown built-in likelihood".into()),
    }
}
#[derive(Deserialize)]
pub struct Input {
    pub builtin_likelihood_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub builtin_likelihood: Value,
    arguments: Vec<String>,
    lower: Vec<String>,
    upper: Vec<String>,
}
pub fn run(c: &Case) -> Result<[Vec<f64>; 4], Box<dyn Error>> {
    let l = model(&c.builtin_likelihood)?;
    let x = numbers(&c.arguments)?;
    let a = numbers(&c.lower)?;
    let b = numbers(&c.upper)?;
    if a.len() != x.len() || b.len() != x.len() {
        return Err("likelihood queries require matching lengths".into());
    }
    let mut output = [vec![], vec![], vec![], vec![]];
    for i in 0..x.len() {
        let l = l.for_observation(i, x.len())?;
        output[0].push(l.log_cdf(x[i]));
        output[1].push(l.log_density(x[i]));
        output[2].push(l.truncated_moment(a[i], b[i], 1));
        output[3].push(l.truncated_moment(a[i], b[i], 2));
    }
    Ok(output)
}
pub fn evaluate(input: Input) -> Value {
    json!({"cases":input.builtin_likelihood_cases.iter().map(|c|match run(c){
    Ok([cdf,pdf,first,second])=>json!({"name":c.name,"ok":true,"value":{"log_cdf":cdf.into_iter().map(bits).collect::<Vec<_>>(),"log_density":pdf.into_iter().map(bits).collect::<Vec<_>>(),"first":first.into_iter().map(bits).collect::<Vec<_>>(),"second":second.into_iter().map(bits).collect::<Vec<_>>()}}),
    Err(_)=>json!({"name":c.name,"ok":false})}).collect::<Vec<_>>()})
}
