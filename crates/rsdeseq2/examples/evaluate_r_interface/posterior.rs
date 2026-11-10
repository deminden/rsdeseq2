use rsdeseq2::shrinkage::{NormalMixture, svalues};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;

#[derive(Deserialize)]
pub struct Input {
    posterior_cases: Vec<Case>,
    svalue_cases: Vec<SvalueCase>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    weights: Vec<String>,
    means: Vec<String>,
    sds: Vec<String>,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
}
#[derive(Deserialize)]
struct SvalueCase {
    name: String,
    rates: Vec<String>,
}
fn numbers(values: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    values
        .iter()
        .map(|v| match v.as_str() {
            "NA" | "NaN" => Ok(f64::NAN),
            _ => Ok(v.parse()?),
        })
        .collect()
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let mut cases = vec![];
    for case in input.posterior_cases {
        let prior = NormalMixture::new(
            numbers(&case.weights)?,
            numbers(&case.means)?,
            numbers(&case.sds)?,
        )?;
        let estimates = numbers(&case.estimates)?;
        let errors = numbers(&case.standard_errors)?;
        if estimates.len() != errors.len() {
            return Err("posterior input lengths differ".into());
        }
        let thresholds = numbers(&case.thresholds)?;
        let mut means = vec![];
        let mut sds = vec![];
        let mut weights = vec![];
        let mut cdf = vec![];
        let mut upper = vec![];
        for (&estimate, &se) in estimates.iter().zip(&errors) {
            let p = prior.posterior(estimate, se)?;
            means.push(super::bits(p.mean()));
            sds.push(super::bits(p.standard_deviation()));
            weights.extend(p.weights().iter().copied().map(super::bits));
            for &threshold in &thresholds {
                cdf.push(super::bits(p.cdf(threshold)));
                upper.push(super::bits(p.survival(threshold)));
            }
        }
        cases.push(
            json!({"name":case.name,"posterior_mean":means,"posterior_sd":sds,
            "posterior_weights":weights,"cdf":cdf,"upper_tail":upper}),
        );
    }
    let mut values = vec![];
    for case in input.svalue_cases {
        let rates = numbers(&case.rates)?
            .into_iter()
            .map(|v| (!v.is_nan()).then_some(v))
            .collect::<Vec<_>>();
        values.push(json!({"name":case.name,"values":svalues(&rates)?.into_iter().map(|v|v.map(super::bits)).collect::<Vec<_>>()}));
    }
    Ok(json!({"cases":cases,"svalues":values}))
}
