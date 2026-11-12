//! Data-only transport for independent accelerated-EM conformance panels.
use super::{bits, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{MixtureFit, MixtureFitOptions, ashr_mixture_fit_options, fit_mixture_proportions},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub em_fit_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub optimizer: String,
    pub observations: usize,
    pub components: usize,
    pub log_likelihoods: Vec<String>,
    pub observation_weights: Option<Vec<String>>,
    pub dirichlet: Vec<String>,
    pub initial_weights: Option<Vec<String>>,
    pub control: Value,
}
fn primitive(v: &Value) -> Result<rsdeseq2::design::expression::Value, Box<dyn Error>> {
    use rsdeseq2::design::expression::{Data, Value as Primitive};
    Ok(match v {
        Value::Null => Primitive::new(Data::Null),
        Value::Bool(v) => Primitive::new(Data::Logical(vec![Some(*v)])),
        Value::Number(v) => Primitive::numbers(vec![v.as_f64().ok_or("invalid numeric control")?]),
        Value::String(v) if matches!(v.as_str(), "Inf" | "-Inf" | "NaN") => {
            Primitive::numbers(vec![number(v)?])
        }
        Value::String(v) => Primitive::new(Data::Text(vec![Some(v.clone())])),
        Value::Array(v) if !v.is_empty() && v.iter().all(Value::is_number) => Primitive::numbers(
            v.iter()
                .map(|n| n.as_f64().ok_or("invalid numeric control"))
                .collect::<Result<Vec<_>, _>>()?,
        ),
        Value::Array(v) => Primitive::new(Data::List(
            v.iter().map(primitive).collect::<Result<Vec<_>, _>>()?,
        )),
        Value::Object(v) => {
            let mut result = Primitive::new(Data::List(
                v.values().map(primitive).collect::<Result<Vec<_>, _>>()?,
            ));
            result.names = Some(v.keys().cloned().collect());
            result
        }
    })
}
pub(crate) fn fit_options(
    method: &str,
    c: &Value,
    fixed_prior: bool,
) -> Result<MixtureFitOptions, Box<dyn Error>> {
    let control = if c.is_null() {
        None
    } else {
        Some(primitive(c)?)
    };
    ashr_mixture_fit_options(method, control.as_ref(), fixed_prior).map_err(|error| match error {
        rsdeseq2::DeseqError::InvalidOptions { reason } => reason.into(),
        error => error.into(),
    })
}
pub fn run(c: &Case) -> Result<MixtureFit, Box<dyn Error>> {
    let ll = numbers(&c.log_likelihoods)?;
    let weights = c
        .observation_weights
        .as_ref()
        .map(|v| numbers(v))
        .transpose()?;
    let prior = numbers(&c.dirichlet)?;
    let initial = c.initial_weights.as_ref().map(|v| numbers(v)).transpose()?;
    let options = fit_options(&c.optimizer, &c.control, false)?;
    Ok(fit_mixture_proportions(
        MatrixView::new(c.observations, c.components, &ll, MatrixOrder::RowMajor)?,
        weights.as_deref(),
        &prior,
        initial.as_deref(),
        options,
    )?)
}
pub fn encode(f: MixtureFit) -> Value {
    let array = |v: Vec<f64>| v.into_iter().map(bits).collect::<Vec<_>>();
    json!({"weights":array(f.weights),"log_likelihood":bits(f.log_likelihood),"penalized_log_likelihood":bits(f.penalized_log_likelihood),"converged":f.converged,"iterations":f.iterations,"evaluations":f.function_evaluations,"residual":bits(f.dual_residual),"weight_sum":bits(f.optimizer_weight_sum),"posterior_concentrations":f.posterior_concentrations.map(array),"optimizer_objective":f.optimizer_objective.map(bits),"message":f.message})
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let cases = input
        .em_fit_cases
        .iter()
        .map(|c| match run(c) {
            Ok(f) => json!({"name":c.name,"value":encode(f)}),
            Err(e) => json!({"name":c.name,"error":e.to_string()}),
        })
        .collect::<Vec<_>>();
    Ok(json!({"cases":cases}))
}
