//! Data-only transport for independent accelerated-EM conformance panels.
use super::{bits, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{
        InteriorPointOptions, MixtureFit, MixtureFitOptions, MixtureOptimizer,
        PolynomialExtrapolation, SquaremOptions, fit_mixture_proportions,
    },
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
fn control_number(c: &Value, n: &str, d: f64) -> Result<f64, Box<dyn Error>> {
    match c.get(n) {
        None => Ok(d),
        Some(Value::Number(v)) => v.as_f64().ok_or("invalid numeric control".into()),
        Some(Value::String(v)) => number(v),
        _ => Err("invalid numeric control".into()),
    }
}
pub(crate) fn sqp_controls(c: &Value, o: &mut MixtureFitOptions) -> Result<(), Box<dyn Error>> {
    let Some(values) = c.as_object() else {
        if c.is_null() || c.as_array().is_some_and(|x| x.is_empty()) {
            return Ok(());
        }
        return Err("SQP control requires a named list".into());
    };
    for (name, value) in values {
        let count = || -> Result<usize, Box<dyn Error>> {
            let x = control_number(c, name, 0.0)?;
            if !x.is_finite() || x < 0.0 || x.fract() != 0.0 || x > i32::MAX as f64 {
                return Err("invalid SQP integer control".into());
            }
            Ok(x as usize)
        };
        match name.as_str() {
            "verbose" | "normalize.rows" => {
                if !value.is_boolean() {
                    return Err("SQP flag requires logical storage".into());
                }
            }
            "eps" => o.likelihood_offset = control_number(c, name, 1e-6)?,
            "numiter.em" => o.em_iterations = count()?,
            "maxiter.sqp" => o.max_iterations = count()?,
            "maxiter.activeset" => {
                if value.is_null() {
                    o.sqp.automatic_active_iterations = true;
                } else {
                    o.active_iterations = count()?;
                    o.sqp.automatic_active_iterations = false;
                }
            }
            "convtol.sqp" => o.gradient_tolerance = control_number(c, name, 1e-8)?,
            "zero.threshold.solution" => o.zero_threshold = control_number(c, name, 1e-8)?,
            "convtol.activeset" => o.sqp.active_tolerance = control_number(c, name, 1e-10)?,
            "zero.threshold.searchdir" => {
                o.sqp.search_direction_threshold = Some(control_number(c, name, 1e-14)?)
            }
            "suffdecr.linesearch" => o.sqp.sufficient_decrease = control_number(c, name, 0.01)?,
            "stepsizereduce" => o.sqp.step_reduction = control_number(c, name, 0.75)?,
            "minstepsize" => o.sqp.minimum_step = control_number(c, name, 1e-8)?,
            "identity.contrib.increase" => o.sqp.ridge_multiplier = control_number(c, name, 10.0)?,
            "tol.svd" => o.sqp.svd_tolerance = Some(control_number(c, name, 0.0)?),
            _ => return Err(format!("unknown SQP control '{name}'").into()),
        }
    }
    Ok(())
}
pub(crate) fn policy(method: &str, c: &Value) -> Result<MixtureOptimizer, Box<dyn Error>> {
    let count = |n: &str, d: usize| -> Result<usize, Box<dyn Error>> {
        let v = control_number(c, n, d as f64)?;
        if !v.is_finite() || v < 0.0 || v.fract() != 0.0 || v >= usize::MAX as f64 {
            return Err("invalid integer EM control".into());
        }
        Ok(v as usize)
    };
    if method == "mixIP" {
        return Ok(MixtureOptimizer::InteriorPoint(InteriorPointOptions {
            relative_tolerance: control_number(c, "rtol", 1e-6)?,
            ..Default::default()
        }));
    }
    if method == "cxxMixSquarem" && matches!(c.get("method"), Some(Value::String(_))) {
        return Err("compiled-style EM requires a numeric method control".into());
    }
    let method_code = if matches!(c.get("method"), Some(Value::String(_))) {
        3.0
    } else {
        control_number(c, "method", 3.0)?
    };
    let options = SquaremOptions {
        tolerance: c
            .get("tol")
            .map(|_| control_number(c, "tol", 1e-7))
            .transpose()?,
        max_evaluations: count("maxiter", 5000)?,
        method: if [1.0, 2.0, 3.0].contains(&method_code) {
            method_code as u8
        } else {
            3
        },
        minimum_step: control_number(c, "step.min0", 1.0)?,
        maximum_step: control_number(c, "step.max0", 1.0)?,
        step_multiplier: control_number(c, "mstep", 4.0)?,
        objective_increase: control_number(c, "objfn.inc", 1.0)?,
        square: c.get("square").and_then(Value::as_bool).unwrap_or(true),
        polynomial_order: count("K", 1)?,
        polynomial_method: if c.get("method").and_then(Value::as_str) == Some("mpe") {
            PolynomialExtrapolation::MinimalPolynomial
        } else {
            PolynomialExtrapolation::ReducedRank
        },
    };
    Ok(match method {
        "mixEM" => MixtureOptimizer::Em(options),
        "w_mixEM" => MixtureOptimizer::WeightedEm(options),
        "mixVBEM" => MixtureOptimizer::VariationalEm(options),
        "cxxMixSquarem" => MixtureOptimizer::CxxEm(options),
        _ => return Err("unsupported EM method".into()),
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
    let mut options = if c.optimizer == "mixSQP" {
        rsdeseq2::shrinkage::LocationShrinkageOptions::default().fit_options
    } else {
        MixtureFitOptions {
            optimizer: policy(&c.optimizer, &c.control)?,
            ..Default::default()
        }
    };
    if c.optimizer == "mixSQP" {
        sqp_controls(&c.control, &mut options)?;
    }
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
