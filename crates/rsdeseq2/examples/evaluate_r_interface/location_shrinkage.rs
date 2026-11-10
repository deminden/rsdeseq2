//! Whole-batch standalone conformance transport. Numerical work is in the core.
use super::{bits, parse_numeric_input as number, parse_numeric_inputs as numbers};
use rsdeseq2::shrinkage::{
    LocationLikelihood, LocationShrinkage, LocationShrinkageOptions, MixtureComponents,
    MixtureGrid, MixtureGridOptions, MixturePriorPenalty, shrink_location_batch,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    location_cases: Vec<Case>,
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
struct Case {
    name: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
    family: String,
    point_mass: bool,
    nullweight: Option<String>,
    scales: Option<Vec<String>>,
    mode: String,
    range: Vec<String>,
    alpha: String,
    prune_threshold: String,
    observation_weights: Option<Vec<String>>,
    fixed_prior: bool,
    supplied_prior: Option<Prior>,
    df: Option<String>,
    custom_likelihood: Option<Value>,
}
fn array(values: Vec<f64>) -> Value {
    json!(values.into_iter().map(bits).collect::<Vec<_>>())
}
fn mixture(components: MixtureComponents, weights: Vec<f64>) -> Value {
    let (class, means, sds, lower, upper) = match components {
        MixtureComponents::Normal { means, sds } => {
            ("normalmix", Some(means), Some(sds), None, None)
        }
        MixtureComponents::Uniform { lower, upper } => {
            ("unimix", None, None, Some(lower), Some(upper))
        }
        MixtureComponents::TruncatedNormal {
            means,
            sds,
            lower,
            upper,
        } => (
            "tnormalmix",
            Some(means),
            Some(sds),
            Some(lower),
            Some(upper),
        ),
    };
    json!({"component_class":[class],"weights":array(weights),"means":means.map(array),"sds":sds.map(array),"lower":lower.map(array),"upper":upper.map(array)})
}
fn encode(r: LocationShrinkage) -> Value {
    let fit=r.fit.map(|f|json!({"weights":array(f.weights),"log_likelihood":array(vec![f.log_likelihood]),"penalized_log_likelihood":array(vec![f.penalized_log_likelihood]),"dual_residual":array(vec![f.dual_residual]),"optimizer_weight_sum":array(vec![f.optimizer_weight_sum]),"function_evaluations":array(vec![f.function_evaluations as f64]),"converged":[f.converged],"message":[f.message]}));
    json!({"fitted_prior":mixture(r.fitted_grid.components,r.fitted_weights),"fit":fit,"posterior_mean":array(r.posterior_mean),"posterior_sd":array(r.posterior_sd),"negative_probability":array(r.negative_probability),"positive_probability":array(r.positive_probability),"local_false_discovery_rate":array(r.local_false_discovery_rate),"local_false_sign_rate":array(r.local_false_sign_rate),"svalues":array(r.svalues),"qvalues":array(r.qvalues),"cdf":array(r.cdf),"log_likelihood":array(vec![r.log_likelihood]),"warnings":r.warnings})
}
fn run(c: &Case) -> Result<LocationShrinkage, Box<dyn Error>> {
    if c.range.len() != 2 {
        return Err("location range requires two bounds".into());
    }
    let penalty = c
        .nullweight
        .as_ref()
        .map(|s| number(s).map(MixturePriorPenalty::NullBiased))
        .transpose()?
        .unwrap_or(MixturePriorPenalty::Uniform);
    let supplied_grid = if let Some(p) = &c.supplied_prior {
        let weights = numbers(&p.weights)?;
        let mut dirichlet = vec![1.0; weights.len()];
        if let Some(w) = &c.nullweight
            && let Some(first) = dirichlet.first_mut()
        {
            *first = number(w)?;
        }
        let required = |v: &Option<Vec<String>>| {
            v.as_ref()
                .ok_or_else(|| -> Box<dyn Error> {
                    "missing supplied-prior component field".into()
                })
                .and_then(|v| numbers(v))
        };
        let components = match p.component_class.as_str() {
            "normalmix" => MixtureComponents::Normal {
                means: required(&p.means)?,
                sds: required(&p.sds)?,
            },
            "unimix" => MixtureComponents::Uniform {
                lower: required(&p.lower)?,
                upper: required(&p.upper)?,
            },
            "tnormalmix" => MixtureComponents::TruncatedNormal {
                means: required(&p.means)?,
                sds: required(&p.sds)?,
                lower: required(&p.lower)?,
                upper: required(&p.upper)?,
            },
            _ => return Err("unknown supplied-prior class".into()),
        };
        Some(MixtureGrid {
            components,
            initial_weights: weights,
            dirichlet,
            warnings: vec![],
        })
    } else {
        None
    };
    let options = LocationShrinkageOptions {
        grid: MixtureGridOptions {
            family: c.family.parse()?,
            scales: c.scales.as_ref().map(|v| numbers(v)).transpose()?,
            point_mass: c.point_mass,
            mode: number(&c.mode)?,
            range: (number(&c.range[0])?, number(&c.range[1])?),
            penalty,
            ..Default::default()
        },
        supplied_grid,
        fixed_prior: c.fixed_prior,
        likelihood: if let Some(value) = &c.custom_likelihood {
            LocationLikelihood::Restricted(std::sync::Arc::new(super::likelihood::restricted(
                value,
            )?))
        } else {
            c.df.as_ref()
                .map(|s| number(s).map(LocationLikelihood::StudentT))
                .transpose()?
                .unwrap_or(LocationLikelihood::Normal)
        },
        alpha: number(&c.alpha)?,
        prune_threshold: number(&c.prune_threshold)?,
        observation_weights: c
            .observation_weights
            .as_ref()
            .map(|v| numbers(v))
            .transpose()?,
        thresholds: numbers(&c.thresholds)?,
        ..Default::default()
    };
    Ok(shrink_location_batch(
        &numbers(&c.estimates)?,
        &numbers(&c.standard_errors)?,
        &options,
    )?)
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    Ok(
        json!({"cases":input.location_cases.iter().map(|c|match run(c){Ok(r)=>json!({"name":c.name,"ok":true,"value":encode(r)}),Err(_)=>json!({"name":c.name,"ok":false})}).collect::<Vec<_>>()}),
    )
}
