use rsdeseq2::shrinkage::{
    MixtureComponents, MixtureGridOptions, MixturePriorPenalty, prepare_mixture_grid,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    grid_cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    family: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    scales: Option<Vec<String>>,
    point_mass: bool,
    mode: String,
    range: Vec<String>,
    nullweight: Option<String>,
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let mut result = vec![];
    for c in input.grid_cases {
        if c.range.len() != 2 {
            return Err("grid range requires two bounds".into());
        }
        let options = MixtureGridOptions {
            family: c.family.parse()?,
            scales: c
                .scales
                .as_ref()
                .map(|v| super::parse_numeric_inputs(v))
                .transpose()?,
            point_mass: c.point_mass,
            mode: super::parse_numeric_input(&c.mode)?,
            range: (
                super::parse_numeric_input(&c.range[0])?,
                super::parse_numeric_input(&c.range[1])?,
            ),
            penalty: match c.nullweight {
                Some(v) => MixturePriorPenalty::NullBiased(super::parse_numeric_input(&v)?),
                None => MixturePriorPenalty::Uniform,
            },
            ..Default::default()
        };
        let g = prepare_mixture_grid(
            &super::parse_numeric_inputs(&c.estimates)?,
            &super::parse_numeric_inputs(&c.standard_errors)?,
            &options,
        )?;
        let (class, means, sds, lower, upper) = match g.components {
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
        let encode =
            |v: Option<Vec<f64>>| v.map(|v| v.into_iter().map(super::bits).collect::<Vec<_>>());
        result.push(json!({"name":c.name,"component_class":class,"initial_weights":encode(Some(g.initial_weights)),"dirichlet":encode(Some(g.dirichlet)),"means":encode(means),"sds":encode(sds),"lower":encode(lower),"upper":encode(upper),"warnings":g.warnings}));
    }
    Ok(json!({"cases":result}))
}
