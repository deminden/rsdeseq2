use rsdeseq2::prelude::*;
use serde::Deserialize;
use serde_json::{Value as Json, json};
use std::error::Error;

#[derive(Deserialize)]
pub struct Input {
    pub frozen_local_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    name: String,
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    size_factors: Option<Vec<String>>,
    normalization_factors: Option<Vec<String>>,
    model: Json,
}
fn numbers(value: &Json) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().parse().unwrap())
        .collect()
}
pub fn evaluate(input: Input) -> Result<Json, Box<dyn Error>> {
    let mut cases = vec![];
    for case in input.frozen_local_cases {
        let counts = CountMatrix::from_row_major_u32(case.genes, case.samples, case.counts)?;
        let parse = |v: Vec<String>| {
            v.iter()
                .map(|v| v.parse::<f64>())
                .collect::<Result<Vec<_>, _>>()
        };
        let (normalized, inverse) = if let Some(nf) = case.normalization_factors {
            let nf = RowMajorMatrix::from_row_major(case.genes, case.samples, parse(nf)?)?;
            (
                normalized_counts_with_factors(&counts, &nf)?,
                local_vst_inverse_size_factor_mean_from_normalization_factors(&nf)?,
            )
        } else {
            let sf = parse(case.size_factors.ok_or("size factors required")?)?;
            (
                normalized_counts(&counts, &sf)?,
                local_vst_inverse_size_factor_mean(&sf)?,
            )
        };
        let model = &case.model;
        let curve = if model["type"] == "saved" {
            FrozenLocalDispersion::Saved(SavedLocalDispersionFit::new(
                numbers(&model["knots"]),
                numbers(&model["residual_values"]),
                numbers(&model["residual_slopes"]),
                numbers(&model["polynomial"])
                    .try_into()
                    .map_err(|_| "four global coefficients required")?,
            )?)
        } else {
            FrozenLocalDispersion::Expression {
                argument: model["argument"].as_str().unwrap().into(),
                expression: super::formula::expression(&model["expression"]),
                constants: model["constants"]
                    .as_object()
                    .map(|v| {
                        v.iter()
                            .map(|(name, value)| (name.clone(), super::formula::value(value)))
                            .collect()
                    })
                    .unwrap_or_default(),
            }
        };
        let (mut grid, dispersions) = curve.grid(&normalized)?;
        grid.extend(dispersions);
        let transformed = vst_frozen_local(&normalized, &curve, inverse)?;
        cases.push(json!({"name":case.name,"curve":grid.into_iter().map(super::bits).collect::<Vec<_>>(),
            "transformed":transformed.as_slice().iter().copied().map(super::bits).collect::<Vec<_>>()}));
    }
    Ok(json!({"cases":cases}))
}
