//! Test-only actual string-parser conformance panel.
use rsdeseq2::design::{
    ExpandedFactorSpec, ExpandedNumericSpec, expanded_formula_design_with_offsets,
};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
pub struct Input {
    pub native_formula_cases: Vec<Case>,
    numeric: Vec<Numeric>,
    factors: Vec<Factor>,
}
#[derive(Deserialize)]
struct Numeric {
    name: String,
    values: Vec<String>,
}
#[derive(Deserialize)]
struct Factor {
    name: String,
    sample_levels: Vec<String>,
    levels: Vec<String>,
    reference: String,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    formula: String,
}

pub fn evaluate(input: &Input) -> Value {
    let numeric = input
        .numeric
        .iter()
        .map(|v| {
            v.values
                .iter()
                .map(|s| f64::from_bits(u64::from_str_radix(s, 16).unwrap()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let numeric = input
        .numeric
        .iter()
        .zip(&numeric)
        .map(|(v, values)| ExpandedNumericSpec {
            name: &v.name,
            values,
        })
        .collect::<Vec<_>>();
    let factors = input
        .factors
        .iter()
        .map(|v| ExpandedFactorSpec {
            factor: &v.name,
            sample_levels: &v.sample_levels,
            levels: Some(&v.levels),
            reference: &v.reference,
        })
        .collect::<Vec<_>>();
    let cases = input.native_formula_cases.iter().map(|case| {
        match expanded_formula_design_with_offsets(&case.formula, &factors, &numeric) {
            Ok(result) => json!({"name":case.name,"ok":true,
                "names":result.design.standard_design.coefficient_names(),
                "values":result.design.standard_design.matrix().as_slice().iter().map(|v|format!("{:016x}",v.to_bits())).collect::<Vec<_>>(),
                "offsets":result.offsets.iter().map(|v|format!("{:016x}",v.to_bits())).collect::<Vec<_>>() }),
            Err(_) => json!({"name":case.name,"ok":false}),
        }
    }).collect::<Vec<_>>();
    json!({"cases":cases})
}
