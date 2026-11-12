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
    #[serde(default)]
    factor_expansions: Vec<super::factor_expansion::Case>,
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
    #[serde(default)]
    details: bool,
    #[serde(default)]
    frame_names: Vec<String>,
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
            Ok(result) => {
                let bits = |v: &[f64]| v.iter().map(|v| format!("{:016x}",v.to_bits())).collect::<Vec<_>>();
                let mut output = json!({"name":case.name,"ok":true,
                    "names":result.design.standard_design.coefficient_names(),
                    "values":bits(result.design.standard_design.matrix().as_slice()),
                    "offsets":bits(&result.offsets)});
                if case.details {
                    output["expanded_names"] = json!(result.design.expanded_design.coefficient_names());
                    output["expanded_values"] = json!(bits(result.design.expanded_design.matrix().as_slice()));
                    output["groups"] = json!(result.design.coefficient_groups);
                    output["frames"] = json!(case.frame_names.iter().map(|name| {
                        result.model_frame.factors.iter().find(|f|f.name==*name).map(|f|json!({
                            "name":f.name,"levels":f.levels,"reference":f.reference,"sample_levels":f.sample_levels,
                            "ordered":result.design.standard_design.origin().ordered_factors().contains(&f.name),
                        })).unwrap_or(Value::Null)
                    }).collect::<Vec<_>>());
                }
                output
            },
            Err(_) => json!({"name":case.name,"ok":false}),
        }
    }).collect::<Vec<_>>();
    let mut output = json!({"cases":cases});
    if !input.factor_expansions.is_empty() {
        output["expansions"] = super::factor_expansion::evaluate(&input.factor_expansions);
    }
    output
}
