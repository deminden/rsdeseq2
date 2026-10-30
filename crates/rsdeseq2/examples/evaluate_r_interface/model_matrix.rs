//! JSON transport for direct interface conformance; numerical work stays in core.
use rsdeseq2::design::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
pub struct Case {
    pub samples: usize,
    pub variables: Vec<Variable>,
    pub terms: Vec<Vec<i32>>,
    pub intercept: bool,
    #[serde(default)]
    pub response: usize,
}
#[derive(Deserialize)]
pub struct Variable {
    pub name: String,
    #[serde(rename = "type")]
    pub kind: String,
    #[serde(default)]
    pub columns: usize,
    #[serde(default)]
    pub values: Vec<String>,
    pub column_names: Option<Vec<String>>,
    #[serde(default)]
    pub codes: Vec<i32>,
    #[serde(default)]
    pub levels: Vec<String>,
    #[serde(default)]
    pub coding: String,
    pub contrast_columns: Option<usize>,
    pub contrast_values: Option<Vec<String>>,
    pub contrast_names: Option<Vec<String>>,
}
fn numbers(values: &[String]) -> Vec<f64> {
    values
        .iter()
        .map(|s| f64::from_bits(u64::from_str_radix(s, 16).unwrap()))
        .collect()
}
fn storage(values: &[f64], rows: usize, cols: usize, order: MatrixOrder) -> Vec<f64> {
    match order {
        MatrixOrder::RowMajor => values.to_vec(),
        MatrixOrder::ColumnMajor => (0..cols)
            .flat_map(|c| (0..rows).map(move |r| values[r * cols + c]))
            .collect(),
    }
}
pub fn evaluate(case: &Case, order: MatrixOrder) -> Value {
    let numeric = case
        .variables
        .iter()
        .map(|v| storage(&numbers(&v.values), case.samples, v.columns, order))
        .collect::<Vec<_>>();
    let custom = case
        .variables
        .iter()
        .map(|v| {
            storage(
                &numbers(v.contrast_values.as_deref().unwrap_or(&[])),
                v.levels.len(),
                v.contrast_columns.unwrap_or(0),
                order,
            )
        })
        .collect::<Vec<_>>();
    let mut variables = Vec::new();
    for (i, v) in case.variables.iter().enumerate() {
        variables.push(if v.kind == "numeric" {
            ModelMatrixVariable::Numeric {
                name: &v.name,
                values: MatrixView::new(case.samples, v.columns, &numeric[i], order).unwrap(),
                names: v.column_names.as_deref(),
            }
        } else {
            let coding = match v.coding.as_str() {
                "contr.treatment" => FactorCoding::Treatment { base: 0 },
                "contr.SAS" => FactorCoding::Treatment {
                    base: v.levels.len() - 1,
                },
                "contr.sum" => FactorCoding::Sum,
                "contr.helmert" => FactorCoding::Helmert,
                "contr.poly" => FactorCoding::Polynomial,
                "matrix" => FactorCoding::Matrix {
                    values: MatrixView::new(
                        v.levels.len(),
                        v.contrast_columns.unwrap(),
                        &custom[i],
                        order,
                    )
                    .unwrap(),
                    names: v.contrast_names.as_deref(),
                },
                _ => panic!("unknown reference coding"),
            };
            ModelMatrixVariable::Factor {
                name: &v.name,
                codes: &v.codes,
                levels: &v.levels,
                coding,
            }
        });
    }
    let plan = ModelMatrixPlan::with_response(
        case.samples,
        &variables,
        &case.terms,
        case.intercept,
        case.response.checked_sub(1),
    )
    .unwrap();
    json!({"names":plan.column_names(),"assignments":plan.assignments(),"values":plan.values().iter().map(|v|format!("{:016x}",v.to_bits())).collect::<Vec<_>>()})
}
