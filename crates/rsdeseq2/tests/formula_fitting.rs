//! Independent original-R matrix and offset semantics for complete workflows.
use rsdeseq2::{design::formula_design_for_fitting_from_model_frame, prelude::*};
use serde_json::Value;

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|v| {
            f64::from_bits(
                u64::from_str_radix(v.as_str().unwrap().strip_prefix("bits:").unwrap(), 16)
                    .unwrap(),
            )
        })
        .collect()
}
fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_owned())
        .collect()
}

#[test]
fn fitting_matrices_follow_original_term_order_and_matrix_variable_products() {
    // CI also supplies a fresh original-R oracle from the executing platform.
    // The ordinary workspace run always keeps its committed reference panel.
    let fresh = std::env::var_os("RSDESEQ2_FORMULA_FITTING_REFERENCE")
        .map(|path| std::fs::read_to_string(path).expect("read fresh formula reference"));
    let r: Value = serde_json::from_str(
        fresh
            .as_deref()
            .unwrap_or(include_str!("data/formula_fitting_r462.json")),
    )
    .unwrap();
    let frame = FormulaModelFrame {
        factors: r["factors"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| FormulaFactorColumn {
                name: f["name"].as_str().unwrap().into(),
                sample_levels: strings(&f["sample_levels"]),
                levels: Some(strings(&f["levels"])),
                reference: Some(f["reference"].as_str().unwrap().into()),
            })
            .collect(),
        numeric_covariates: r["numeric"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| FormulaNumericColumn {
                name: v["name"].as_str().unwrap().into(),
                values: numbers(&v["values"]),
            })
            .collect(),
    };
    let cases = r["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 32);
    let bits = |v: &[f64]| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    let mut failures = Vec::new();
    for c in cases {
        let formula = c["formula"].as_str().unwrap();
        let parsed = formula_design_for_fitting_from_model_frame(formula, &frame).unwrap();
        let matrix = parsed.design.standard_design;
        assert_eq!(matrix.n_samples(), r["samples"].as_u64().unwrap() as usize);
        if matrix.n_coefficients() != c["columns"].as_u64().unwrap() as usize {
            failures.push(format!(
                "{} columns: native {}, original {}",
                c["name"],
                matrix.n_coefficients(),
                c["columns"]
            ));
        }
        if bits(matrix.matrix().as_slice()) != bits(&numbers(&c["values"])) {
            failures.push(format!("{} matrix", c["name"]));
        }
        if bits(&parsed.offsets) != bits(&numbers(&c["offsets"])) {
            failures.push(format!("{} offsets", c["name"]));
        }
    }
    assert!(failures.is_empty(), "original R disagreement: {failures:?}");
}
