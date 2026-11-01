#[macro_use]
#[path = "support/formula_reference.rs"]
mod reference;
use rsdeseq2::design::numeric_transform::{
    PolynomialCoefficients, ScalePolicy, integer_covariate, polynomial_basis, scale_covariates,
};
use rsdeseq2::design::{ExpandedNumericSpec, expanded_formula_design};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use serde::Deserialize;

#[derive(Deserialize)]
struct Reference {
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Coefficients {
    alpha: Vec<String>,
    norm2: Vec<String>,
}
#[derive(Deserialize)]
struct Policy {
    mode: String,
    values: Option<Vec<String>>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    operation: String,
    input: Vec<String>,
    rows: usize,
    columns: usize,
    degree: Option<usize>,
    raw: Option<bool>,
    supplied_coefficients: Option<Coefficients>,
    center_policy: Option<Policy>,
    scale_policy: Option<Policy>,
    expected: Option<Vec<String>>,
    coefficients: Option<Coefficients>,
    center: Option<Vec<String>>,
    scale: Option<Vec<String>>,
}
fn floats(hex: &[String]) -> Vec<f64> {
    hex.iter()
        .map(|v| f64::from_bits(u64::from_str_radix(v, 16).unwrap()))
        .collect()
}
fn bits(values: &[f64], expected: &[String], name: &str) {
    assert_eq!(values.len(), expected.len(), "{name}: shape");
    for (i, (v, e)) in values.iter().zip(expected).enumerate() {
        assert_eq!(
            format!("{:016x}", v.to_bits()),
            *e,
            "{name}: element {i}: {v:?}"
        );
    }
}
fn controls(policy: &Policy) -> (u8, Vec<f64>) {
    (
        match policy.mode.as_str() {
            "automatic" => 0,
            "disabled" => 1,
            "supplied" => 2,
            _ => panic!(),
        },
        policy
            .values
            .as_ref()
            .map(|v| floats(v))
            .unwrap_or_default(),
    )
}
fn policy(control: &(u8, Vec<f64>)) -> ScalePolicy<'_> {
    match control.0 {
        0 => ScalePolicy::Automatic,
        1 => ScalePolicy::Disabled,
        _ => ScalePolicy::Supplied(&control.1),
    }
}

// Exercise the production parser as well as the reusable numerical leaf.
fn check_formula(case: &Case, input: &[f64]) {
    let control = |p: &Policy| match p.mode.as_str() {
        "automatic" => "TRUE".to_string(),
        "disabled" => "FALSE".to_string(),
        "supplied" => format!("{:?}", floats(p.values.as_ref().unwrap())[0]),
        _ => panic!(),
    };
    let formula = match case.operation.as_str() {
        "poly" if case.supplied_coefficients.is_none() && case.degree.unwrap() <= 16 => format!(
            "~poly(dose,{},raw={})",
            case.degree.unwrap(),
            if case.raw.unwrap() { "TRUE" } else { "FALSE" }
        ),
        "scale" if case.columns == 1 => format!(
            "~scale(dose,center={},scale={})",
            control(case.center_policy.as_ref().unwrap()),
            control(case.scale_policy.as_ref().unwrap())
        ),
        "integer" => "~as.integer(dose)".into(),
        _ => return,
    };
    let result = expanded_formula_design(
        &formula,
        &[],
        &[ExpandedNumericSpec {
            name: "dose",
            values: input,
        }],
    );
    if let Some(expected) = &case.expected {
        let design = result.unwrap_or_else(|e| panic!("{} parser: {e}", case.name));
        let matrix = design.standard_design.matrix();
        let values = matrix
            .as_slice()
            .chunks_exact(matrix.n_cols())
            .flat_map(|r| r[1..].iter().copied())
            .collect::<Vec<_>>();
        bits(&values, expected, &format!("{} parser", case.name));
    } else {
        assert!(result.is_err(), "{} parser should reject", case.name);
    }
}

#[test]
fn formula_transforms_match_independent_r_bits_and_invalid_output_decisions() {
    let reference: Reference =
        serde_json::from_str(formula_reference!("numeric_transform_r462")).unwrap();
    for case in reference.cases {
        let input = floats(&case.input);
        check_formula(&case, &input);
        match case.operation.as_str() {
            "poly" => {
                let supplied =
                    case.supplied_coefficients
                        .as_ref()
                        .map(|c| PolynomialCoefficients {
                            alpha: floats(&c.alpha),
                            norm2: floats(&c.norm2),
                        });
                let result = polynomial_basis(
                    &input,
                    case.degree.unwrap(),
                    case.raw.unwrap(),
                    supplied.as_ref(),
                );
                if let Some(expected) = case.expected {
                    let observed = result.unwrap_or_else(|e| panic!("{}: {e}", case.name));
                    bits(observed.values.as_slice(), &expected, &case.name);
                    match (observed.coefficients, case.coefficients) {
                        (Some(o), Some(e)) => {
                            bits(&o.alpha, &e.alpha, &format!("{} alpha", case.name));
                            bits(&o.norm2, &e.norm2, &format!("{} norm2", case.name));
                        }
                        (None, None) => {}
                        _ => panic!("{} coefficient presence", case.name),
                    }
                } else {
                    assert!(result.is_err(), "{} should reject", case.name);
                }
            }
            "scale" => {
                let center = controls(case.center_policy.as_ref().unwrap());
                let scale = controls(case.scale_policy.as_ref().unwrap());
                for order in [MatrixOrder::RowMajor, MatrixOrder::ColumnMajor] {
                    let storage = match order {
                        MatrixOrder::RowMajor => input.clone(),
                        MatrixOrder::ColumnMajor => (0..case.columns)
                            .flat_map(|c| (0..case.rows).map(move |r| (r, c)))
                            .map(|(r, c)| input[r * case.columns + c])
                            .collect(),
                    };
                    let result = scale_covariates(
                        MatrixView::new(case.rows, case.columns, &storage, order).unwrap(),
                        policy(&center),
                        policy(&scale),
                    );
                    if let Some(expected) = &case.expected {
                        let observed = result.unwrap_or_else(|e| panic!("{}: {e}", case.name));
                        bits(observed.values.as_slice(), expected, &case.name);
                        for (actual, expected) in [
                            (observed.center, &case.center),
                            (observed.scale, &case.scale),
                        ] {
                            match (actual, expected) {
                                (Some(a), Some(e)) => bits(&a, e, &case.name),
                                (None, None) => {}
                                _ => panic!("{} scale metadata", case.name),
                            }
                        }
                    } else {
                        assert!(result.is_err(), "{} should reject", case.name);
                    }
                }
            }
            "integer" => match &case.expected {
                Some(expected) => bits(
                    &[integer_covariate(input[0]).unwrap()],
                    expected,
                    &case.name,
                ),
                None => assert!(integer_covariate(input[0]).is_err(), "{}", case.name),
            },
            _ => panic!("unknown transform"),
        }
    }
}

#[test]
fn polynomial_and_scaling_validate_controls_before_using_buffers() {
    assert!(polynomial_basis(&[1.0, 2.0], usize::MAX, false, None).is_err());
    assert!(polynomial_basis(&[1.0, 2.0], usize::MAX - 1, false, None).is_err());
    for coefficients in [
        PolynomialCoefficients {
            alpha: vec![],
            norm2: vec![1.0, 2.0, 1.0],
        },
        PolynomialCoefficients {
            alpha: vec![1.0],
            norm2: vec![1.0, 2.0],
        },
        PolynomialCoefficients {
            alpha: vec![1.0],
            norm2: vec![1.0, 2.0, -1.0],
        },
        PolynomialCoefficients {
            alpha: vec![f64::NAN],
            norm2: vec![1.0, 2.0, 1.0],
        },
    ] {
        assert!(polynomial_basis(&[1.0, 2.0], 1, false, Some(&coefficients)).is_err());
    }
    let input = MatrixView::new(2, 1, &[1.0, 2.0], MatrixOrder::ColumnMajor).unwrap();
    for supplied in [&[][..], &[1.0, 2.0][..], &[f64::INFINITY][..]] {
        assert!(
            scale_covariates(
                input,
                ScalePolicy::Supplied(supplied),
                ScalePolicy::Disabled
            )
            .is_err()
        );
        assert!(
            scale_covariates(
                input,
                ScalePolicy::Disabled,
                ScalePolicy::Supplied(supplied)
            )
            .is_err()
        );
    }
}
