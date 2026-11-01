#[macro_use]
#[path = "support/formula_reference.rs"]
mod reference;
// Frozen independent R matrices, including attributes/order and signed zeros.
#[path = "../examples/evaluate_r_interface/model_matrix.rs"]
mod input;
use rsdeseq2::design::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};

#[test]
fn formula_term_assembly_matches_independent_r_in_both_storage_orders() {
    let r: serde_json::Value =
        serde_json::from_str(formula_reference!("model_matrix_r462")).unwrap();
    let cases = r["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 48);
    for value in cases {
        let case: input::Case = serde_json::from_value(value.clone()).unwrap();
        for order in [MatrixOrder::RowMajor, MatrixOrder::ColumnMajor] {
            let actual = input::evaluate(&case, order);
            assert_eq!(
                actual["names"], value["expected_names"],
                "{} {order:?}",
                value["name"]
            );
            assert_eq!(
                actual["assignments"], value["assignments"],
                "{} {order:?}",
                value["name"]
            );
            assert_eq!(
                actual["values"], value["expected"],
                "{} {order:?}",
                value["name"]
            );
        }
    }
}
#[test]
fn model_matrix_rejects_invalid_shapes_codes_names_and_contrasts() {
    let numeric = [1.0, 2.0, 3.0];
    let matrix = MatrixView::new(3, 1, &numeric, MatrixOrder::RowMajor).unwrap();
    let variable = ModelMatrixVariable::Numeric {
        name: "x",
        values: matrix,
        names: None,
    };
    for terms in [vec![vec![3]], vec![vec![]], vec![vec![0]]] {
        assert!(ModelMatrixPlan::new(3, &[variable], &terms, true).is_err());
    }
    assert!(ModelMatrixPlan::new(2, &[variable], &[vec![1]], true).is_err());
    assert!(
        ModelMatrixPlan::with_response(3, &[variable], &[vec![1]], true, Some(1))
            .err()
            .unwrap()
            .to_string()
            .contains("response")
    );
    let levels = vec!["a".into(), "b".into()];
    for codes in [[0, 1, 2], [1, 3, 2], [i32::MIN, 1, 2]] {
        let factor = ModelMatrixVariable::Factor {
            name: "f",
            codes: &codes,
            levels: &levels,
            coding: FactorCoding::Treatment { base: 0 },
        };
        assert!(ModelMatrixPlan::new(3, &[factor], &[vec![1]], true).is_err());
    }
    let codes = [1, 2, 1];
    let factor = ModelMatrixVariable::Factor {
        name: "f",
        codes: &codes,
        levels: &levels,
        coding: FactorCoding::Treatment { base: 2 },
    };
    assert!(ModelMatrixPlan::new(3, &[factor], &[vec![1]], true).is_err());
    let valid_factor = ModelMatrixVariable::Factor {
        name: "f",
        codes: &codes,
        levels: &levels,
        coding: FactorCoding::Treatment { base: 0 },
    };
    assert!(
        ModelMatrixPlan::with_response(3, &[valid_factor], &[vec![1]], true, Some(0))
            .err()
            .unwrap()
            .to_string()
            .contains("response")
    );
    let custom = ModelMatrixVariable::Factor {
        name: "f",
        codes: &codes,
        levels: &levels,
        coding: FactorCoding::Matrix {
            values: matrix,
            names: None,
        },
    };
    assert!(ModelMatrixPlan::new(3, &[custom], &[vec![1]], true).is_err());
}

#[test]
fn expanded_list_contrasts_match_independent_r_group_means_and_rejections() {
    #[derive(serde::Deserialize)]
    struct Case {
        name: String,
        coefficients: usize,
        positive: Vec<usize>,
        negative: Vec<usize>,
        expected: Option<Vec<String>>,
    }
    let reference: serde_json::Value =
        serde_json::from_str(formula_reference!("model_matrix_r462")).unwrap();
    let cases = reference["average_contrasts"].as_array().unwrap();
    assert_eq!(cases.len(), 10);
    for case in cases {
        let case: Case = serde_json::from_value(case.clone()).unwrap();
        let actual = rsdeseq2::contrasts::averaged_list_contrast(
            case.coefficients,
            &case.positive,
            &case.negative,
        );
        if let Some(expected) = case.expected {
            assert_eq!(
                actual
                    .unwrap()
                    .iter()
                    .map(|v| format!("{:016x}", v.to_bits()))
                    .collect::<Vec<_>>(),
                expected,
                "{}",
                case.name
            );
        } else {
            assert!(actual.is_err(), "{}", case.name)
        }
    }
}
