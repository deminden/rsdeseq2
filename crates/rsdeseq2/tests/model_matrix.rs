use crate::support::formula_reference;
// Frozen independent R matrices, including attributes/order and signed zeros.
#[path = "../tools/validation/evaluate_r_interface/model_matrix.rs"]
mod input;
use rsdeseq2::design::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
use rsdeseq2::errors::ModelMatrixResponseProblem;
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use rsdeseq2::{DeseqError, NumericProfile, numeric_profile};
use sha2::{Digest, Sha256};
use std::assert_matches;

// R 4.6.2 on macOS: qr(outer(seq_len(n) - mean(seq_len(n)), 0:(n-1), `^`))$rank.
// The portable hashes cover row-major binary64 values, serialized big-endian.
const PORTABLE_ORDERED_BASES: [(usize, usize, &str); 7] = [
    (
        32,
        26,
        "9e4f56f1c1043e4c32b607dcf7e7eafafeba7e06172482e0779fa72cb4d07ef2",
    ),
    (
        53,
        33,
        "5be99cfcf0be6ffe94ea65401a6a276a5d35d3ac52952ff9980dc995b92de1d6",
    ),
    (
        64,
        36,
        "40bae3b3220dbd75a9ff34fff76bd1ffe38520e324c60a1070734a6d41c379f2",
    ),
    (
        68,
        37,
        "0a950f85647af4b0589ba5253f74644dfaa5322e2303af41bf9637a2bcd97d97",
    ),
    (
        85,
        40,
        "66edf3766ac98a3eea536c8f819abab08bc7e31754aa042ca0000a1d92545d15",
    ),
    (
        91,
        41,
        "dcb15cd6cb491d4de9c0b984978066dd8d1568203fe943eabd315c756ddfe1a4",
    ),
    (
        95,
        42,
        "5f277f764562db7650fe63400e09383516c8cf1c0303071767edac5e367fc80f",
    ),
];

#[test]
fn model_matrix_materialization_checks_output_capacity() {
    let max_cells = isize::MAX as usize / size_of::<f64>();
    let plan = ModelMatrixPlan::new(max_cells + 1, &[], &[], true).unwrap();
    assert_matches!(
        plan.try_values(),
        Err(DeseqError::InvalidDimensions { context, expected, actual })
            if context == "model matrix values allocation overflow"
                && expected == max_cells
                && actual == max_cells + 1
    );

    let no_columns = ModelMatrixPlan::new(usize::MAX, &[], &[], false).unwrap();
    assert!(no_columns.try_values().unwrap().is_empty());
    assert_eq!(no_columns.column_major_values().next(), None);
}

#[test]
fn formula_term_assembly_matches_independent_r_in_both_storage_orders() {
    let Some(text) = formula_reference!("model_matrix_r462") else {
        return;
    };
    let r: serde_json::Value = serde_json::from_str(text).unwrap();
    let cases = r["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 48);
    let portable = numeric_profile() == NumericProfile::Portable;
    for value in cases {
        let case: input::Case = serde_json::from_value(value.clone()).unwrap();
        let name = value["name"].as_str().unwrap();
        let portable_basis = if portable && name.starts_with("ordered_") {
            PORTABLE_ORDERED_BASES
                .iter()
                .find(|&&(levels, _, _)| levels == case.samples)
        } else {
            None
        };
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
            if let Some(&(levels, rank, hash)) = portable_basis {
                check_portable_ordered_basis(
                    &actual["values"],
                    &value["expected"],
                    levels,
                    rank,
                    hash,
                );
            } else {
                assert_eq!(actual["values"], value["expected"], "{name} {order:?}");
            }
        }
    }
}

fn check_portable_ordered_basis(
    actual: &serde_json::Value,
    reference: &serde_json::Value,
    levels: usize,
    rank: usize,
    expected_hash: &str,
) {
    let bits = |matrix: &serde_json::Value| {
        matrix
            .as_array()
            .unwrap()
            .iter()
            .map(|value| u64::from_str_radix(value.as_str().unwrap(), 16).unwrap())
            .collect::<Vec<_>>()
    };
    let actual = bits(actual);
    let reference = bits(reference);
    assert_eq!(actual.len(), levels * levels);
    assert_eq!(reference.len(), actual.len());
    let dot = |matrix: &[u64], column: usize, other: &[u64], other_column: usize| {
        (0..levels)
            .map(|row| {
                f64::from_bits(matrix[row * levels + column])
                    * f64::from_bits(other[row * levels + other_column])
            })
            .sum::<f64>()
    };
    for column in 0..levels {
        let agreement = dot(&actual, column, &reference, column);
        // Beyond numerical rank, qr.qy supplies an orthogonal completion.
        // The nearly singular R diagonal may orient those columns either way.
        let orientation = if column >= rank && agreement < 0.0 {
            -1.0
        } else {
            1.0
        };
        let bound = if column <= 20 { 2e-10 } else { 5e-5 };
        for row in 0..levels {
            let index = row * levels + column;
            let a = f64::from_bits(actual[index]);
            let b = f64::from_bits(reference[index]);
            assert!(
                (a - orientation * b).abs() <= bound,
                "{levels} levels, row {row}, column {column}: {a} vs {b}"
            );
        }
        for other in 0..=column {
            let inner = dot(&actual, column, &actual, other);
            let expected = if column == 0 {
                levels as f64
            } else {
                f64::from(u8::from(column == other))
            };
            assert!(
                (inner - expected).abs() <= 1e-14,
                "{levels} levels, inner product ({column}, {other}): {inner} vs {expected}"
            );
        }
    }
    let mut hash = Sha256::new();
    for cell in actual {
        hash.update(cell.to_be_bytes());
    }
    let actual_hash = hash
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(actual_hash, expected_hash, "{levels} levels");
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
    assert_matches!(
        ModelMatrixPlan::with_response(3, &[variable], &[vec![1]], true, Some(1))
            .err()
            .unwrap(),
        DeseqError::InvalidModelMatrixResponse {
            index: 1,
            variables: 1,
            problem: ModelMatrixResponseProblem::OutOfBounds,
        }
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
    assert_matches!(
        ModelMatrixPlan::with_response(3, &[valid_factor], &[vec![1]], true, Some(0))
            .err()
            .unwrap(),
        DeseqError::InvalidModelMatrixResponse {
            index: 0,
            variables: 1,
            problem: ModelMatrixResponseProblem::NonNumeric,
        }
    );
    let empty = MatrixView::new(0, usize::MAX, &[], MatrixOrder::RowMajor).unwrap();
    for values in [matrix, empty] {
        let custom = ModelMatrixVariable::Factor {
            name: "f",
            codes: &codes,
            levels: &levels,
            coding: FactorCoding::Matrix {
                values,
                names: None,
            },
        };
        let error = ModelMatrixPlan::new(3, &[custom], &[vec![1]], true)
            .err()
            .unwrap();
        assert_eq!(
            error.to_string(),
            "invalid options: explicit factor contrasts must be a finite level-by-column matrix"
        );
    }
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
    let Some(text) = formula_reference!("model_matrix_r462") else {
        return;
    };
    let reference: serde_json::Value = serde_json::from_str(text).unwrap();
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
