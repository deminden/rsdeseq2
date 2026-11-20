use crate::support::formula_reference;
use rsdeseq2::DeseqError;
use rsdeseq2::design::expression::{
    Argument, Data, EvaluationOptions, Expression, Value, evaluate,
};
use rsdeseq2::errors::FormulaResource;
use serde_json::Value as Json;
use std::assert_matches;

#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/formula.rs"]
mod formula;
use formula::{expression, value};

fn equal_value(actual: &Value, expected: &Value, path: &str) {
    assert_eq!(actual.dimensions, expected.dimensions, "{path}: dimensions");
    assert_eq!(actual.names, expected.names, "{path}: names");
    assert_eq!(actual.dimnames, expected.dimnames, "{path}: dimnames");
    assert_eq!(
        actual.dimname_names, expected.dimname_names,
        "{path}: label names"
    );
    match (&actual.data, &expected.data) {
        (
            Data::Number {
                values: a,
                integer: ai,
            },
            Data::Number {
                values: b,
                integer: bi,
            },
        ) => {
            assert_eq!(ai, bi, "{path}: integer type");
            assert_eq!(a.len(), b.len(), "{path}: length");
            for (i, (a, b)) in a.iter().zip(b).enumerate() {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "{path}: element {i}, {a:?} vs {b:?}"
                )
            }
        }
        (Data::List(a), Data::List(b)) => {
            assert_eq!(a.len(), b.len(), "{path}: list length");
            for (i, (a, b)) in a.iter().zip(b).enumerate() {
                equal_value(a, b, &format!("{path}[[{i}]]"))
            }
        }
        (a, b) => assert_eq!(a, b, "{path}: data"),
    }
    let mut an = actual.attributes.iter().map(|(n, _)| n).collect::<Vec<_>>();
    let mut en = expected
        .attributes
        .iter()
        .map(|(n, _)| n)
        .collect::<Vec<_>>();
    an.sort();
    en.sort();
    assert_eq!(an, en, "{path}: attributes");
    for (name, expected) in &expected.attributes {
        equal_value(
            actual.attribute(name).unwrap(),
            expected,
            &format!("{path}@{name}"),
        )
    }
}
#[test]
fn built_in_expressions_match_independent_r_values_and_metadata() {
    compare_reference(formula_reference!("formula_expression_r462"));
}
#[test]
fn extended_builtin_arguments_and_matrix_polynomials_match_r() {
    compare_reference(formula_reference!("formula_expression_extended_r462"));
    compare_reference(include_str!("data/formula_elementary_boundaries_r462.json"));
}

#[test]
fn sum_preserves_r_missing_precedence_across_nan_orders() {
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    let profile = rsdeseq2::numeric_profile();
    let linux = matches!(
        profile,
        rsdeseq2::NumericProfile::RLinuxX86 | rsdeseq2::NumericProfile::RLinuxArm64
    );
    let generated_nan = match profile {
        rsdeseq2::NumericProfile::RLinuxX86 => 0xfff8_0000_0000_0000,
        rsdeseq2::NumericProfile::RLinuxArm64 => 0x7fff_ffff_ffff_ffff,
        _ => f64::NAN.to_bits(),
    };
    let options = EvaluationOptions::default();
    for (values, expected) in [
        (vec![na, f64::NAN], 0x7ff8_0000_0000_07a2),
        (
            vec![f64::NAN, na],
            if linux {
                f64::NAN.to_bits()
            } else {
                0x7ff8_0000_0000_07a2
            },
        ),
        (
            vec![f64::INFINITY, f64::NEG_INFINITY, na],
            if linux {
                generated_nan
            } else {
                0x7ff8_0000_0000_07a2
            },
        ),
    ] {
        let input = Expression::Literal(Value::numbers(values));
        let result = evaluate(&call("sum", vec![input]), &[], &options).unwrap();
        let Data::Number { values, .. } = result.data else {
            panic!("sum must return a number");
        };
        assert_eq!(values[0].to_bits(), expected);
    }
    let input = Expression::Literal(Value::numbers(vec![f64::INFINITY, f64::NEG_INFINITY]));
    let result = evaluate(&call("sum", vec![input]), &[], &options).unwrap();
    let Data::Number { values, .. } = result.data else {
        panic!("sum must return a number");
    };
    assert_eq!(values[0].to_bits(), generated_nan);
    let payload = f64::from_bits(0x7ff8_0000_0000_1234);
    for (input, expected) in [
        (vec![payload, f64::NAN], payload.to_bits()),
        (vec![f64::NAN, payload], f64::NAN.to_bits()),
        (
            vec![f64::INFINITY, f64::NEG_INFINITY, payload],
            generated_nan,
        ),
    ] {
        let input = Expression::Literal(Value::numbers(input));
        let result = evaluate(&call("sum", vec![input]), &[], &options).unwrap();
        let Data::Number { values, .. } = result.data else {
            panic!("sum must return a number");
        };
        assert_eq!(values[0].to_bits(), expected);
    }
}

#[test]
fn sum_preserves_r_integer_promotion_and_missing_type() {
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    let max = f64::from(i32::MAX);
    let integer = |values| {
        Expression::Literal(Value::new(Data::Number {
            values,
            integer: true,
        }))
    };
    let number = |x| Expression::Literal(Value::scalar(x));
    let logical_na = Expression::Literal(Value::new(Data::Logical(vec![None])));
    let options = EvaluationOptions::default();

    // Type and payload observed in original R 4.6.2 on macOS and Linux arm64.
    for (arguments, expected, integer_type) in [
        (vec![number(f64::NAN), logical_na], na, false),
        (vec![integer(vec![na]), number(f64::NAN)], na, false),
        (
            vec![integer(vec![max]), integer(vec![1.0]), integer(vec![-1.0])],
            max,
            false,
        ),
        (vec![integer(vec![max, 1.0, -1.0])], max, true),
        (
            vec![integer(vec![max]), integer(vec![1.0]), integer(vec![na])],
            na,
            false,
        ),
        (
            vec![integer(vec![na]), integer(vec![max]), integer(vec![1.0])],
            na,
            true,
        ),
    ] {
        let result = evaluate(&call("sum", arguments), &[], &options).unwrap();
        let Data::Number { values, integer } = result.data else {
            panic!("sum must return a number");
        };
        assert_eq!(values[0].to_bits(), expected.to_bits());
        assert_eq!(integer, integer_type);
    }
}

#[test]
fn mean_preserves_r_missing_and_nan_order() {
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    let options = EvaluationOptions::default();
    let linux_x86 = rsdeseq2::numeric_profile() == rsdeseq2::NumericProfile::RLinuxX86;
    for (values, expected) in [
        (vec![na, f64::NAN], 0x7ff8_0000_0000_07a2),
        (
            vec![f64::NAN, na],
            if linux_x86 {
                0x7ff8_0000_0000_07a2
            } else {
                f64::NAN.to_bits()
            },
        ),
    ] {
        let input = Expression::Literal(Value::numbers(values));
        let result = evaluate(&call("mean", vec![input]), &[], &options).unwrap();
        let Data::Number { values, .. } = result.data else {
            panic!("mean must return a number");
        };
        assert_eq!(values[0].to_bits(), expected);
    }
}
#[test]
fn linux_arm64_scalar_reductions_follow_r_overflow_and_nan_bits() {
    if rsdeseq2::numeric_profile() != rsdeseq2::NumericProfile::RLinuxArm64 {
        return;
    }
    let literal = |values| Expression::Literal(Value::numbers(values));
    let options = EvaluationOptions::default();
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    // Original R 4.6.2 on Linux arm64, with a 113-bit long double.
    for (name, inputs, expected) in [
        (
            "sum",
            vec![literal(vec![f64::MAX, 2.0_f64.powi(961)])],
            f64::INFINITY,
        ),
        ("sum", vec![literal(vec![1e308, 1e308, -1e308])], 1e308),
        ("mean", vec![literal(vec![1e308, 1e308])], 1e308),
        (
            "sum",
            vec![
                literal(vec![f64::INFINITY]),
                literal(vec![f64::NEG_INFINITY]),
                literal(vec![na]),
            ],
            f64::NAN,
        ),
        (
            "mean",
            vec![literal(vec![f64::INFINITY, f64::NEG_INFINITY, na])],
            f64::from_bits(0x7fff_ffff_ffff_ffff),
        ),
    ] {
        let result = evaluate(&call(name, inputs), &[], &options).unwrap();
        let Data::Number { values, .. } = result.data else {
            panic!("{name} must return a number");
        };
        assert_eq!(values[0].to_bits(), expected.to_bits());
    }
}

fn compare_reference(contents: &str) {
    let reference: Json = serde_json::from_str(contents).unwrap();
    let variables = reference["variables"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(n, v)| (n.clone(), value(v)))
        .collect::<Vec<_>>();
    let options = EvaluationOptions {
        collation: reference["metadata"]["collation"].as_str().unwrap().into(),
    };
    let mut failures = Vec::new();
    for case in reference["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let result = evaluate(&expression(&case["expression"]), &variables, &options)
            .map_err(|e| e.to_string());
        let comparison = std::panic::catch_unwind(|| {
            if case["expected"].is_null() {
                assert!(result.is_err(), "{name}: expected rejection")
            } else {
                equal_value(
                    &result.unwrap_or_else(|e| panic!("{name}: {e}")),
                    &value(&case["expected"]),
                    name,
                )
            }
        });
        if comparison.is_err() {
            failures.push(name.to_owned())
        }
    }
    assert!(
        failures.is_empty(),
        "independent expression failures: {failures:?}"
    );
}

#[test]
fn primitive_wire_preserves_reference_bits_and_rejects_truncation() {
    use rsdeseq2::design::expression::wire;
    let reference: Json =
        serde_json::from_str(formula_reference!("formula_expression_r462")).unwrap();
    for (name, input) in reference["variables"].as_object().unwrap() {
        let original = value(input);
        let bytes = wire::encode(&original).unwrap();
        equal_value(&wire::decode(&bytes).unwrap(), &original, name);
    }
    let nested = value(&reference["variables"]["coefficients"]);
    let bytes = wire::encode(&nested).unwrap();
    for end in 0..bytes.len() {
        assert!(
            wire::decode(&bytes[..end]).is_err(),
            "accepted truncation at {end}"
        );
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(wire::decode(&trailing).is_err());
    let mut oversized = wire::encode(&Value::scalar(1.0)).unwrap();
    oversized[5..13].copy_from_slice(&(usize::MAX as u64).to_le_bytes());
    assert_matches!(
        wire::decode(&oversized).unwrap_err(),
        DeseqError::InvalidFormulaTransportLength {
            count: usize::MAX, element_bytes: 8, remaining_bytes,
        } if remaining_bytes == oversized.len() - 13
    );
    let huge = call(
        "seq_len",
        vec![Expression::Literal(Value::scalar(i32::MAX as f64))],
    );
    assert_matches!(
        evaluate(&huge, &[], &EvaluationOptions::default()).unwrap_err(),
        DeseqError::FormulaLimitExceeded {
            resource: FormulaResource::OutputElements, limit, observed: Some(count),
        } if limit == wire::MAX_BYTES / size_of::<f64>() && count == i32::MAX as usize
    );
}

fn call(name: &str, inputs: Vec<Expression>) -> Expression {
    Expression::Call {
        function: name.into(),
        arguments: inputs
            .into_iter()
            .map(|expression| Argument {
                name: None,
                expression,
            })
            .collect(),
    }
}

#[test]
fn rejects_custom_calls_even_in_unused_branches_and_validates_owned_inputs() {
    let options = EvaluationOptions::default();
    let custom = call(
        "custom_function",
        vec![Expression::Symbol("unknown".into())],
    );
    let expression = call(
        "ifelse",
        vec![
            Expression::Literal(Value::scalar(1.0)),
            Expression::Literal(Value::scalar(2.0)),
            custom,
        ],
    );
    assert_matches!(
        evaluate(&expression, &[], &options).unwrap_err(),
        DeseqError::UnsupportedFormulaFunction { function, expected_namespace: None }
            if function == "custom_function"
    );
    for (name, namespace) in [
        ("base::poly", Some("stats::")),
        ("stats::log", Some("base::")),
        ("other::log", None),
    ] {
        assert_matches!(
            evaluate(&call(name, vec![]), &[], &options).unwrap_err(),
            DeseqError::UnsupportedFormulaFunction { function, expected_namespace }
                if function == name && expected_namespace == namespace
        );
    }
    let mut malformed = Vec::new();
    let mut v = Value::numbers(vec![1.0]);
    v.dimensions = Some([usize::MAX, 2]);
    malformed.push(v);
    let mut v = Value::numbers(vec![1.0]);
    v.dimnames[0] = Some(vec!["a".into()]);
    malformed.push(v);
    let mut v = Value::numbers(vec![1.0]);
    v.names = Some(vec![]);
    malformed.push(v);
    let mut v = Value::numbers(vec![1.0]);
    v.dimensions = Some([1, 1]);
    v.dimnames[1] = Some(vec![]);
    malformed.push(v);
    malformed.push(Value::new(Data::Number {
        values: vec![1.5],
        integer: true,
    }));
    let mut v = Value::new(Data::Number {
        values: vec![2.0],
        integer: true,
    });
    v.set_attribute("class", Value::new(Data::Text(vec![Some("factor".into())])));
    v.set_attribute("levels", Value::new(Data::Text(vec![Some("only".into())])));
    malformed.push(v);
    let mut v = Value::scalar(1.0);
    v.attributes = vec![
        ("a".into(), Value::scalar(1.0)),
        ("a".into(), Value::scalar(2.0)),
    ];
    malformed.push(v);
    for v in malformed {
        assert!(evaluate(&Expression::Literal(v.clone()), &[], &options).is_err());
        assert!(
            evaluate(
                &Expression::Symbol("x".into()),
                &[("x".into(), v)],
                &options
            )
            .is_err()
        );
    }
    assert!(
        evaluate(
            &Expression::Symbol("x".into()),
            &[
                ("x".into(), Value::scalar(1.0)),
                ("x".into(), Value::scalar(2.0))
            ],
            &options
        )
        .is_err()
    );
    let mut nested = Expression::Literal(Value::scalar(1.0));
    for _ in 0..128 {
        nested = call("(", vec![nested]);
    }
    assert_eq!(
        evaluate(&nested, &[], &options).unwrap(),
        Value::scalar(1.0)
    );
    let mut nested = Expression::Literal(Value::scalar(1.0));
    for _ in 0..130 {
        nested = call("I", vec![nested]);
    }
    let error = evaluate(&nested, &[], &options).unwrap_err();
    assert_eq!(
        error.to_string(),
        "invalid options: formula expression nesting exceeds 128"
    );
    assert_matches!(
        error,
        DeseqError::FormulaLimitExceeded {
            resource: FormulaResource::ExpressionDepth,
            limit: 128,
            observed: Some(129),
        }
    );
    let mut value = Value::scalar(1.0);
    for _ in 0..129 {
        value = Value::new(Data::List(vec![value]));
    }
    assert_matches!(
        evaluate(&Expression::Literal(value), &[], &options).unwrap_err(),
        DeseqError::FormulaLimitExceeded {
            resource: FormulaResource::ValueDepth,
            limit: 128,
            observed: Some(129),
        }
    );
}

#[test]
fn sum_cancellation_retains_r_argument_rounding_boundaries() {
    if !(rsdeseq2::numeric_profile() == rsdeseq2::NumericProfile::RLinuxX86) {
        // These boundaries describe R's x87 accumulation only.
        return;
    }
    let options = EvaluationOptions::default();
    // Original R rsum accumulates one vector in LDOUBLE, then do_summary
    // combines the individually rounded arguments in binary64.
    let vector = Expression::Literal(Value::numbers(vec![1e16, 1.0, -1e16]));
    let single = evaluate(&call("sum", vec![vector]), &[], &options).unwrap();
    assert_eq!(single, Value::scalar(1.0));
    let arguments = [1e16, 1.0, -1e16]
        .into_iter()
        .map(|x| Expression::Literal(Value::scalar(x)))
        .collect();
    let separate = evaluate(&call("sum", arguments), &[], &options).unwrap();
    assert_eq!(separate, Value::scalar(0.0));
}
