#[macro_use]
#[path = "support/formula_reference.rs"]
mod reference;
use rsdeseq2::design::expression::{
    Argument, Data, EvaluationOptions, Expression, Value, evaluate,
};
use serde_json::Value as Json;

#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/formula.rs"]
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
}

#[test]
fn sum_preserves_r_missing_precedence_across_nan_orders() {
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    let linux_x86 = cfg!(all(target_os = "linux", target_arch = "x86_64"));
    let generated_nan = if linux_x86 {
        0xfff8_0000_0000_0000
    } else {
        f64::NAN.to_bits()
    };
    let options = EvaluationOptions::default();
    for (values, expected) in [
        (vec![na, f64::NAN], 0x7ff8_0000_0000_07a2),
        (
            vec![f64::NAN, na],
            if linux_x86 {
                f64::NAN.to_bits()
            } else {
                0x7ff8_0000_0000_07a2
            },
        ),
        (
            vec![f64::INFINITY, f64::NEG_INFINITY, na],
            if linux_x86 {
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
fn mean_preserves_r_missing_and_nan_order() {
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    let options = EvaluationOptions::default();
    let linux_x86 = cfg!(all(target_os = "linux", target_arch = "x86_64"));
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
    let mut oversized = bytes;
    oversized[5..13].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(wire::decode(&oversized).is_err());
    let huge = call(
        "seq_len",
        vec![Expression::Literal(Value::scalar(i32::MAX as f64))],
    );
    assert!(
        evaluate(&huge, &[], &EvaluationOptions::default())
            .unwrap_err()
            .to_string()
            .contains("size limit")
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
    assert!(
        evaluate(&expression, &[], &options)
            .unwrap_err()
            .to_string()
            .contains("custom or unsupported formula function 'custom_function'")
    );
    for namespace in ["base::poly", "stats::log", "other::log"] {
        assert!(evaluate(&call(namespace, vec![]), &[], &options).is_err());
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
    assert!(
        evaluate(&nested, &[], &options)
            .unwrap_err()
            .to_string()
            .contains("nesting exceeds 128")
    );
}

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[test]
fn sum_cancellation_retains_r_argument_rounding_boundaries() {
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
