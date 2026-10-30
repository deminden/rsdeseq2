use rsdeseq2::design::expression::{
    Argument, Data, EvaluationOptions, Expression, Value, evaluate,
};
use serde_json::Value as Json;

fn strings(v: &Json) -> Option<Vec<String>> {
    v.as_array()
        .map(|a| a.iter().map(|v| v.as_str().unwrap().into()).collect())
}
fn value(v: &Json) -> Value {
    let list = v["values"].as_array();
    let data = match v["type"].as_str().unwrap() {
        "number" => Data::Number {
            integer: v["integer"].as_bool().unwrap(),
            values: list
                .unwrap()
                .iter()
                .map(|v| f64::from_bits(u64::from_str_radix(v.as_str().unwrap(), 16).unwrap()))
                .collect(),
        },
        "text" => Data::Text(
            list.unwrap()
                .iter()
                .map(|v| v.as_str().map(str::to_owned))
                .collect(),
        ),
        "logical" => Data::Logical(list.unwrap().iter().map(Json::as_bool).collect()),
        "list" => Data::List(list.unwrap().iter().map(value).collect()),
        "null" => Data::Null,
        _ => panic!(),
    };
    let mut result = Value::new(data);
    result.names = strings(&v["names"]);
    result.dimnames = [strings(&v["dimnames"][0]), strings(&v["dimnames"][1])];
    result.dimensions = v["dimensions"].as_array().map(|d| {
        [
            d[0].as_u64().unwrap() as usize,
            d[1].as_u64().unwrap() as usize,
        ]
    });
    if let Some(a) = v["attributes"].as_object() {
        result.attributes = a.iter().map(|(n, v)| (n.clone(), value(v))).collect()
    }
    result
}
fn expression(v: &Json) -> Expression {
    match v["type"].as_str().unwrap() {
        "symbol" => Expression::Symbol(v["name"].as_str().unwrap().into()),
        "literal" => Expression::Literal(value(&v["value"])),
        "call" => Expression::Call {
            function: v["function_name"].as_str().unwrap().into(),
            arguments: v["arguments"]
                .as_array()
                .unwrap()
                .iter()
                .map(|a| Argument {
                    name: a["name"].as_str().map(str::to_owned),
                    expression: expression(&a["expression"]),
                })
                .collect(),
        },
        _ => panic!(),
    }
}
fn equal_value(actual: &Value, expected: &Value, path: &str) {
    assert_eq!(actual.dimensions, expected.dimensions, "{path}: dimensions");
    assert_eq!(actual.names, expected.names, "{path}: names");
    assert_eq!(actual.dimnames, expected.dimnames, "{path}: dimnames");
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
    let reference: Json =
        serde_json::from_str(include_str!("data/formula_expression_r462.json")).unwrap();
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
