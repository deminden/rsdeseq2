//! Independent original-apeglm profiles, retaining the public numerical gates.
use std::error::Error;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/formula.rs"]
mod formula;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/heavy_tail.rs"]
mod heavy_tail;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/likelihood.rs"]
mod likelihood;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/regression.rs"]
mod regression;
fn bits(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}
fn parse_numeric_input(s: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(s) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(s, 16)?)
    } else if s == "NA" {
        f64::from_bits(0x7ff0_0000_0000_07a2)
    } else {
        s.parse()?
    })
}
fn parse_numeric_inputs(v: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    v.iter().map(|s| parse_numeric_input(s)).collect()
}
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/regression_posterior_r462.json")).unwrap()
}
fn compare(
    a: &[f64],
    b: &serde_json::Value,
    absolute: f64,
    relative: f64,
    context: &str,
    failures: &mut Vec<String>,
) {
    let b = b.as_array().unwrap();
    assert_eq!(a.len(), b.len(), "{context}");
    for (i, (&a, b)) in a.iter().zip(b).enumerate() {
        let b = parse_numeric_input(b.as_str().unwrap()).unwrap();
        if a == b
            || (a.is_nan()
                && b.is_nan()
                && (a.to_bits() as u32 == 1954) == (b.to_bits() as u32 == 1954))
        {
            continue;
        }
        if a.is_finite() && b.is_finite() && (a - b).abs() <= absolute + relative * b.abs() {
            continue;
        }
        failures.push(format!(
            "{context}[{i}]: {a:.17e} != {b:.17e}; error {:.3e}",
            (a - b).abs()
        ));
    }
}
#[test]
fn regression_density_expressions_preserve_independent_precision() {
    let mut failures = vec![];
    let truth: serde_json::Value =
        serde_json::from_str(include_str!("data/regression_density_truth.json")).unwrap();
    let mut reference_differences = vec![];
    for c in reference()["leaves"].as_array().unwrap() {
        let input: regression::Leaf = serde_json::from_value(c.clone()).unwrap();
        let actual = regression::leaf(&input).unwrap();
        let high_precision = truth["leaves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|v| v["name"] == c["name"]);
        let expected = high_precision
            .map(|v| &v["expected"])
            .unwrap_or(&c["expected"]);
        compare(&actual, expected, 2e-10, 2e-12, &input.name, &mut failures);
        if let Some(t) = high_precision {
            let accurate = t["expected"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap())
                .collect::<Vec<_>>();
            let mut differences = vec![];
            compare(
                &accurate,
                &c["expected"],
                2e-10,
                2e-12,
                &input.name,
                &mut differences,
            );
            reference_differences.extend(differences);
        }
    }
    // The original R values remain frozen. Nine observations exceed the new
    // leaf gate because R loses tail precision or cancels large NB terms.
    // An independent 90-decimal integral/log-gamma reference certifies these
    // leaves at the SAME gate; public posterior gates remain unchanged.
    assert_eq!(reference_differences.len(), 9, "{reference_differences:?}");
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
#[test]
fn general_and_analytic_nb_regression_match_original_apeglm() {
    let mut failures = vec![];
    for c in reference()["cases"].as_array().unwrap() {
        let input: regression::Case = serde_json::from_value(c.clone()).unwrap();
        let actual = regression::run(&input);
        if c["expected_error"].is_string() {
            let expected = match input.name.as_str() {
                "invalid_shape" => "nonconformable likelihood matrix product",
                "null_param" => "parameters must be supplied",
                _ => panic!("missing rejection category for {}", input.name),
            };
            match actual {
                Ok(_) => failures.push(format!("{} accepted invalid source request", input.name)),
                Err(e) if !e.to_string().contains(expected) => {
                    failures.push(format!("{} rejected for wrong reason: {e}", input.name))
                }
                Err(_) => {}
            }
            continue;
        }
        let r = match actual {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{} unexpected error {e}", input.name));
                continue;
            }
        };
        let reorder = |m: &rsdeseq2::matrix::RowMajorMatrix<f64>| {
            (0..input.design.columns)
                .flat_map(|j| (0..input.observations.rows).map(move |i| *m.get(i, j).unwrap()))
                .collect::<Vec<_>>()
        };
        for (field, values, bound) in [
            ("mode", reorder(&r.mode), 0.01 * std::f64::consts::LN_2),
            (
                "sd",
                reorder(r.standard_deviation.as_ref().unwrap()),
                0.01 * std::f64::consts::LN_2,
            ),
            (
                "fsr",
                r.false_sign_rate
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
                    .collect(),
                0.02,
            ),
            (
                "sign_svalues",
                r.sign_svalues
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
                    .collect(),
                0.02,
            ),
            (
                "fsos",
                r.false_sign_or_small_rate
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
                    .collect(),
                0.05,
            ),
            (
                "threshold_svalues",
                r.threshold_svalues
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
                    .collect(),
                0.05,
            ),
        ] {
            compare(
                &values,
                &c[field],
                bound,
                0.,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        compare(
            &[r.prior.scale],
            &serde_json::json!([c["fitted_scale"]]),
            1e-12,
            0.,
            &format!("{} prior scale", input.name),
            &mut failures,
        );
        if let Some(v) = r.prior.adapted_variance {
            compare(
                &[v],
                &serde_json::json!([c["prior_variance"]]),
                1e-12,
                0.,
                &format!("{} prior variance", input.name),
                &mut failures,
            );
        } else {
            assert!(c["prior_variance"].is_null());
        }
    }
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
