use std::error::Error;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/adaptive.rs"]
mod adaptive;
fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}
fn parse_numeric_input(s: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(s) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(s, 16)?)
    } else {
        match s {
            "NA" => f64::from_bits(0x7ff0_0000_0000_07a2),
            "NaN" => f64::NAN,
            _ => s.parse()?,
        }
    })
}
fn parse_numeric_inputs(v: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    v.iter().map(|s| parse_numeric_input(s)).collect()
}
#[test]
fn original_result_queries_and_recycling_match_original_wrapper() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/adaptive_result_r462.json")).unwrap();
    let cases = reference["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 10);
    let mut failures = vec![];
    for c in cases {
        let input: adaptive::Case = serde_json::from_value(c.clone()).unwrap();
        let r = adaptive::run(&input).unwrap();
        for (field, actual) in [
            ("posterior_mean", Some(r.posterior_mean)),
            ("posterior_sd", Some(r.posterior_sd)),
            ("svalues", r.svalues),
        ] {
            let expected = c[field].as_array();
            if actual.is_none() != expected.is_none() {
                failures.push(format!("{} {field}: missing field mismatch", input.name));
                continue;
            }
            if let (Some(actual), Some(expected)) = (actual, expected) {
                assert_eq!(actual.len(), expected.len(), "{} {field}", input.name);
                for (i, (a, b)) in actual.into_iter().zip(expected).enumerate() {
                    let b = parse_numeric_input(b.as_str().unwrap()).unwrap();
                    let same_missing = a.is_nan()
                        && b.is_nan()
                        && (a.to_bits() as u32 == 1954) == (b.to_bits() as u32 == 1954);
                    if a != b
                        && !same_missing
                        && !(a.is_finite() && b.is_finite() && (a - b).abs() <= 2e-6)
                    {
                        failures.push(format!("{} {field}[{i}]: {a:?} != {b:?}", input.name));
                    }
                }
            }
        }
        let warnings: Vec<String> = serde_json::from_value(c["warnings"].clone()).unwrap();
        if r.fit.warnings != warnings {
            failures.push(format!(
                "{} warnings: {:?} != {warnings:?}",
                input.name, r.fit.warnings
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
