use crate::support as common;
use common::hex_bits as bits;
use rsdeseq2::{NumericProfile, numeric_profile};
use std::error::Error;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/adaptive.rs"]
mod adaptive;
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
        let portable_uniform =
            numeric_profile() == NumericProfile::Portable && input.name == "undefined_uniform";
        for (field, actual) in [
            ("posterior_mean", Some(r.posterior_mean)),
            ("posterior_sd", Some(r.posterior_sd)),
            ("svalues", r.svalues),
        ] {
            let (expected, tolerance) = if portable_uniform {
                // At the two extremes the widest component dominates, and the
                // distance from its nearest endpoint is effectively exponential.
                let tail_scale = 0.001_f64.powi(2) / (1e6 - 4.0);
                // At zero both uniforms cover the entire observable Gaussian.
                let continuous = (0.3 / 2.0 + 0.5 / 8.0) * 0.001 * std::f64::consts::TAU.sqrt();
                let (values, tolerance) = match field {
                    "posterior_mean" => ([-4.0 + tail_scale, 0.0, 4.0 - tail_scale], 1e-15),
                    "posterior_sd" => (
                        [
                            tail_scale,
                            0.001 * (continuous / (0.2 + continuous)).sqrt(),
                            tail_scale,
                        ],
                        1e-18,
                    ),
                    // Only the central estimate has a local false-sign-or-small rate of 1.
                    "svalues" => ([0.0, 1.0 / 3.0, 0.0], f64::EPSILON),
                    _ => unreachable!(),
                };
                (Some(values.to_vec()), tolerance)
            } else {
                (
                    c[field].as_array().map(|values| {
                        values
                            .iter()
                            .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap())
                            .collect()
                    }),
                    2e-6,
                )
            };
            if actual.is_none() != expected.is_none() {
                failures.push(format!("{} {field}: missing field mismatch", input.name));
                continue;
            }
            if let (Some(actual), Some(expected)) = (actual, expected) {
                assert_eq!(actual.len(), expected.len(), "{} {field}", input.name);
                for (i, (a, b)) in actual.into_iter().zip(expected).enumerate() {
                    let same_missing = a.is_nan()
                        && b.is_nan()
                        && (a.to_bits() as u32 == 1954) == (b.to_bits() as u32 == 1954);
                    if a != b
                        && !same_missing
                        && !(a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance)
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
