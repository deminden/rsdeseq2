//! Shared original-package posterior, missingness, state and diagnostic gates.
use crate::support::*;
pub fn compare_posteriors(reference: serde_json::Value, compare_warning_events: bool) {
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
        let input: regression::Case = serde_json::from_value(c.clone()).unwrap();
        let actual = regression::run(&input);
        if c["expected_error"].is_string() {
            let expected =
                c["native_expected_error"]
                    .as_str()
                    .unwrap_or_else(|| match input.name.as_str() {
                        "invalid_shape" => "nonconformable likelihood matrix product",
                        "null_param" => "parameters must be supplied",
                        _ => panic!("missing rejection category for {}", input.name),
                    });
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
        if let Some(state) = &r.random_state {
            assert_eq!(
                random::encode_state(state),
                c["expected_random_state"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_str().unwrap().to_owned())
                    .collect::<Vec<_>>(),
                "{} random state",
                input.name
            );
        }
        let warning_events = |v: &[String]| {
            v.iter()
                .map(|v| v.split(',').next().unwrap().to_owned())
                .collect::<Vec<_>>()
        };
        let expected_warnings = c["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_owned())
            .collect::<Vec<_>>();
        // Source counts vary near floating-point plateaus for random starts.
        // Those observations are retained in the workflow report; native and
        // wrapped warning data must still match exactly. Keep the pre-existing
        // deterministic-profile event gates intact.
        if compare_warning_events {
            assert_eq!(
                warning_events(&r.warnings),
                warning_events(&expected_warnings),
                "{} warning events",
                input.name
            );
        }
        let reorder = |m: &rsdeseq2::matrix::RowMajorMatrix<f64>| {
            (0..input.design.columns)
                .flat_map(|j| (0..input.observations.rows).map(move |i| *m.get(i, j).unwrap()))
                .collect::<Vec<_>>()
        };
        for (field, values, bound) in [
            ("mode", reorder(&r.mode), 0.01 * std::f64::consts::LN_2),
            (
                "sd",
                r.standard_deviation
                    .as_ref()
                    .map(reorder)
                    .unwrap_or_else(|| {
                        vec![
                            f64::from_bits(0x7ff0_0000_0000_07a2);
                            input.observations.rows * input.design.columns
                        ]
                    }),
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
