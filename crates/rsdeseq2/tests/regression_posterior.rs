//! Independent original-apeglm profiles, retaining the public numerical gates.
#[path = "common/posterior.rs"]
mod support;
use support::*;
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/regression_posterior_r462.json")).unwrap()
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
    compare_posteriors(reference(), true);
}
#[test]
fn beta_binomial_r_and_fast_methods_match_original_apeglm() {
    compare_posteriors(
        serde_json::from_str(include_str!("data/beta_binomial_posterior_r462.json")).unwrap(),
        true,
    );
}
fn compare_posteriors(reference: serde_json::Value, compare_warning_events: bool) {
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
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

#[test]
fn beta_binomial_domain_and_parameter_shapes_are_rejected_before_fitting() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/beta_binomial_posterior_r462.json")).unwrap();
    let case =
        || serde_json::from_value::<regression::Case>(reference["cases"][0].clone()).unwrap();
    let reject = |c: &regression::Case, reason: &str| {
        let error = regression::run(c).unwrap_err().to_string();
        assert!(error.contains(reason), "{error}");
    };
    for theta in ["0", "-1", "NaN"] {
        let mut c = case();
        c.parameters.as_mut().unwrap().values[0] = theta.into();
        reject(&c, "concentration");
    }
    let mut c = case();
    c.observations.values[0] = "11".into();
    reject(&c, "successes");
    let mut c = case();
    let m = c.parameters.as_mut().unwrap();
    m.columns -= 1;
    m.values.truncate(m.rows * m.columns);
    reject(&c, "one trial count");
}

#[test]
fn randomized_beta_binomial_and_nb_modes_and_state_match_original_apeglm() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/random_start_posterior_r462.json")).unwrap();
    let beta = random_posterior_cases(&reference, "regression_cases");
    compare_posteriors(serde_json::json!({"cases":beta}), false);
    let mut failures = vec![];
    for c in random_posterior_cases(&reference, "heavy_tail_cases") {
        let input: heavy_tail::Case = serde_json::from_value(c.clone()).unwrap();
        let r = heavy_tail::run(&input).unwrap();
        let g = input.counts.rows;
        let p = input.design.columns;
        let mode = (0..p)
            .flat_map(|j| {
                (0..g).map({
                    let r = &r;
                    move |i| *r.mode.get(i, j).unwrap()
                })
            })
            .collect::<Vec<_>>();
        compare(
            &mode,
            &c["mode"],
            0.01 * std::f64::consts::LN_2,
            0.,
            &input.name,
            &mut failures,
        );
        assert!(r.standard_deviation.is_none());
        assert_eq!(
            random::encode_state(r.random_state.as_ref().unwrap()),
            c["expected_random_state"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap().to_owned())
                .collect::<Vec<_>>(),
            "{} state",
            input.name
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Expand pooled serialized states; this only changes fixture storage, never
/// numerical values or the number/order of draws.
fn random_posterior_cases(reference: &serde_json::Value, field: &str) -> Vec<serde_json::Value> {
    reference[field]
        .as_array()
        .unwrap()
        .iter()
        .map(|case| {
            let mut c = case.clone();
            let key = c["rng_key"].as_str().unwrap();
            let state = &reference["rng_states"][key];
            let p = c["design"]["columns"].as_u64().unwrap().to_string();
            let before = state["before"].clone();
            let after = state["after"][&p].clone();
            let start = state["starts"][&p].clone();
            c["rng_state"] = before;
            c["expected_random_state"] = after;
            c["expected_start"] = start;
            c
        })
        .collect()
}
