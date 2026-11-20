//! Independent original-ashr search traces and analytical location estimates.
use crate::support::posterior as support;
use support::*;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/ashr_mode.rs"]
mod ashr_mode;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
use rsdeseq2::shrinkage::*;
use rsdeseq2::{NumericProfile, numeric_profile};

#[test]
fn estimated_modes_and_complete_summaries_match_original_profiles() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/ashr_mode_r462.json")).unwrap();
    check_reference(&reference["cases"]);
}
#[test]
fn mode_searches_preserve_original_recycled_and_expression_warnings() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/likelihood_warnings_r462.json")).unwrap();
    check_reference(&reference["mode_cases"]);
}
fn check_reference(cases: &serde_json::Value) {
    let mut failures = vec![];
    for c in cases.as_array().unwrap() {
        let input: ashr_mode::Case = serde_json::from_value(c.clone()).unwrap();
        let result = ashr_mode::run(&input);
        if c["expected_error"].is_string() || c["native_expected_error"].is_string() {
            assert!(result.is_err(), "{}", input.data.name);
            continue;
        }
        let (mode, result) = result.unwrap_or_else(|e| panic!("{}: {e}", input.data.name));
        // Small changes in the objective can take Brent's interpolation down
        // a different path. Compare the two affected portable summaries with
        // independent ashr fits at fixed coordinates near their optima.
        let portable = c
            .get("portable_profile")
            .filter(|_| numeric_profile() == NumericProfile::Portable);
        let (reference, expected_mode) =
            portable.map_or((c, &c["expected_mode"]), |p| (p, &p["mode"]));
        let expected_mode = parse_numeric_input(expected_mode.as_str().unwrap()).unwrap();
        assert!(
            (mode.mode - expected_mode).abs() < 2e-6,
            "{}: mode {} vs {}",
            input.data.name,
            mode.mode,
            expected_mode
        );
        let bounds: Vec<f64> = c["search_bounds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap())
            .collect();
        assert_eq!(
            (mode.bounds.0.to_bits(), mode.bounds.1.to_bits()),
            (bounds[0].to_bits(), bounds[1].to_bits()),
            "{}",
            input.data.name
        );
        let trials = mode.evaluations + 1;
        let trial_warnings = reference["trial_warnings"].as_array();
        if let Some(trial_warnings) = trial_warnings {
            assert_eq!(trials, 23, "{}", input.data.name);
            let expected: Vec<_> = (0..trials)
                .flat_map(|_| trial_warnings.iter().cloned())
                .collect();
            assert_eq!(
                serde_json::json!(mode.warnings),
                serde_json::json!(expected)
            );
        } else {
            assert_eq!(
                trials,
                c["trials"].as_array().unwrap().len(),
                "{}",
                input.data.name
            );
        }
        for (field, actual) in [
            ("posterior_mean", &result.posterior_mean),
            ("posterior_sd", &result.posterior_sd),
            ("svalues", &result.svalues),
        ] {
            compare(
                actual,
                &reference[field],
                2e-6,
                0.0,
                &format!("{} {field}", input.data.name),
                &mut failures,
            );
        }
        let ll = parse_numeric_input(reference["log_likelihood"].as_str().unwrap()).unwrap();
        assert!(
            (mode.log_likelihood - ll).abs() < 2e-6 && (result.log_likelihood - ll).abs() < 2e-6,
            "{}: likelihood",
            input.data.name
        );
        let warnings = reference
            .get("warnings")
            .unwrap_or(&c["warnings"])
            .as_array()
            .unwrap()
            .iter()
            .filter(|v| !v.as_str().unwrap().starts_with("subscript out of bounds"))
            .cloned()
            .collect::<Vec<_>>();
        let summary_warnings = if trial_warnings.is_some() {
            assert!(result.warnings.starts_with(&mode.warnings));
            &result.warnings[mode.warnings.len()..]
        } else {
            &result.warnings
        };
        assert_eq!(
            serde_json::json!(summary_warnings),
            serde_json::json!(warnings),
            "{}",
            input.data.name
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn point_prior_profile_matches_precision_weighted_normal_mle_and_validates_bounds() {
    let x = [-2.0, 1.0, 3.0];
    let se = [1.0, 2.0, 0.5];
    let options = LocationShrinkageOptions {
        grid: MixtureGridOptions {
            scales: Some(vec![0.0]),
            ..Default::default()
        },
        ..Default::default()
    };
    let fit = estimate_location_mode(&x, &se, &options, PriorModeSearch::ObservationRange).unwrap();
    assert!((fit.mode - (10.25 / 5.25)).abs() < 1e-8);
    let boundary =
        estimate_location_mode(&x, &se, &options, PriorModeSearch::Interval(-1.0, 1.0)).unwrap();
    assert!(boundary.mode < 1.0 && 1.0 - boundary.mode < 0.0003);
    for search in [
        PriorModeSearch::Interval(1.0, 1.0),
        PriorModeSearch::Interval(2.0, 1.0),
        PriorModeSearch::Interval(f64::NEG_INFINITY, 1.0),
        PriorModeSearch::ExpandedObservationRange(f64::NAN),
    ] {
        assert!(estimate_location_mode(&x, &se, &options, search).is_err());
    }
    assert!(
        estimate_location_mode(&[1.0; 3], &se, &options, PriorModeSearch::ObservationRange)
            .is_err()
    );
    assert!(
        estimate_location_mode(
            &[f64::NAN, 1.0, 2.0],
            &se,
            &options,
            PriorModeSearch::ObservationRange
        )
        .is_err()
    );
    assert!(
        estimate_location_mode(
            &x,
            &[1.0, -1.0, 0.5],
            &options,
            PriorModeSearch::ObservationRange
        )
        .is_err()
    );
    assert!(
        estimate_location_mode(
            &x,
            &se,
            &LocationShrinkageOptions {
                fixed_prior: true,
                ..options
            },
            PriorModeSearch::ObservationRange
        )
        .is_err()
    );
}

#[test]
fn likelihood_only_search_does_not_evaluate_invalid_posterior_moments() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/ashr_mode_r462.json")).unwrap();
    let c = reference["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["name"] == "custom_logistic")
        .unwrap();
    let input: ashr_mode::Case = serde_json::from_value(c.clone()).unwrap();
    let mut options = location_shrinkage::options(&input.data).unwrap();
    let x = parse_numeric_inputs(&input.data.estimates).unwrap();
    let se = parse_numeric_inputs(&input.data.standard_errors).unwrap();
    let good =
        estimate_location_mode(&x, &se, &options, PriorModeSearch::ObservationRange).unwrap();
    let mut model = c["custom_likelihood"].clone();
    model["first"] = serde_json::json!({"arguments":["a","b"],"constants":{},"expression":{
        "type":"call","function_name":"log","arguments":[{"name":null,"expression":{
            "type":"literal","value":{"type":"number","integer":false,"values":[bits(-1.0)],
                "dimensions":null,"names":null,"dimnames":[null,null],"attributes":[]}}}]}});
    options.likelihood = LocationLikelihood::Restricted(std::sync::Arc::new(
        likelihood::restricted(&model).unwrap(),
    ));
    let selected =
        estimate_location_mode(&x, &se, &options, PriorModeSearch::ObservationRange).unwrap();
    assert_eq!(selected.mode.to_bits(), good.mode.to_bits());
    assert_eq!(
        selected.log_likelihood.to_bits(),
        good.log_likelihood.to_bits()
    );
    assert_eq!(selected.warnings, good.warnings);
}
