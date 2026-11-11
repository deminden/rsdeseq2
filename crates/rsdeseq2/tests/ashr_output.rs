//! Original-package selection/missingness references and independent contracts.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/ashr_output.rs"]
mod ashr_output;
use rsdeseq2::shrinkage::*;

#[test]
fn selected_fields_and_missing_rates_match_original_ashr() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/ashr_output_r462.json")).unwrap();
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
        let input: ashr_output::Case = serde_json::from_value(c.clone()).unwrap();
        let result = ashr_output::run(&input);
        if c["expected_error"].is_string() {
            assert!(result.is_err(), "{}", input.name);
            continue;
        }
        let r = result.unwrap_or_else(|e| panic!("{}: {e}", input.name));
        for (field, actual) in [
            ("posterior_mean", r.posterior_mean.as_deref()),
            ("posterior_sd", r.posterior_sd.as_deref()),
            ("svalues", r.svalues.as_deref()),
        ] {
            assert_eq!(
                actual.is_none(),
                c[field].is_null(),
                "{} {field}",
                input.name
            );
            if let Some(actual) = actual {
                compare(
                    actual,
                    &c[field],
                    2e-6,
                    0.0,
                    &format!("{} {field}", input.name),
                    &mut failures,
                );
            }
        }
        assert_eq!(
            r.fitted_prior,
            c["fitted_prior"].as_bool().unwrap(),
            "{}",
            input.name
        );
        assert_eq!(
            serde_json::json!(r.fit.warnings),
            c["warnings"],
            "{}",
            input.name
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn no_selected_summary_does_not_assign_into_a_different_population() {
    let x = [-1.0, 0.5, 2.0];
    let query = AdaptiveShrinkageQuery {
        estimates: &[0.0],
        standard_errors: &[1.0],
        lfc_threshold: 0.0,
        svalues: false,
    };
    let options = LocationShrinkageOptions::default();
    let fit = adaptive_shrinkage_with_output(
        &x,
        &[1.0],
        &options,
        query,
        AdaptiveOutputSelection::level(0.0).unwrap(),
    )
    .unwrap();
    assert!(fit.posterior_mean.is_none());
    assert!(fit.fitted_prior);
    assert!(
        adaptive_shrinkage_with_output(
            &x,
            &[1.0],
            &options,
            query,
            AdaptiveOutputSelection::names(["PosteriorSD"])
        )
        .is_err()
    );
    assert!(
        adaptive_shrinkage_with_output(
            &x,
            &[1.0],
            &options,
            AdaptiveShrinkageQuery {
                lfc_threshold: f64::NAN,
                ..query
            },
            AdaptiveOutputSelection::default()
        )
        .is_err()
    );
}
