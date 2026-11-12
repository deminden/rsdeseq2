//! Original built-in likelihood kernels, complete workflows and analytical checks.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
use rsdeseq2::shrinkage::BuiltinLocationLikelihood;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/ashr_mode.rs"]
mod ashr_mode;
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/count_likelihood_r462.json")).unwrap()
}
#[path = "common/builtin_checks.rs"]
mod checks;
#[test]
fn count_kernels_match_original_scalar_and_row_parameters() {
    checks::kernels(reference());
}
#[test]
fn complete_count_fits_preserve_original_summaries_and_missingness() {
    checks::complete_fits(reference());
}
#[test]
fn count_kernels_have_analytic_exponential_and_uniform_beta_solutions() {
    use rsdeseq2::shrinkage::{BinomialLink, PoissonLink};
    let l =
        BuiltinLocationLikelihood::poisson(vec![0.0], vec![2.0], PoissonLink::Identity).unwrap();
    let e = (-2.0_f64).exp();
    assert!((l.log_cdf(1.0).exp() - (1.0 - e) / 2.0).abs() < 1e-14);
    assert!((l.log_density(-1.0) + 2.0).abs() < 1e-14);
    assert!((l.truncated_moment(-1.0, 0.0, 1) + (1.0 - 3.0 * e) / (2.0 * (1.0 - e))).abs() < 1e-14);
    assert!(
        (l.truncated_moment(-1.0, 0.0, 2) - (2.0 - 10.0 * e) / (4.0 * (1.0 - e))).abs() < 1e-14
    );
    let l =
        BuiltinLocationLikelihood::binomial(vec![0.0], vec![0.0], BinomialLink::Identity).unwrap();
    assert!((l.log_cdf(-0.3).exp() - 0.3).abs() < 1e-14);
    assert!(l.log_density(0.5).abs() < 1e-14);
    assert!((l.truncated_moment(-1.0, 0.0, 1) + 0.5).abs() < 1e-14);
    assert!((l.truncated_moment(-1.0, 0.0, 2) - 1.0 / 3.0).abs() < 1e-14);
    assert!(
        BuiltinLocationLikelihood::binomial(vec![2.0], vec![1.0], BinomialLink::Identity).is_err()
    );
    assert!(
        BuiltinLocationLikelihood::poisson(vec![0.0], vec![0.0], PoissonLink::Identity).is_err()
    );
}

#[test]
fn count_mode_bounds_searches_and_complete_summaries_match_original() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/count_mode_r462.json")).unwrap();
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
        let input: ashr_mode::Case = serde_json::from_value(c.clone()).unwrap();
        let result = ashr_mode::run(&input);
        if c["expected_error"].is_string() {
            assert!(result.is_err(), "{}", input.data.name);
            continue;
        }
        let (mode, r) = result.unwrap_or_else(|e| panic!("{}: {e}", input.data.name));
        compare(
            &[mode.mode],
            &serde_json::json!([c["expected_mode"]]),
            2e-6,
            0.0,
            &format!("{} mode", input.data.name),
            &mut failures,
        );
        let bounds = c["search_bounds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap().to_bits())
            .collect::<Vec<_>>();
        assert_eq!(
            [mode.bounds.0.to_bits(), mode.bounds.1.to_bits()],
            bounds.as_slice(),
            "{}",
            input.data.name
        );
        assert_eq!(
            mode.evaluations + 1,
            c["trials"].as_array().unwrap().len(),
            "{}",
            input.data.name
        );
        for (field, v) in [
            ("posterior_mean", &r.posterior_mean),
            ("posterior_sd", &r.posterior_sd),
            ("svalues", &r.svalues),
        ] {
            compare(
                v,
                &c[field],
                2e-6,
                0.0,
                &format!("{} {field}", input.data.name),
                &mut failures,
            );
        }
        compare(
            &[r.log_likelihood],
            &serde_json::json!([c["log_likelihood"]]),
            2e-6,
            0.0,
            &format!("{} log likelihood", input.data.name),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn count_batch_orientation_is_shared_by_the_direct_posterior_api() {
    let mut failures = vec![];
    for c in reference()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| c["name"].as_str().unwrap().contains("global_orientation"))
    {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let options = location_shrinkage::options(&input).unwrap();
        let grid = options.supplied_grid.as_ref().unwrap();
        let prior = grid
            .prior_with_likelihood(&grid.initial_weights, &options.likelihood)
            .unwrap();
        let posteriors = prior
            .posteriors(
                &parse_numeric_inputs(&input.estimates).unwrap(),
                &parse_numeric_inputs(&input.standard_errors).unwrap(),
            )
            .unwrap();
        let means = posteriors.iter().map(|p| p.mean()).collect::<Vec<_>>();
        compare(
            &means,
            &c["posterior_mean"],
            2e-6,
            0.0,
            input.name.as_str(),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
