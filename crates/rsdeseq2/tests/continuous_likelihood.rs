//! Original built-in likelihood kernels, complete workflows and analytical checks.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::BuiltinLocationLikelihood,
};
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/continuous_likelihood_r462.json")).unwrap()
}
#[test]
fn continuous_kernels_match_original_scalar_and_row_parameters() {
    let mut failures = vec![];
    for c in reference()["kernel_cases"].as_array().unwrap() {
        let input: builtin_likelihood::Case = serde_json::from_value(c.clone()).unwrap();
        let result = builtin_likelihood::run(&input);
        if c["expected_error"].is_string() {
            assert!(result.is_err(), "{}", input.name);
            continue;
        }
        let r = result.unwrap();
        for (i, name) in ["log_cdf", "log_density", "first", "second"]
            .iter()
            .enumerate()
        {
            compare(
                &r[i],
                &c[name],
                if i < 2 { 2e-8 } else { 2e-6 },
                2e-12,
                &format!("{} {name}", input.name),
                &mut failures,
            );
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn complete_continuous_fits_preserve_original_summaries_and_missingness() {
    let mut failures = vec![];
    for c in reference()["cases"].as_array().unwrap() {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let r = location_shrinkage::run(&input);
        if c["expected_error"].is_string() {
            assert!(r.is_err(), "{}", input.name);
            continue;
        }
        let r = r.unwrap_or_else(|e| panic!("{}: {e}", input.name));
        for (name, v) in [
            ("posterior_mean", &r.posterior_mean),
            ("posterior_sd", &r.posterior_sd),
            ("svalues", &r.svalues),
            ("cdf", &r.cdf),
        ] {
            compare(
                v,
                &c[name],
                2e-6,
                0.0,
                &format!("{} {name}", input.name),
                &mut failures,
            );
        }
        compare(
            &[r.log_likelihood],
            &serde_json::json!([c["log_likelihood"]]),
            2e-6,
            0.0,
            &format!("{} likelihood", input.name),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn logistic_identity_and_observation_geometry_have_independent_solutions() {
    let l = BuiltinLocationLikelihood::log_f(vec![2.0], vec![2.0]).unwrap();
    assert!((l.log_cdf(0.0) + 2.0_f64.ln()).abs() < 1e-13);
    assert!((l.log_density(0.0) + 4.0_f64.ln()).abs() < 1e-13);
    assert!(l.truncated_moment(-2.0, 2.0, 1).abs() < 1e-12);
    assert_eq!(l.truncated_moment(-2.0, 2.0, 2).to_bits() as u32, 1954);
    let scales = [1.0, 3.0, 2.0, 4.0];
    let n = BuiltinLocationLikelihood::normal_mixture_by_observation(
        vec![0.5, 0.5],
        MatrixView::new(2, 2, &scales, MatrixOrder::ColumnMajor).unwrap(),
    )
    .unwrap();
    let phi = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
    assert!((n.for_observation(0, 2).unwrap().log_density(0.0).exp() - phi * 0.75).abs() < 1e-15);
    assert!(
        (n.for_observation(1, 2).unwrap().log_density(0.0).exp() - phi * 7.0 / 24.0).abs() < 1e-15
    );
    assert!(n.for_observation(0, 3).is_err());
    assert!(n.for_observation(2, 2).is_err());
    assert!(BuiltinLocationLikelihood::normal_mixture(vec![1.0], vec![-1.0]).is_err());
    assert!(BuiltinLocationLikelihood::log_f(vec![0.0], vec![2.0]).is_err());
    assert!(n.truncated_moment(f64::NAN, 1.0, 1).is_nan());
}
