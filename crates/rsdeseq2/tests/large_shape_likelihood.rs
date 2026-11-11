//! Original-package tails/moments plus independent high-precision density truths.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "common/builtin_checks.rs"]
mod checks;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/large_shape_likelihood_r462.json")).unwrap()
}
#[test]
fn large_shape_tails_and_moments_match_original_without_iteration_exhaustion() {
    checks::kernels(reference());
}
#[test]
fn large_shape_density_normalization_matches_independent_hundred_digit_truths() {
    let reference = reference();
    let truths: serde_json::Value =
        serde_json::from_str(include_str!("data/large_shape_density_truth.json")).unwrap();
    let mut failures = vec![];
    let inputs = reference["kernel_cases"].as_array().unwrap();
    let rows = truths["cases"].as_array().unwrap();
    assert_eq!(inputs.len(), rows.len());
    for (c, truth) in inputs.iter().zip(rows) {
        assert_eq!(c["name"], truth["name"]);
        let input: builtin_likelihood::Case = serde_json::from_value(c.clone()).unwrap();
        let r = builtin_likelihood::run(&input).unwrap();
        compare(
            &r[1],
            &truth["log_density"],
            2e-8,
            2e-12,
            input.name.as_str(),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn symmetric_large_beta_cdf_has_an_exact_half_mass() {
    use rsdeseq2::shrinkage::{BinomialLink, BuiltinLocationLikelihood};
    for shape in [1e5, 1e8, 1e10, 1e14] {
        let likelihood = BuiltinLocationLikelihood::binomial(
            vec![shape - 1.0],
            vec![2.0 * (shape - 1.0)],
            BinomialLink::Identity,
        )
        .unwrap();
        let expected = -std::f64::consts::LN_2 - (2.0 * shape - 1.0).ln();
        assert!((likelihood.log_cdf(0.5) - expected).abs() < 2e-12);
    }
}
