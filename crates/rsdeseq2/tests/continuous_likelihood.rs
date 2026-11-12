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
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::BuiltinLocationLikelihood,
};
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/continuous_likelihood_r462.json")).unwrap()
}
#[path = "common/builtin_checks.rs"]
mod checks;
#[test]
fn continuous_kernels_match_original_scalar_and_row_parameters() {
    checks::kernels(reference());
}
#[test]
fn complete_continuous_fits_preserve_original_summaries_and_missingness() {
    checks::complete_fits(reference());
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
