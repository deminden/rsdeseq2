//! Independent canonical capture layout, component order and missing-row policies.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[path = "common/builtin_checks.rs"]
mod checks;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/normalmix_capture_r462.json")).unwrap()
}
#[test]
fn capture_kernels_match_original_matrix_coercion_and_recycling() {
    checks::kernels(reference());
}
#[test]
fn capture_fits_match_original_fitted_fixed_and_missing_row_summaries() {
    checks::complete_fits(reference());
}
#[test]
fn equal_capture_lengths_flatten_columns_without_row_weight_interpretation() {
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::BuiltinLocationLikelihood,
    };
    let weights = [0.05, 0.2, 0.1, 0.15, 0.25, 0.25];
    let scales = [0.3, 1.2, 0.5, 1.6, 0.7, 2.0];
    let w = MatrixView::new(3, 2, &weights, MatrixOrder::ColumnMajor).unwrap();
    let s = MatrixView::new(2, 3, &scales, MatrixOrder::ColumnMajor).unwrap();
    let likelihood = BuiltinLocationLikelihood::normal_mixture_from_captures(w, s, true).unwrap();
    let phi = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
    let expected: f64 = weights.iter().zip(scales).map(|(w, s)| w * phi / s).sum();
    assert!((likelihood.log_density(0.0).exp() - expected).abs() < 2e-15);
    assert_eq!(
        likelihood
            .for_observation(0, 2)
            .unwrap()
            .log_density(0.0)
            .to_bits(),
        likelihood
            .for_observation(1, 2)
            .unwrap()
            .log_density(0.0)
            .to_bits()
    );
}
