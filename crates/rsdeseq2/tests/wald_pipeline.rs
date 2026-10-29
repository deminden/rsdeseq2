//! Numerical coverage grouped by responsibility; shared fixtures stay in this test crate.

use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn assert_wald_likelihood_state(fit: &DeseqFit, counts: &CountMatrix) {
    let log_like = fit.log_like.as_ref().unwrap();
    let full_deviance = fit.full_deviance.as_ref().unwrap();
    assert_eq!(log_like.len(), counts.n_genes());
    assert_eq!(full_deviance.len(), counts.n_genes());
    for (gene, (log_like, deviance)) in log_like.iter().zip(full_deviance).enumerate() {
        if log_like.is_nan() {
            assert!(
                deviance.is_nan(),
                "full deviance for gene {gene} should be NaN when log_like is NaN"
            );
        } else {
            assert_relative_eq!(*deviance, -2.0 * *log_like, epsilon = 1e-12);
        }
    }
}

fn assert_float_close_or_nan(actual: f64, expected: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
    } else {
        assert_relative_eq!(actual, expected, epsilon = 1e-12);
    }
}

fn assert_slice_close_or_nan(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_float_close_or_nan(*actual, *expected, &format!("{label}[{index}]"));
    }
}

fn assert_matrix_close_or_nan(
    actual: &RowMajorMatrix<f64>,
    expected: &RowMajorMatrix<f64>,
    label: &str,
) {
    assert_eq!(actual.n_rows(), expected.n_rows(), "{label}: row mismatch");
    assert_eq!(
        actual.n_cols(),
        expected.n_cols(),
        "{label}: column mismatch"
    );
    assert_slice_close_or_nan(actual.as_slice(), expected.as_slice(), label);
}

fn assert_wald_fit_intermediates_match(actual: &DeseqFit, expected: &DeseqFit, label: &str) {
    assert_eq!(actual.beta_converged, expected.beta_converged);
    assert_eq!(actual.beta_iter, expected.beta_iter);
    assert_matrix_close_or_nan(
        actual.beta.as_ref().unwrap(),
        expected.beta.as_ref().unwrap(),
        &format!("{label} beta"),
    );
    assert_matrix_close_or_nan(
        actual.beta_se.as_ref().unwrap(),
        expected.beta_se.as_ref().unwrap(),
        &format!("{label} beta_se"),
    );
    assert_matrix_close_or_nan(
        actual.beta_covariance.as_ref().unwrap(),
        expected.beta_covariance.as_ref().unwrap(),
        &format!("{label} beta_covariance"),
    );
    assert_matrix_close_or_nan(
        actual.mu.as_ref().unwrap(),
        expected.mu.as_ref().unwrap(),
        &format!("{label} mu"),
    );
    assert_matrix_close_or_nan(
        actual.hat_diagonal.as_ref().unwrap(),
        expected.hat_diagonal.as_ref().unwrap(),
        &format!("{label} hat_diagonal"),
    );
    assert_slice_close_or_nan(
        actual.log_like.as_ref().unwrap(),
        expected.log_like.as_ref().unwrap(),
        &format!("{label} log_like"),
    );
    assert_slice_close_or_nan(
        actual.full_deviance.as_ref().unwrap(),
        expected.full_deviance.as_ref().unwrap(),
        &format!("{label} full_deviance"),
    );
}

#[path = "wald_pipeline/aliases.rs"]
mod aliases;

#[path = "wald_pipeline/contrasts.rs"]
mod contrasts;

#[path = "wald_pipeline/filtering.rs"]
mod filtering;

#[path = "wald_pipeline/fitting.rs"]
mod fitting;

#[path = "wald_pipeline/formula_metadata.rs"]
mod formula_metadata;

#[path = "wald_pipeline/replacement.rs"]
mod replacement;

#[path = "wald_pipeline/replacement_metadata.rs"]
mod replacement_metadata;

#[path = "wald_pipeline/wald_options.rs"]
mod wald_options;

#[path = "wald_pipeline/zeros_weights.rs"]
mod zeros_weights;
