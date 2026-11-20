//! Numerical coverage grouped by responsibility; shared fixtures stay in this test crate.

use crate::support as common;
use approx::assert_relative_eq;
use common::{assert_float_close_or_nan, assert_matrix_close_or_nan, assert_slice_close_or_nan};
use rsdeseq2::prelude::*;

fn full_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap()
}

fn full_design_with_r_cleaned_coefficient() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition.B.1".into()]),
    )
    .unwrap()
}

fn reduced_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap()
}

fn counts_with_zero_row() -> CountMatrix {
    CountMatrix::from_row_major_u32_with_names(
        8,
        8,
        vec![
            0, 0, 0, 0, 0, 0, 0, 0, //
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29, //
            30, 50, 25, 55, 35, 45, 28, 52, //
            55, 105, 60, 100, 50, 110, 65, 95, //
            120, 200, 130, 190, 115, 205, 125, 195, //
            240, 400, 250, 390, 230, 410, 260, 380, //
            15, 18, 12, 17, 45, 50, 40, 55,
        ],
        Some(
            [
                "zero", "up", "flat", "variable", "high_up", "stable", "low_up", "broad",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap()
}

fn native_lrt_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 8])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
}

fn assert_option_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close_or_nan(actual, expected, label),
        (None, None) => {}
        _ => panic!("{label}: actual={actual:?}, expected={expected:?}"),
    }
}

fn assert_lrt_fit_state_matches(actual: &DeseqFit, expected: &DeseqFit, label: &str) {
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
    assert_slice_close_or_nan(
        actual.log_like.as_ref().unwrap(),
        expected.log_like.as_ref().unwrap(),
        &format!("{label} full log_like"),
    );
    assert_slice_close_or_nan(
        actual.full_deviance.as_ref().unwrap(),
        expected.full_deviance.as_ref().unwrap(),
        &format!("{label} full_deviance"),
    );
    assert_slice_close_or_nan(
        actual.reduced_log_like.as_ref().unwrap(),
        expected.reduced_log_like.as_ref().unwrap(),
        &format!("{label} reduced log_like"),
    );
    assert_matrix_close_or_nan(
        actual.mu.as_ref().unwrap(),
        expected.mu.as_ref().unwrap(),
        &format!("{label} full mu"),
    );
    assert_matrix_close_or_nan(
        actual.hat_diagonal.as_ref().unwrap(),
        expected.hat_diagonal.as_ref().unwrap(),
        &format!("{label} full hat"),
    );
    assert_matrix_close_or_nan(
        actual.reduced_mu.as_ref().unwrap(),
        expected.reduced_mu.as_ref().unwrap(),
        &format!("{label} reduced mu"),
    );
    assert_matrix_close_or_nan(
        actual.reduced_hat_diagonal.as_ref().unwrap(),
        expected.reduced_hat_diagonal.as_ref().unwrap(),
        &format!("{label} reduced hat"),
    );
}

#[path = "native_lrt_pipeline/contrasts.rs"]
mod contrasts;

#[path = "native_lrt_pipeline/fitted_state.rs"]
mod fitted_state;

#[path = "native_lrt_pipeline/fitting.rs"]
mod fitting;

#[path = "native_lrt_pipeline/formula_aliases.rs"]
mod formula_aliases;

#[path = "native_lrt_pipeline/formula_designs.rs"]
mod formula_designs;

#[path = "native_lrt_pipeline/replacement.rs"]
mod replacement;

#[path = "native_lrt_pipeline/replacement_formula.rs"]
mod replacement_formula;
