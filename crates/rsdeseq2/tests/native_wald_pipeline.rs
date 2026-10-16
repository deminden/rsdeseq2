//! Integration coverage grouped by workflow; shared fixtures remain in this test crate.

use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn two_group_design() -> DesignMatrix {
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

fn two_group_design_with_r_cleaned_coefficient() -> DesignMatrix {
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

fn native_wald_counts_with_zero_row() -> CountMatrix {
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

fn native_wald_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
}

fn glm_mu_native_wald_builder() -> DeseqBuilder {
    native_wald_builder().gene_wise_dispersion_options(GeneWiseDispersionOptions {
        use_cox_reid: false,
        fit_method: GeneWiseDispersionFitMethod::Grid,
        niter: 2,
        ..GeneWiseDispersionOptions::default()
    })
}

fn unit_weights_for(counts: &CountMatrix) -> RowMajorMatrix<f64> {
    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 1.0).unwrap()
}

fn nonunit_weights_for(counts: &CountMatrix) -> RowMajorMatrix<f64> {
    let sample_pattern = [1.0, 0.85, 1.0, 0.8, 1.0, 0.75, 1.0, 0.7];
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let gene_scale = 1.0 - 0.03 * (gene % 3) as f64;
        for weight in sample_pattern.iter().take(counts.n_samples()) {
            values.push(*weight * gene_scale);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values).unwrap()
}

fn assert_float_close_or_nan(actual: f64, expected: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
        return;
    }
    let diff = (actual - expected).abs();
    let allowed = 1e-9 + 1e-9 * expected.abs().max(1.0);
    assert!(
        diff <= allowed,
        "{label}: actual={actual}, expected={expected}, diff={diff}, allowed={allowed}"
    );
}

fn assert_optional_float_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close_or_nan(actual, expected, label),
        _ => assert_eq!(actual, expected, "{label}: optional value mismatch"),
    }
}

fn assert_result_rows_close(actual: &[DeseqResultRow], expected: &[DeseqResultRow], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: row count mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        let row_label = format!("{label} row {index}");
        assert_eq!(actual.gene, expected.gene, "{row_label}: gene mismatch");
        assert_float_close_or_nan(
            actual.base_mean,
            expected.base_mean,
            &format!("{row_label} baseMean"),
        );
        assert_optional_float_close(
            actual.log2_fold_change,
            expected.log2_fold_change,
            &format!("{row_label} log2FoldChange"),
        );
        assert_optional_float_close(
            actual.lfc_se,
            expected.lfc_se,
            &format!("{row_label} lfcSE"),
        );
        assert_optional_float_close(actual.stat, expected.stat, &format!("{row_label} stat"));
        assert_optional_float_close(
            actual.pvalue,
            expected.pvalue,
            &format!("{row_label} pvalue"),
        );
        assert_optional_float_close(actual.padj, expected.padj, &format!("{row_label} padj"));
        assert_optional_float_close(
            actual.dispersion,
            expected.dispersion,
            &format!("{row_label} dispersion"),
        );
        assert_eq!(
            actual.converged, expected.converged,
            "{row_label}: converged mismatch"
        );
        assert_optional_float_close(
            actual.max_cooks,
            expected.max_cooks,
            &format!("{row_label} maxCooks"),
        );
        assert_eq!(
            actual.cooks_outlier, expected.cooks_outlier,
            "{row_label}: cooks outlier mismatch"
        );
        assert_eq!(
            actual.filtered, expected.filtered,
            "{row_label}: filtered mismatch"
        );
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

fn assert_wald_fit_state_matches(actual: &DeseqFit, expected: &DeseqFit, label: &str) {
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
        &format!("{label} log_like"),
    );
    assert_slice_close_or_nan(
        actual.full_deviance.as_ref().unwrap(),
        expected.full_deviance.as_ref().unwrap(),
        &format!("{label} full_deviance"),
    );
    assert_matrix_close_or_nan(
        actual.mu.as_ref().unwrap(),
        expected.mu.as_ref().unwrap(),
        &format!("{label} fitted mu"),
    );
    assert_matrix_close_or_nan(
        actual.hat_diagonal.as_ref().unwrap(),
        expected.hat_diagonal.as_ref().unwrap(),
        &format!("{label} hat diagonal"),
    );
}

fn assert_option_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close_or_nan(actual, expected, label),
        (None, None) => {}
        _ => panic!("{label}: actual={actual:?}, expected={expected:?}"),
    }
}

fn assert_option_slice_close(actual: &[Option<f64>], expected: &[Option<f64>], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_option_close(*actual, *expected, &format!("{label}[{index}]"));
    }
}

#[path = "native_wald_pipeline/entry_points.rs"]
mod entry_points;

#[path = "native_wald_pipeline/stored_contrasts.rs"]
mod stored_contrasts;

#[path = "native_wald_pipeline/replacement.rs"]
mod replacement;

#[path = "native_wald_pipeline/cooks_gates.rs"]
mod cooks_gates;

#[path = "native_wald_pipeline/test_dispatch.rs"]
mod test_dispatch;

#[path = "native_wald_pipeline/lrt_formula.rs"]
mod lrt_formula;

#[path = "native_wald_pipeline/linear_fits.rs"]
mod linear_fits;

#[path = "native_wald_pipeline/contrast_requests.rs"]
mod contrast_requests;

#[path = "native_wald_pipeline/dispersion_fits.rs"]
mod dispersion_fits;

#[path = "native_wald_pipeline/weights_and_local_trends.rs"]
mod weights_and_local_trends;

#[path = "native_wald_pipeline/transforms.rs"]
mod transforms;

#[path = "native_wald_pipeline/fast_vst.rs"]
mod fast_vst;
