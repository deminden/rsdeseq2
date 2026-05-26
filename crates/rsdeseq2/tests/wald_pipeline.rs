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

#[test]
fn fixed_dispersion_wald_pipeline_uses_intercept_shortcut() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        3,
        vec![2, 4, 6, 10, 10, 10],
        Some(vec!["gene_a".into(), "gene_b".into()]),
        None,
    )
    .unwrap();
    let design =
        DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], Some(vec!["Intercept".into()]))
            .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1, 0.2], 0)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap().n_coefficients(), 1);
    assert_eq!(fit.dispersion.as_deref(), Some(&[0.1, 0.2][..]));
    assert_eq!(fit.beta_converged.as_ref().unwrap(), &vec![true, true]);
    assert_eq!(fit.wald.as_ref().unwrap().stat.len(), 2);
    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[0],
        4.0_f64.log2(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        4.0_f64.log2(),
        epsilon = 1e-12
    );
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(fit.hat_diagonal.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.hat_diagonal.as_ref().unwrap().n_cols(), 3);
    assert_wald_likelihood_state(&fit, &counts);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_rows(), 2);
    assert_eq!(beta_covariance.n_cols(), 1);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[0].sqrt(),
        fit.beta_se.as_ref().unwrap().row(0).unwrap()[0],
        epsilon = 1e-12
    );
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.max_cooks.as_ref().unwrap().len(), 2);
    assert!(results.rows[0].max_cooks.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_uses_normalization_factors_for_offsets() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        3,
        vec![10, 20, 40],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 2.0, 4.0]).unwrap();
    let design =
        DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], Some(vec!["Intercept".into()]))
            .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![100.0, 100.0, 100.0])
        .normalization_factors(normalization_factors.clone())
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 0)
        .unwrap();

    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_eq!(fit.size_factors, vec![100.0, 100.0, 100.0]);
    assert_relative_eq!(fit.base_mean[0], 10.0, epsilon = 1e-12);
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[0],
        10.0_f64.log2(),
        epsilon = 1e-12
    );
    for (actual, expected) in fit
        .mu
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .zip([10.0, 20.0, 40.0])
    {
        assert_relative_eq!(*actual, expected, epsilon = 1e-12);
    }
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_wald_likelihood_state(&fit, &counts);
    assert_relative_eq!(results.rows[0].base_mean, 10.0, epsilon = 1e-12);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        10.0_f64.log2(),
        epsilon = 1e-12
    );
}

#[test]
fn fixed_dispersion_wald_pipeline_uses_irls_for_two_group_design() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), 2);
    assert!(fit.beta_converged.as_ref().unwrap()[0]);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_rows(), 1);
    assert_eq!(beta_covariance.n_cols(), 4);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[3].sqrt(),
        fit.beta_se.as_ref().unwrap().row(0).unwrap()[1],
        epsilon = 1e-8
    );
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[1],
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_wald_likelihood_state(&fit, &counts);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(results.rows[0].pvalue, fit.wald.as_ref().unwrap().pvalue[0]);
    assert_eq!(results.rows[0].dispersion, Some(0.05));
    assert_eq!(fit.cooks.as_ref().unwrap().n_cols(), 4);
    assert_eq!(results.rows[0].max_cooks, None);
}

#[test]
fn fixed_dispersion_wald_contrast_matches_selected_coefficient_for_two_group_design() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (coefficient_fit, coefficient_results) = builder
        .clone()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();
    let (contrast_fit, contrast_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.05], &[0.0, 1.0])
        .unwrap();

    assert_wald_fit_intermediates_match(
        &contrast_fit,
        &coefficient_fit,
        "coefficient-equivalent contrast",
    );
    assert_relative_eq!(
        contrast_results.rows[0].log2_fold_change.unwrap(),
        coefficient_results.rows[0].log2_fold_change.unwrap(),
        epsilon = 1e-10
    );
    assert_relative_eq!(
        contrast_results.rows[0].lfc_se.unwrap(),
        coefficient_results.rows[0].lfc_se.unwrap(),
        epsilon = 1e-10
    );
    assert_eq!(
        contrast_fit.wald.as_ref().unwrap().stat[0],
        coefficient_fit.wald.as_ref().unwrap().stat[0]
    );
    assert_eq!(
        contrast_results.rows[0].pvalue,
        coefficient_results.rows[0].pvalue
    );
    assert_eq!(contrast_results.rows[0].dispersion, Some(0.05));
    assert_eq!(contrast_fit.cooks.as_ref().unwrap().n_cols(), 4);
}

#[test]
fn fixed_dispersion_wald_contrast_spec_resolves_coefficient_name() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (fit, results) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let (primitive_fit, _primitive_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.05], &[0.0, 1.0])
        .unwrap();

    assert_wald_fit_intermediates_match(&fit, &primitive_fit, "coefficient-name contrast spec");
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn fixed_dispersion_wald_contrast_spec_resolves_factor_level_shape() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (fit, results) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::factor_level("condition", "B", "A"),
        )
        .unwrap();
    let (primitive_fit, _primitive_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.05], &[0.0, 1.0])
        .unwrap();

    assert_wald_fit_intermediates_match(&fit, &primitive_fit, "factor-level contrast spec");
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn fixed_dispersion_wald_contrast_spec_infers_shared_reference_factor_levels() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        6,
        vec![10, 10, 20, 20, 40, 40],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "condition_C_vs_A".into(),
        ]),
    )
    .unwrap();

    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (fit, results) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::factor_level("condition", "C", "B"),
        )
        .unwrap();
    let (primitive_fit, _primitive_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.05], &[0.0, -1.0, 1.0])
        .unwrap();

    assert_wald_fit_intermediates_match(
        &fit,
        &primitive_fit,
        "shared-reference factor-level contrast spec",
    );
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_C_vs_B")
    );
}

#[test]
fn fixed_dispersion_wald_contrast_spec_resolves_name_lists() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 20, 20, 40],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "batch_Y_vs_X".into(),
        ]),
    )
    .unwrap();

    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (fit, results) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::list(vec!["condition_B_vs_A".into()], vec!["batch_Y_vs_X".into()]),
        )
        .unwrap();
    let (primitive_fit, _primitive_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.05], &[0.0, 1.0, -1.0])
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let expected = beta.row(0).unwrap()[1] - beta.row(0).unwrap()[2];
    assert_wald_fit_intermediates_match(&fit, &primitive_fit, "list contrast spec");
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        expected,
        epsilon = 1e-10
    );
    assert_eq!(results.metadata.result_name.as_deref(), Some("contrast"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: condition_B_vs_A vs batch_Y_vs_X")
    );
}

#[test]
fn original_zero_zero_list_contrast_zeroes_lfc_like_numeric_contrast() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            100, 110, 0, 0, 100, 110, 0, 0, //
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, 1.0, //
            1.0, 0.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "condition_C_vs_A".into(),
            "condition_D_vs_A".into(),
        ]),
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 2.0, 2.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (list_fit, list_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.1, 0.1],
            &ContrastSpec::list(
                vec!["condition_D_vs_A".into()],
                vec!["condition_B_vs_A".into()],
            ),
        )
        .unwrap();
    let (_numeric_fit, numeric_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1, 0.1], &[0.0, -1.0, 0.0, 1.0])
        .unwrap();

    assert_eq!(list_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(numeric_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(list_results.rows[0].pvalue, numeric_results.rows[0].pvalue);
    assert_eq!(list_results.rows[1].log2_fold_change, None);
    assert_eq!(
        list_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        list_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: condition_D_vs_A vs condition_B_vs_A")
    );
    assert_eq!(
        list_fit.wald.as_ref().unwrap().stat,
        numeric_results
            .rows
            .iter()
            .map(|row| row.stat)
            .collect::<Vec<_>>()
    );
}

#[test]
fn original_zero_intercept_factor_level_contrasts_return_signed_lfcs() {
    let counts = CountMatrix::from_row_major_u32(
        1,
        12,
        vec![100, 100, 100, 100, 200, 200, 200, 200, 400, 400, 400, 400],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        12,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 1.0,
        ],
        Some(vec![
            "condition1".into(),
            "condition2".into(),
            "condition3".into(),
        ]),
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 12])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        });

    let (_fit_21, result_21) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::factor_level("condition", "2", "1"),
        )
        .unwrap();
    let (_fit_32, result_32) = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::factor_level("condition", "3", "2"),
        )
        .unwrap();
    let (_fit_13, result_13) = builder
        .fit_fixed_dispersion_wald_contrast_spec(
            &counts,
            &design,
            &[0.05],
            &ContrastSpec::factor_level("condition", "1", "3"),
        )
        .unwrap();

    assert_relative_eq!(
        result_21.rows[0].log2_fold_change.unwrap(),
        1.0,
        epsilon = 1e-8
    );
    assert_relative_eq!(
        result_32.rows[0].log2_fold_change.unwrap(),
        1.0,
        epsilon = 1e-8
    );
    assert_relative_eq!(
        result_13.rows[0].log2_fold_change.unwrap(),
        -2.0,
        epsilon = 1e-8
    );
    assert_eq!(
        result_21.metadata.result_name.as_deref(),
        Some("condition_2_vs_1")
    );
}

#[test]
fn fixed_dispersion_wald_contrast_applies_explicit_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.2], &[1.0])
        .unwrap();

    assert_relative_eq!(
        results.rows[0].max_cooks.unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
}

#[test]
fn fixed_dispersion_wald_contrast_applies_contrast_all_zero_numeric() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        6,
        vec![
            0, 0, 0, 0, 50, 60, //
            10, 12, 30, 36, 50, 60,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 1.0,
        ],
        Some(vec![
            "conditionA".into(),
            "conditionB".into(),
            "conditionC".into(),
        ]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1, 0.1], &[-1.0, 1.0, 0.0])
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    assert_eq!(results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(results.rows[0].stat, Some(0.0));
    assert_eq!(results.rows[0].pvalue, Some(1.0));
    assert_eq!(wald.stat[0], Some(0.0));
    assert_eq!(wald.pvalue[0], Some(1.0));
    assert!(results.rows[0].lfc_se.is_some());
    assert!(results.rows[1].pvalue.is_some());
    assert_ne!(results.rows[1].pvalue, Some(1.0));
}

#[test]
fn fixed_dispersion_wald_factor_level_contrast_applies_character_contrast_all_zero() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        6,
        vec![
            0, 0, 0, 0, 50, 60, //
            10, 12, 30, 36, 50, 60,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "condition_C_vs_A".into(),
        ]),
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "C", "C"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (fit, results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (primitive_fit, _primitive_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1, 0.1], &[0.0, 1.0, 0.0])
        .unwrap();

    assert_wald_fit_intermediates_match(&fit, &primitive_fit, "factor-level contrast helper");
    assert_eq!(fit.wald.as_ref().unwrap().stat[0], Some(0.0));
    assert_eq!(fit.wald.as_ref().unwrap().pvalue[0], Some(1.0));
    assert_eq!(results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(results.rows[0].stat, Some(0.0));
    assert_eq!(results.rows[0].pvalue, Some(1.0));
    assert!(results.rows[0].lfc_se.is_some());
    assert!(results.rows[1].pvalue.is_some());
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    assert_ne!(results.rows[1].pvalue, Some(1.0));
}

#[test]
fn original_weighted_contrast_lfc_se_does_not_depend_on_contrast_type() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            10, 12, 30, 36, //
            20, 25, 40, 48,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 0.8, 1.0, 0.7, //
            0.9, 1.0, 0.85, 1.0,
        ],
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 4])
        .observation_weights(weights)
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_factor_fit, factor_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (_numeric_fit, numeric_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1, 0.1], &[0.0, 1.0])
        .unwrap();

    for (factor, numeric) in factor_results.rows.iter().zip(numeric_results.rows.iter()) {
        assert_eq!(factor.log2_fold_change, numeric.log2_fold_change);
        assert_eq!(factor.lfc_se, numeric.lfc_se);
        assert_eq!(factor.stat, numeric.stat);
        assert_eq!(factor.pvalue, numeric.pvalue);
    }
}

#[test]
fn fixed_dispersion_wald_contrast_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[], &[1.0])
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1], &[1.0, 0.0])
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1], &[0.0])
        .is_err());
}

#[test]
fn fixed_dispersion_wald_contrast_marks_weight_failed_rows_as_skipped() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            10, 20, 30, 40, //
            50, 60, 70, 80,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 0.0, 0.0,
        ],
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .observation_weights(weights)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1, 0.1], &[0.0, 1.0])
        .unwrap();

    assert_eq!(fit.weights_fail, Some(vec![false, true]));
    assert_eq!(fit.all_zero, vec![false, true]);
    assert!(results.rows[0].pvalue.is_some());
    assert_eq!(results.rows[1].log2_fold_change, None);
    assert_eq!(results.rows[1].pvalue, None);
    assert!(fit.beta.as_ref().unwrap().row(1).unwrap()[0].is_nan());
}

#[test]
fn original_zero_weighted_sample_matches_removed_sample_fit() {
    let weighted_counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 12, 80, 120]).unwrap();
    let weighted_design =
        DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 1.0, 0.0, 1.0]).unwrap();

    let (weighted_fit, weighted_results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .observation_weights(weights)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&weighted_counts, &weighted_design, &[0.05], 0)
        .unwrap();

    let subset_counts = CountMatrix::from_row_major_u32(1, 3, vec![10, 12, 120]).unwrap();
    let subset_design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();
    let (subset_fit, subset_results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&subset_counts, &subset_design, &[0.05], 0)
        .unwrap();

    assert_eq!(weighted_fit.weights_fail, Some(vec![false]));
    assert_relative_eq!(
        weighted_fit.beta.as_ref().unwrap().row(0).unwrap()[0],
        subset_fit.beta.as_ref().unwrap().row(0).unwrap()[0],
        epsilon = 1e-10
    );
    assert_relative_eq!(
        weighted_fit.beta_se.as_ref().unwrap().row(0).unwrap()[0],
        subset_fit.beta_se.as_ref().unwrap().row(0).unwrap()[0],
        epsilon = 1e-10
    );
    assert_relative_eq!(
        weighted_fit.log_like.as_ref().unwrap()[0],
        subset_fit.log_like.as_ref().unwrap()[0],
        epsilon = 1e-10
    );
    assert_eq!(
        weighted_results.rows[0].pvalue,
        subset_results.rows[0].pvalue
    );
}

#[test]
fn fixed_dispersion_wald_pipeline_can_use_qr_irls_solver() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            solver: IrlsSolver::Qr,
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    assert!(fit.beta_converged.as_ref().unwrap()[0]);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(results.rows[0].pvalue, fit.wald.as_ref().unwrap().pvalue[0]);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_use_t_pvalues() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .wald_t_degrees_of_freedom(4.0)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    assert_eq!(wald.degrees_of_freedom.as_ref().unwrap(), &vec![Some(4.0)]);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(wald.stat[0].unwrap(), 4.0).unwrap(),
        epsilon = 1e-15
    );
    assert_eq!(results.rows[0].pvalue, wald.pvalue[0]);
    assert_eq!(results.rows[0].padj, wald.pvalue[0]);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_use_lfc_threshold() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .wald_lfc_threshold(0.5, WaldAlternative::Greater)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    assert!(wald.stat[0].unwrap() >= 0.0);
    assert_eq!(results.rows[0].pvalue, wald.pvalue[0]);
    assert_eq!(results.rows[0].padj, wald.pvalue[0]);
    assert!(wald.pvalue[0].unwrap() > 0.0);
    assert!(wald.pvalue[0].unwrap() < 0.5);
    assert_eq!(results.metadata.lfc_threshold, 0.5);
    assert_eq!(results.metadata.alt_hypothesis.as_deref(), Some("greater"));
}

#[test]
fn fixed_dispersion_wald_pipeline_validates_per_gene_t_df_length() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 6, 10, 10, 10]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let err = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .wald_t_per_gene_degrees_of_freedom(vec![4.0])
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1, 0.2], 0)
        .unwrap_err();

    assert!(err.to_string().contains("degrees of freedom"));
}

#[test]
fn fixed_dispersion_wald_pipeline_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();
    let rank_deficient = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald(&counts, &design, &[], 0)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_wald(&counts, &rank_deficient, &[0.1], 1)
        .is_err());
}

#[test]
fn fixed_dispersion_wald_pipeline_applies_explicit_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2], 0)
        .unwrap();

    assert_relative_eq!(
        results.rows[0].max_cooks.unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_disable_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .disable_cooks_cutoff()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2], 0)
        .unwrap();

    assert_eq!(results.rows[0].cooks_outlier, None);
    assert!(results.rows[0].pvalue.is_some());
    assert!(results.rows[0].padj.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_applies_independent_filtering() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 1, 1, 100, 100, 100]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .independent_filtering_theta(vec![1.0, 1.0])
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2, 0.2], 0)
        .unwrap();

    let filtering = results.independent_filtering.as_ref().unwrap();
    assert!(filtering.enabled);
    assert_eq!(filtering.filter_threshold, Some(100.0));
    assert_eq!(results.rows[0].filtered, Some(true));
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[1].filtered, Some(false));
    assert!(results.rows[1].padj.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_can_disable_independent_filtering() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 1, 1, 100, 100, 100]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2, 0.2], 0)
        .unwrap();

    let filtering = results.independent_filtering.as_ref().unwrap();
    assert!(!filtering.enabled);
    assert_eq!(results.rows[0].filtered, None);
    assert_eq!(results.rows[1].filtered, None);
    assert!(results.rows[0].padj.is_some());
    assert!(results.rows[1].padj.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_expands_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        3,
        vec![0, 0, 0, 2, 4, 6],
        Some(vec!["zero_gene".into(), "signal_gene".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .wald_t_residual_degrees_of_freedom()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1, 0.2], 0)
        .unwrap();

    assert_eq!(fit.all_zero, vec![true, false]);
    assert!(fit.dispersion.as_ref().unwrap()[0].is_nan());
    assert_eq!(fit.dispersion.as_ref().unwrap()[1], 0.2);
    assert!(fit.beta.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().row(1).unwrap()[0],
        4.0_f64.log2(),
        epsilon = 1e-12
    );
    assert!(fit.beta_covariance.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.beta_covariance.as_ref().unwrap().row(1).unwrap()[0].is_finite());
    assert!(fit.mu.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_eq!(fit.beta_converged.as_ref().unwrap(), &vec![false, true]);
    assert_eq!(fit.wald.as_ref().unwrap().pvalue[0], None);
    assert!(fit.wald.as_ref().unwrap().pvalue[1].is_some());
    assert_eq!(
        fit.wald
            .as_ref()
            .unwrap()
            .degrees_of_freedom
            .as_ref()
            .unwrap(),
        &vec![None, Some(2.0)]
    );
    assert!(fit.cooks.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_relative_eq!(
        fit.cooks.as_ref().unwrap().row(1).unwrap()[0],
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
    assert_eq!(fit.max_cooks.as_ref().unwrap()[0], None);
    assert_relative_eq!(
        fit.max_cooks.as_ref().unwrap()[1].unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );

    assert_eq!(results.rows[0].gene.as_deref(), Some("zero_gene"));
    assert_eq!(results.rows[0].base_mean, 0.0);
    assert_eq!(results.rows[0].log2_fold_change, None);
    assert_eq!(results.rows[0].lfc_se, None);
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].dispersion, None);
    assert_eq!(results.rows[0].converged, None);
    assert_eq!(results.rows[0].max_cooks, None);
    assert_eq!(results.rows[0].cooks_outlier, None);
    assert_eq!(results.rows[1].gene.as_deref(), Some("signal_gene"));
    assert_eq!(results.rows[1].dispersion, Some(0.2));
    assert_relative_eq!(
        results.rows[1].max_cooks.unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
}

#[test]
fn fixed_dispersion_wald_pipeline_handles_all_zero_only_with_supplied_size_factors() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![0, 0, 0]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 0)
        .unwrap();

    assert_eq!(fit.base_mean, vec![0.0]);
    assert!(fit.beta.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.beta_se.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.beta_covariance.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.mu.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.hat_diagonal.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.cooks.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_eq!(fit.max_cooks.as_ref().unwrap()[0], None);
    assert_eq!(fit.wald.as_ref().unwrap().stat[0], None);
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].converged, None);
    assert_eq!(results.rows[0].max_cooks, None);
    assert_eq!(results.rows[0].cooks_outlier, None);
}
