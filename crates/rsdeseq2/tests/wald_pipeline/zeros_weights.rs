use super::*;

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
        Some("condition_X2_vs_X1")
    );
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
