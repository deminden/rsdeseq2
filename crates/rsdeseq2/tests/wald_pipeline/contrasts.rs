use super::*;

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
fn fixed_dispersion_wald_results_contrast_routes_deseq2_contrast_forms() {
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

    let (_fit, character_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
        )
        .unwrap();
    assert_eq!(character_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(character_results.rows[0].pvalue, Some(1.0));
    assert_eq!(
        character_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (_fit, numeric_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::numeric(vec![0.0, 1.0, 0.0]),
            None::<&[String]>,
        )
        .unwrap();
    assert!(numeric_results.rows[0].pvalue.is_some());
    assert_ne!(numeric_results.rows[0].pvalue, Some(1.0));

    let (_fit, list_results) = builder
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::list(
                vec!["condition_C_vs_A".into()],
                vec!["condition_B_vs_A".into()],
            ),
            None::<&[String]>,
        )
        .unwrap();
    assert_eq!(
        list_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        list_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: condition_C_vs_A vs condition_B_vs_A")
    );
}

#[test]
fn fixed_dispersion_wald_results_character_contrast_requires_sample_levels() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![1, 2]).unwrap();
    let design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let err = DeseqBuilder::new()
        .size_factors(vec![1.0; 2])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1],
            &ResultsContrast::character("condition", "B", "A"),
            None::<&[String]>,
        )
        .unwrap_err();

    assert!(
        err.to_string()
            .contains("requires sample levels for contrastAllZero")
    );
}

#[test]
fn fixed_dispersion_wald_contrast_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald_contrast(&counts, &design, &[], &[1.0])
            .is_err()
    );
    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1], &[1.0, 0.0])
            .is_err()
    );
    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.1], &[0.0])
            .is_err()
    );
}
