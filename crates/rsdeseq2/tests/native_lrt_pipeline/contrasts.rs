use super::*;

#[test]
fn native_glm_mu_lrt_contrast_keeps_lrt_pvalues_and_reports_contrast_effect() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let contrast = [0.0, 1.0];

    let (coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let (contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();
    let invalid = ResultRequest::new(
        ResultsContrast::Numeric(contrast.to_vec()),
        WaldTestOptions::t_degrees_of_freedom(6.0),
        CooksCutoff::Disabled,
        IndependentFilteringOptions::default(),
    )
    .unwrap();
    assert!(
        contrast_fit
            .results_for_request(&counts, &invalid, None)
            .unwrap_err()
            .to_string()
            .contains("cannot change a stored LRT")
    );
    let final_dispersions = coefficient_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_contrast_fit, fixed_contrast_results) = builder
        .fit_fixed_dispersion_lrt_contrast(&counts, &full, &reduced, &final_dispersions, &contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&contrast_fit, &coefficient_fit, "LRT contrast");
    assert_lrt_fit_state_matches(
        &fixed_contrast_fit,
        &coefficient_fit,
        "fixed-dispersion LRT contrast",
    );
    assert_eq!(
        contrast_results.metadata.test_type,
        Some(TestType::Lrt),
        "contrast result remains an LRT table"
    );
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    for gene in 0..counts.n_genes() {
        assert_eq!(
            contrast_results.rows[gene].stat,
            coefficient_results.rows[gene].stat
        );
        assert_eq!(
            contrast_results.rows[gene].pvalue,
            coefficient_results.rows[gene].pvalue
        );
        assert_eq!(
            contrast_results.rows[gene].padj,
            coefficient_results.rows[gene].padj
        );
        assert_option_close(
            contrast_results.rows[gene].log2_fold_change,
            coefficient_results.rows[gene].log2_fold_change,
            &format!("LRT contrast LFC gene {gene}"),
        );
        assert_option_close(
            contrast_results.rows[gene].lfc_se,
            coefficient_results.rows[gene].lfc_se,
            &format!("LRT contrast SE gene {gene}"),
        );
        assert_eq!(
            fixed_contrast_results.rows[gene], contrast_results.rows[gene],
            "fixed/native LRT contrast row {gene}"
        );
    }
}

#[test]
fn lrt_results_contrast_request_routes_character_list_and_numeric_forms() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().automatic_dispersion_means();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"];
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.iter().map(|level| level.to_string()).collect(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let (coefficient_fit, _coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = coefficient_fit.dispersion.as_ref().unwrap().clone();

    let (_character_fit, character_results) = builder
        .fit_lrt_with_results_contrast_request(
            &counts,
            &full,
            &reduced,
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
        )
        .unwrap();
    let (_model_frame_fit, model_frame_results) = builder
        .fit_lrt_with_results_contrast_request_from_model_frame(
            &counts,
            &full,
            &reduced,
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
        )
        .unwrap();
    let (_fixed_character_fit, fixed_character_results) = builder
        .fit_fixed_dispersion_lrt_results_contrast(
            &counts,
            &full,
            &reduced,
            &final_dispersions,
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
        )
        .unwrap();
    let (_list_fit, list_results) = builder
        .fit_lrt_with_results_contrast_request(
            &counts,
            &full,
            &reduced,
            &ResultsContrast::list(vec!["condition_B_vs_A".into()], vec![]),
            Option::<&[&str]>::None,
        )
        .unwrap();
    let (_numeric_fit, numeric_results) = builder
        .clone()
        .test(TestType::Lrt)
        .reduced_design(reduced.clone())
        .fit_with_results_contrast_request(
            &counts,
            &full,
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            Option::<&[&str]>::None,
        )
        .unwrap();

    assert_eq!(character_results.rows, fixed_character_results.rows);
    assert_eq!(model_frame_results, character_results);
    assert_eq!(character_results.rows, list_results.rows);
    assert_eq!(character_results.rows, numeric_results.rows);
    assert_eq!(
        character_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        list_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        numeric_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    let missing_levels = builder.fit_lrt_with_results_contrast_request(
        &counts,
        &full,
        &reduced,
        &ResultsContrast::character("condition", "B", "A"),
        Option::<&[&str]>::None,
    );
    assert!(matches!(
        missing_levels,
        Err(DeseqError::InvalidOptions { .. })
    ));
}

#[test]
fn lrt_results_contrast_request_canonicalizes_cleaned_sample_level_aliases() {
    let counts = counts_with_zero_row();
    let full = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "cell.type_B.cell_vs_A.cell".into(),
        ]),
    )
    .unwrap();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let levels = [
        "A cell", "A cell", "A cell", "A cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();

    let (_factor_fit, factor_results) = builder
        .fit_lrt_with_results_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("cell type", "B cell", "A cell", &levels),
        )
        .unwrap();
    let (_character_fit, character_results) = builder
        .fit_lrt_with_results_contrast_request(
            &counts,
            &full,
            &reduced,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&levels),
        )
        .unwrap();
    let replacement_options = CooksReplacementOptions::new(f64::MAX);
    let factor_replacement = builder
        .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("cell type", "B cell", "A cell", &levels),
            &replacement_options,
        )
        .unwrap();
    let request_replacement = builder
        .fit_lrt_with_results_contrast_request_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&levels),
            &replacement_options,
        )
        .unwrap();

    assert_eq!(character_results.rows, factor_results.rows);
    assert_eq!(
        character_results.metadata.result_name.as_deref(),
        Some("cell.type_B.cell_vs_A.cell")
    );
    assert_eq!(
        character_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: cell.type B cell vs A cell")
    );
    assert_eq!(
        request_replacement.results.rows,
        factor_replacement.results.rows
    );
    assert_eq!(
        request_replacement.results.metadata.result_name.as_deref(),
        Some("cell.type_B.cell_vs_A.cell")
    );
    assert_eq!(
        request_replacement.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: cell.type B cell vs A cell")
    );

    let ambiguous_levels = [
        "A cell", "A-cell", "A cell", "A-cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let ambiguous = builder.fit_lrt_with_results_contrast_request(
        &counts,
        &full,
        &reduced,
        &ResultsContrast::character("cell.type", "B.cell", "A.cell"),
        Some(&ambiguous_levels),
    );
    assert!(
        ambiguous
            .unwrap_err()
            .to_string()
            .contains("denominator level 'A.cell' resolves ambiguously")
    );
}

#[test]
fn native_linear_mu_lrt_contrast_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let contrast = [0.0, 1.0];

    let (coefficient_fit, coefficient_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let (contrast_fit, contrast_results) = builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();
    let final_dispersions = coefficient_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_contrast_fit, fixed_contrast_results) = builder
        .fit_fixed_dispersion_lrt_contrast(&counts, &full, &reduced, &final_dispersions, &contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&contrast_fit, &coefficient_fit, "linear-mu LRT contrast");
    assert_lrt_fit_state_matches(
        &fixed_contrast_fit,
        &coefficient_fit,
        "linear-mu fixed LRT contrast",
    );
    assert_eq!(contrast_results.rows, fixed_contrast_results.rows);
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    for gene in 0..counts.n_genes() {
        assert_eq!(
            contrast_results.rows[gene].stat,
            coefficient_results.rows[gene].stat
        );
        assert_eq!(
            contrast_results.rows[gene].pvalue,
            coefficient_results.rows[gene].pvalue
        );
        assert_eq!(
            contrast_results.rows[gene].padj,
            coefficient_results.rows[gene].padj
        );
    }
}

#[test]
fn native_linear_mu_lrt_contrast_specs_set_metadata_and_factor_levels() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };

    let (numeric_fit, numeric_results) = builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .fit_lrt_linear_mu_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    let (factor_fit, factor_results) = builder
        .fit_lrt_linear_mu_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&named_fit, &numeric_fit, "linear-mu named LRT contrast");
    assert_lrt_fit_state_matches(
        &factor_fit,
        &numeric_fit,
        "linear-mu factor-level LRT contrast",
    );
    assert_eq!(named_results.rows, numeric_results.rows);
    assert_eq!(factor_results.rows, numeric_results.rows);
    assert_eq!(
        named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_lrt_parametric_contrast_helpers_ignore_builder_fit_type() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let mean_builder = native_lrt_builder().fit_type(FitType::Mean);
    let parametric_builder = mean_builder.clone().fit_type(FitType::Parametric);
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };

    let (linear_parametric_fit, linear_parametric_results) = mean_builder
        .fit_lrt_linear_mu_contrast_parametric(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (linear_expected_fit, linear_expected_results) = parametric_builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_parametric_fit,
        &linear_expected_fit,
        "linear-mu parametric LRT contrast",
    );
    assert_eq!(linear_parametric_results, linear_expected_results);

    let (linear_named_fit, linear_named_results) = mean_builder
        .fit_lrt_linear_mu_contrast_spec_parametric(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_named_fit,
        &linear_expected_fit,
        "linear-mu named parametric LRT contrast",
    );
    assert_eq!(linear_named_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (linear_factor_fit, linear_factor_results) = mean_builder
        .fit_lrt_linear_mu_factor_level_contrast_parametric(
            &counts,
            &full,
            &reduced,
            factor_contrast,
        )
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_factor_fit,
        &linear_expected_fit,
        "linear-mu factor parametric LRT contrast",
    );
    assert_eq!(linear_factor_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );

    let (glm_parametric_fit, glm_parametric_results) = mean_builder
        .fit_lrt_glm_mu_contrast_parametric(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (glm_expected_fit, glm_expected_results) = parametric_builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_parametric_fit,
        &glm_expected_fit,
        "GLM-mu parametric LRT contrast",
    );
    assert_eq!(glm_parametric_results, glm_expected_results);

    let (glm_named_fit, glm_named_results) = mean_builder
        .fit_lrt_glm_mu_contrast_spec_parametric(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_named_fit,
        &glm_expected_fit,
        "GLM-mu named parametric LRT contrast",
    );
    assert_eq!(glm_named_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (glm_factor_fit, glm_factor_results) = mean_builder
        .fit_lrt_glm_mu_factor_level_contrast_parametric(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_factor_fit,
        &glm_expected_fit,
        "GLM-mu factor parametric LRT contrast",
    );
    assert_eq!(glm_factor_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_glm_mu_lrt_contrast_all_zero_only_zeroes_lfc() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        5,
        6,
        vec![
            0, 0, 0, 0, 0, 0, //
            20, 22, 0, 0, 0, 0, //
            12, 28, 14, 26, 11, 29, //
            55, 105, 50, 110, 65, 95, //
            120, 200, 115, 205, 125, 195,
        ],
        Some(
            [
                "zero",
                "contrast_groups_zero",
                "variable",
                "high_up",
                "stable",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
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
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_D_vs_A".to_string(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".to_string()]),
    )
    .unwrap();
    let builder = native_lrt_builder().size_factors(vec![1.0; 6]);
    let contrast = [0.0, -1.0, 1.0];

    let (_coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 2).unwrap();
    let (_contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();

    assert_eq!(
        contrast_results.rows[1].gene.as_deref(),
        Some("contrast_groups_zero")
    );
    assert_eq!(contrast_results.rows[1].log2_fold_change, Some(0.0));
    assert_eq!(
        contrast_results.rows[1].stat,
        coefficient_results.rows[1].stat
    );
    assert_eq!(
        contrast_results.rows[1].pvalue,
        coefficient_results.rows[1].pvalue
    );
    assert_eq!(
        contrast_results.rows[1].padj,
        coefficient_results.rows[1].padj
    );
    assert!(contrast_results.rows[1].stat.unwrap().is_finite());
    assert!(contrast_results.rows[1].pvalue.unwrap().is_finite());
}

#[test]
fn native_glm_mu_lrt_factor_level_contrast_uses_character_all_zero_cleanup() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        5,
        6,
        vec![
            0, 0, 0, 0, 0, 0, //
            20, 22, 0, 0, 0, 0, //
            12, 28, 14, 26, 11, 29, //
            55, 105, 50, 110, 65, 95, //
            120, 200, 115, 205, 125, 195,
        ],
        Some(
            [
                "zero",
                "contrast_groups_zero",
                "variable",
                "high_up",
                "stable",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
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
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_D_vs_A".to_string(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".to_string()]),
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "D", "D"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let builder = native_lrt_builder().size_factors(vec![1.0; 6]);
    let contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "D",
        denominator: "B",
        reference: Some("A"),
        sample_levels: &levels,
    };

    let (_coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 2).unwrap();
    let (_contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_factor_level_contrast(&counts, &full, &reduced, contrast)
        .unwrap();

    assert_eq!(
        contrast_results.rows[1].gene.as_deref(),
        Some("contrast_groups_zero")
    );
    assert_eq!(contrast_results.rows[1].log2_fold_change, Some(0.0));
    assert_eq!(
        contrast_results.rows[1].stat,
        coefficient_results.rows[1].stat
    );
    assert_eq!(
        contrast_results.rows[1].pvalue,
        coefficient_results.rows[1].pvalue
    );
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("condition_D_vs_B")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition D vs B")
    );
}

#[test]
fn native_glm_mu_lrt_factor_level_contrast_applies_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
        6,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "A", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();
    let (_numeric_fit, numeric_results) = builder
        .clone()
        .fit_lrt_with_results_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (_fit, results) = builder
        .fit_lrt_with_results_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(numeric_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].cooks_outlier, Some(false));
}

#[test]
fn top_level_lrt_contrast_spec_resolves_named_full_model_effect() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");

    let (numeric_fit, numeric_results) = builder
        .fit_lrt_with_results_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .fit_lrt_with_results_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();

    assert_lrt_fit_state_matches(&named_fit, &numeric_fit, "named LRT contrast");
    assert_eq!(named_results.rows, numeric_results.rows);
    assert_eq!(
        named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}
