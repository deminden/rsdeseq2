use super::*;

#[test]
fn top_level_lrt_formula_contrast_request_uses_cleaned_stored_model_frame_factor_name() {
    let counts = counts_with_zero_row();
    let reduced = reduced_design();
    let levels = ["A 0", "A 0", "A 0", "A 0", "B-1", "B-1", "B-1", "B-1"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels.clone(),
            levels: Some(vec!["A 0".to_string(), "B-1".to_string()]),
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let builder = native_lrt_builder().model_frame(model_frame);
    let full = builder
        .expanded_formula_design("~ `cell type`")
        .unwrap()
        .standard_design;
    let options = CooksReplacementOptions::new(f64::MAX);

    let (_explicit_fit, explicit_results) = builder
        .fit_lrt_with_results_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("cell type", "B-1", "A 0", &levels),
        )
        .unwrap();
    let (formula_fit, formula_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
        )
        .unwrap();
    let offset_values = builder.current_model_frame().unwrap().numeric_covariates[0]
        .values
        .to_vec();
    let mut normalization_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _gene in 0..counts.n_genes() {
        for offset in &offset_values {
            normalization_values.push(offset.exp());
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), normalization_values)
            .unwrap();
    let explicit_offset_builder = builder
        .clone()
        .normalization_factors(normalization_factors.clone());
    let (explicit_offset_fit, explicit_offset_results) = explicit_offset_builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
        )
        .unwrap();
    let (formula_offset_fit, formula_offset_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type` + offset(exposure)",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
        )
        .unwrap();
    assert_eq!(
        formula_offset_fit.normalization_factors,
        Some(normalization_factors.clone())
    );
    assert_lrt_fit_state_matches(
        &formula_offset_fit,
        &explicit_offset_fit,
        "top-level formula LRT offset contrast request",
    );
    assert_eq!(formula_offset_results.rows, explicit_offset_results.rows);
    let formula_replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            &options,
        )
        .unwrap();
    let explicit_replacement = builder
        .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("cell type", "B-1", "A 0", &levels),
            &options,
        )
        .unwrap();
    let explicit_offset_replacement = builder
        .clone()
        .normalization_factors(normalization_factors.clone())
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            &options,
        )
        .unwrap();
    let formula_offset_replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type` + offset(exposure)",
            "~ 1",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            &options,
        )
        .unwrap();

    assert_eq!(formula_results, explicit_results);
    let mut detached_fit = formula_fit.clone();
    detached_fit.model_frame = None;
    let explicit_object_results = detached_fit
        .lrt_results_factor_level_contrast(
            &counts,
            FactorLevelContrast::new("cell type", "B-1", "A 0", &levels),
            counts.gene_names(),
        )
        .unwrap();
    assert_eq!(explicit_object_results.rows, formula_results.rows);
    assert_eq!(explicit_object_results.metadata, formula_results.metadata);
    let missing_metadata = detached_fit
        .lrt_results_contrast_request(
            &counts,
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            None,
        )
        .unwrap_err();
    assert!(
        missing_metadata
            .to_string()
            .contains("use lrt_results_factor_level_contrast")
    );
    assert_eq!(formula_replacement.results, explicit_replacement.results);
    assert_eq!(
        formula_offset_replacement
            .original_fit
            .normalization_factors,
        Some(normalization_factors)
    );
    assert_lrt_fit_state_matches(
        &formula_offset_replacement.original_fit,
        &explicit_offset_replacement.original_fit,
        "formula LRT contrast replacement offset original",
    );
    assert_eq!(
        formula_offset_replacement.original_results,
        explicit_offset_replacement.original_results
    );
    assert_eq!(
        formula_offset_replacement.refit_results,
        explicit_offset_replacement.refit_results
    );
    assert_eq!(
        formula_offset_replacement.results,
        explicit_offset_replacement.results
    );
    assert_eq!(
        formula_results.metadata.result_name.as_deref(),
        Some("cell.type_B.1_vs_A.0")
    );
}

#[test]
fn native_formula_lrt_default_uses_cleaned_stored_reference_for_cooks_gate() {
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
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_A cell".into(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let levels = ["A cell", "A cell", "A cell", "B cell", "B cell", "B cell"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A cell".to_string(), "B cell".to_string()]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
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
        .disable_independent_filtering()
        .reduced_design(reduced.clone());

    let (_generic_fit, generic_results) = builder
        .clone()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let generic_replacement = builder
        .clone()
        .fit_lrt_with_results_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();
    let (_formula_fit, formula_results) = builder
        .clone()
        .model_frame(model_frame.clone())
        .fit_lrt_formula_with_results(&counts, "~ `cell type`", "~ 1")
        .unwrap();
    let formula_replacement = builder
        .model_frame(model_frame)
        .fit_lrt_formula_with_results_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert_eq!(formula_results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(formula_results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
    assert_eq!(
        generic_replacement.results.rows[0].cooks_outlier,
        Some(true)
    );
    assert_eq!(generic_replacement.results.rows[0].pvalue, None);
    assert_eq!(
        formula_replacement.results.rows[0].cooks_outlier,
        Some(false)
    );
    assert!(formula_replacement.results.rows[0].pvalue.is_some());
}

#[test]
fn top_level_lrt_formula_with_results_accepts_r_cleaned_quoted_coefficient_alias() {
    let counts = counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_lrt_formula_with_results_name(
            &counts,
            "~ `cell type`",
            "~ 1",
            "cell.type_B.cell_vs_T.cell",
        )
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_lrt_formula_with_results_name(
            &counts,
            "~ `cell type`",
            "~ 1",
            "cell type_B cell_vs_T cell",
        )
        .unwrap();

    assert_lrt_fit_state_matches(&cleaned_fit, &raw_fit, "quoted formula R-cleaned named LRT");
    assert_eq!(cleaned_results, raw_results);
}

#[test]
fn top_level_lrt_formula_list_contrast_accepts_r_cleaned_quoted_coefficient_alias() {
    let counts = counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
        )
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(vec!["cell type_B cell_vs_T cell".into()], Vec::new()),
        )
        .unwrap();

    assert_lrt_fit_state_matches(&cleaned_fit, &raw_fit, "quoted formula R-cleaned list LRT");
    assert_eq!(cleaned_results.rows, raw_results.rows);
    assert_eq!(
        cleaned_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_B.cell_vs_T.cell effect")
    );
    assert_eq!(
        raw_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_B cell_vs_T cell effect")
    );
}

#[test]
fn top_level_lrt_formula_list_contrast_accepts_cleaned_positive_and_negative_aliases() {
    let counts = counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "NK cell", "NK cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec![
                "T cell".to_string(),
                "B cell".to_string(),
                "NK cell".to_string(),
            ]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
        )
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
        )
        .unwrap();

    assert_lrt_fit_state_matches(
        &cleaned_fit,
        &raw_fit,
        "three-level formula R-cleaned list LRT",
    );
    assert_eq!(cleaned_results.rows, raw_results.rows);
    assert_eq!(
        cleaned_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}
