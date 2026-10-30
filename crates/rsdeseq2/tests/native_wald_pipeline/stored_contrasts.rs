use super::*;

#[test]
fn fitted_wald_object_uses_stored_model_frame_for_character_contrast_results() {
    let counts = native_wald_counts_with_zero_row();
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
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);

    let (fit, mut expected_results) = builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
        )
        .unwrap();
    expected_results.independent_filtering = None;
    let object_results = fit
        .wald_results_contrast_request(
            &counts,
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
            counts.gene_names(),
        )
        .unwrap();
    let dispatched_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
            counts.gene_names(),
        )
        .unwrap();
    let list_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
            counts.gene_names(),
        )
        .unwrap();
    let numeric_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            counts.gene_names(),
        )
        .unwrap();

    assert_eq!(object_results, expected_results);
    assert_eq!(dispatched_results, expected_results);
    assert_eq!(list_results.rows, expected_results.rows);
    assert_eq!(numeric_results.rows, expected_results.rows);

    let transformed_levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let transformed_model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: transformed_levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let transformed_builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(transformed_model_frame);
    let transformed_formula = "~ factor(condition, levels=c('B','A'))";
    let transformed_contrast =
        ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B");
    let (transformed_fit, mut expected_transformed_results) = transformed_builder
        .fit_formula_with_results_contrast_request(
            &counts,
            transformed_formula,
            &transformed_contrast,
        )
        .unwrap();
    expected_transformed_results.independent_filtering = None;
    let object_transformed_results = transformed_fit
        .results_contrast_request(&counts, &transformed_contrast, counts.gene_names())
        .unwrap();

    assert_eq!(object_transformed_results, expected_transformed_results);
    assert_eq!(
        object_transformed_results.metadata.result_name.as_deref(),
        Some("factor.condition..levels...c..B....A..._A_vs_B")
    );

    let ordered_formula = "~ ordered(condition, levels=c('B','A'))";
    let ordered_contrast = ResultsContrast::numeric(vec![0.0, 1.0]);
    let (ordered_fit, mut expected_ordered_results) = transformed_builder
        .fit_formula_with_results_contrast_request(&counts, ordered_formula, &ordered_contrast)
        .unwrap();
    expected_ordered_results.independent_filtering = None;
    let object_ordered_results = ordered_fit
        .results_contrast_request(&counts, &ordered_contrast, counts.gene_names())
        .unwrap();
    assert_eq!(object_ordered_results, expected_ordered_results);
    let invalid = ordered_fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B"),
            counts.gene_names(),
        )
        .unwrap_err();
    assert!(invalid.to_string().contains("ordered factor"));

    assert_eq!(
        object_ordered_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        object_ordered_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );

    let relabeled_builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(FormulaModelFrame {
            factors: vec![FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec![
                    "A".to_string(),
                    "A".to_string(),
                    "A".to_string(),
                    "A".to_string(),
                    "B".to_string(),
                    "B".to_string(),
                    "B".to_string(),
                    "B".to_string(),
                ],
                levels: Some(vec!["A".to_string(), "B".to_string()]),
                reference: Some("A".to_string()),
            }],
            numeric_covariates: Vec::new(),
        });
    let relabeled_formula = "~ factor(condition, levels=c('A','B'), labels=c('control','treated'))";
    let relabeled_contrast = ResultsContrast::character(
        "factor.condition..levels...c..A....B....labels...c..control....treated...",
        "treated",
        "control",
    );
    let (relabeled_fit, mut expected_relabeled_results) = relabeled_builder
        .fit_formula_with_results_contrast_request(&counts, relabeled_formula, &relabeled_contrast)
        .unwrap();
    expected_relabeled_results.independent_filtering = None;
    let object_relabeled_results = relabeled_fit
        .results_contrast_request(&counts, &relabeled_contrast, counts.gene_names())
        .unwrap();
    assert_eq!(object_relabeled_results, expected_relabeled_results);
    assert_eq!(
        object_relabeled_results.metadata.result_name.as_deref(),
        Some(
            "factor.condition..levels...c..A....B....labels...c..control....treated..._treated_vs_control"
        )
    );
    assert_eq!(
        object_relabeled_results.metadata.comparison.as_deref(),
        Some(
            "factor-level contrast: factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\")) treated vs control"
        )
    );

    let error = relabeled_builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ ordered(condition, levels=c('A','B'), labels=c('control','treated'))",
            &ResultsContrast::character(
                "ordered.condition..levels...c..A....B....labels...c..control....treated...",
                "treated",
                "control",
            ),
        )
        .unwrap_err();
    assert!(error.to_string().contains("ordered factor"));
}

#[test]
fn fitted_wald_object_applies_numeric_contrast_all_zero_from_stored_design() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            0, 0, 0, 0, 45, 55, //
            12, 11, 24, 22, 18, 19, //
            0, 0, 0, 0, 0, 0,
        ],
        Some(vec![
            "outside_contrast".into(),
            "tested".into(),
            "all_zero".into(),
        ]),
        None,
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "C", "C"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".into(), "B".into(), "C".into()]),
            reference: Some("A".into()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let design = builder
        .expanded_formula_design("~ 0 + condition")
        .unwrap()
        .standard_design;
    let dispersions = vec![0.2; counts.n_genes()];
    let (fit, character_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &dispersions,
            &ResultsContrast::character("condition", "B", "A"),
            builder.current_model_frame().unwrap(),
        )
        .unwrap();
    let list_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::list(vec!["conditionB".into()], vec!["conditionA".into()]),
            counts.gene_names(),
        )
        .unwrap();
    let numeric_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::numeric(vec![-1.0, 1.0, 0.0]),
            counts.gene_names(),
        )
        .unwrap();

    assert_eq!(list_results.rows, numeric_results.rows);
    assert_eq!(character_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(list_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(list_results.rows[0].stat, Some(0.0));
    assert_eq!(list_results.rows[0].pvalue, Some(1.0));
    assert_eq!(numeric_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(numeric_results.rows[0].stat, Some(0.0));
    assert_eq!(numeric_results.rows[0].pvalue, Some(1.0));
    assert_eq!(numeric_results.rows[2].log2_fold_change, None);

    let formula_design = builder
        .expanded_formula_design_with_offsets("~ factor(condition, levels=c('A','B','C'))")
        .unwrap();
    let transformed_contrast =
        ResultsContrast::character("factor.condition..levels...c..A....B....C...", "B", "A");
    let (formula_fit, mut formula_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &formula_design.design.standard_design,
            &dispersions,
            &transformed_contrast,
            &formula_design.model_frame,
        )
        .unwrap();
    formula_results.independent_filtering = None;
    let formula_object_results = formula_fit
        .results_contrast_request(&counts, &transformed_contrast, counts.gene_names())
        .unwrap();

    assert_eq!(formula_object_results, formula_results);
    assert_eq!(
        formula_fit
            .current_model_frame()
            .unwrap()
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec![
            "condition",
            "factor(condition, levels = c(\"A\", \"B\", \"C\"))",
        ]
    );
}

#[test]
fn top_level_formula_list_contrast_cooks_replacement_accepts_r_cleaned_quoted_coefficient_alias() {
    let counts = native_wald_counts_with_zero_row();
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
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(vec!["cell type_B cell_vs_T cell".into()], Vec::new()),
            &options,
        )
        .unwrap();

    let CooksReplacementTestOutput::Wald(cleaned_output) = cleaned_output else {
        panic!("formula list contrast replacement should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(raw_output) = raw_output else {
        panic!("formula list contrast replacement should route to Wald");
    };
    assert_wald_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "quoted formula R-cleaned list replacement Wald original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_B.cell_vs_T.cell effect")
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_B cell_vs_T cell effect")
    );
}

#[test]
fn top_level_formula_list_contrast_cooks_replacement_accepts_cleaned_numeric_transform_alias() {
    let counts = native_wald_counts_with_zero_row();
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
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            &ResultsContrast::list(vec!["cell.type_as.double".into()], Vec::new()),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            &ResultsContrast::list(vec!["cell type_as.double".into()], Vec::new()),
            &options,
        )
        .unwrap();

    let CooksReplacementTestOutput::Wald(cleaned_output) = cleaned_output else {
        panic!("formula numeric-transform list replacement should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(raw_output) = raw_output else {
        panic!("formula numeric-transform list replacement should route to Wald");
    };
    assert_wald_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "numeric-transform formula R-cleaned list replacement Wald original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output
            .original_fit
            .current_model_frame()
            .unwrap()
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap()
            .name,
        "cell type_as.double"
    );
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_as.double effect")
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_as.double effect")
    );
}

#[test]
fn top_level_formula_numeric_contrast_cooks_replacement_retains_numeric_transform_metadata() {
    let counts = native_wald_counts_with_zero_row();
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
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let numeric_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            &options,
        )
        .unwrap();

    let CooksReplacementTestOutput::Wald(numeric_output) = numeric_output else {
        panic!("formula numeric-transform numeric replacement should route to Wald");
    };
    assert_eq!(numeric_output.refit_plan.n_refit, 0);
    assert_eq!(
        numeric_output
            .original_fit
            .current_model_frame()
            .unwrap()
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap()
            .name,
        "cell type_as.double"
    );
    assert_eq!(
        numeric_output.results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
}

#[test]
fn top_level_formula_list_contrast_cooks_replacement_accepts_cleaned_positive_and_negative_aliases()
{
    let counts = native_wald_counts_with_zero_row();
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
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
            &options,
        )
        .unwrap();

    let CooksReplacementTestOutput::Wald(cleaned_output) = cleaned_output else {
        panic!("formula list contrast replacement should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(raw_output) = raw_output else {
        panic!("formula list contrast replacement should route to Wald");
    };
    assert_wald_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "three-level formula R-cleaned list replacement Wald original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}
