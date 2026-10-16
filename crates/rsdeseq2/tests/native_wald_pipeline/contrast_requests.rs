use super::*;

#[test]
fn top_level_results_contrast_request_routes_character_list_and_numeric_forms() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().disable_independent_filtering();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (factor_fit, factor_results) = builder
        .clone()
        .fit_with_results_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (character_fit, character_results) = builder
        .clone()
        .fit_with_results_contrast_request(
            &counts,
            &design,
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
        )
        .unwrap();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let (model_frame_fit, model_frame_results) = builder
        .clone()
        .fit_with_results_contrast_request_from_model_frame(
            &counts,
            &design,
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
        )
        .unwrap();
    let stored_builder = builder.clone().model_frame(model_frame.clone());
    let stored_design = stored_builder
        .expanded_formula_design("~ condition")
        .unwrap();
    assert_eq!(stored_design.standard_design, design);
    let stored_design_with_offsets = stored_builder
        .expanded_formula_design_with_offsets("~ condition")
        .unwrap();
    assert_eq!(stored_design_with_offsets.design, stored_design);
    assert!(stored_design_with_offsets.offsets.is_empty());
    assert!(builder.expanded_formula_design("~ condition").is_err());

    let (stored_model_frame_fit, stored_model_frame_results) = stored_builder
        .clone()
        .fit_with_results_contrast_request(
            &counts,
            &stored_design.standard_design,
            &ResultsContrast::character("condition", "B", "A"),
            None::<&[String]>,
        )
        .unwrap();
    let (formula_fit, formula_results) = stored_builder
        .clone()
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ condition",
            &ResultsContrast::character("condition", "B", "A"),
        )
        .unwrap();
    let (formula_name_fit, formula_name_results) = stored_builder
        .clone()
        .fit_formula_with_results_name(&counts, "~ condition", "condition_B_vs_A")
        .unwrap();
    let formula_name_fit_only = stored_builder
        .clone()
        .fit_formula_name(&counts, "~ condition", "condition_B_vs_A")
        .unwrap();
    let (formula_default_fit, formula_default_results) = stored_builder
        .clone()
        .fit_formula_with_results(&counts, "~ condition")
        .unwrap();
    let formula_default_fit_only = stored_builder
        .clone()
        .fit_formula(&counts, "~ condition")
        .unwrap();
    let offset_values = model_frame.numeric_covariates[0].values.to_vec();
    let mut normalization_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _gene in 0..counts.n_genes() {
        for offset in &offset_values {
            normalization_values.push(offset.exp());
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), normalization_values)
            .unwrap();
    let explicit_offset_builder = stored_builder
        .clone()
        .size_factors(vec![1.0; counts.n_samples()])
        .normalization_factors(normalization_factors.clone());
    let formula_offset_builder = stored_builder
        .clone()
        .size_factors(vec![1.0; counts.n_samples()]);
    let (explicit_offset_fit, explicit_offset_results) = explicit_offset_builder
        .fit_formula_with_results(&counts, "~ condition")
        .unwrap();
    let (formula_offset_fit, formula_offset_results) = formula_offset_builder
        .fit_formula_with_results(&counts, "~ condition + offset(exposure)")
        .unwrap();
    assert_eq!(
        formula_offset_fit.normalization_factors,
        Some(normalization_factors)
    );
    assert_wald_fit_state_matches(
        &formula_offset_fit,
        &explicit_offset_fit,
        "top-level stored model-frame formula offset request",
    );
    assert_result_rows_close(
        &formula_offset_results.rows,
        &explicit_offset_results.rows,
        "top-level stored model-frame formula offset results",
    );
    let (explicit_offset_contrast_fit, explicit_offset_contrast_results) = explicit_offset_builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ condition",
            &ResultsContrast::character("condition", "B", "A"),
        )
        .unwrap();
    let (formula_offset_contrast_fit, formula_offset_contrast_results) = formula_offset_builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ condition + offset(exposure)",
            &ResultsContrast::character("condition", "B", "A"),
        )
        .unwrap();
    assert_wald_fit_state_matches(
        &formula_offset_contrast_fit,
        &explicit_offset_contrast_fit,
        "top-level stored model-frame formula offset contrast request",
    );
    assert_result_rows_close(
        &formula_offset_contrast_results.rows,
        &explicit_offset_contrast_results.rows,
        "top-level stored model-frame formula offset contrast results",
    );
    assert_wald_fit_state_matches(
        &character_fit,
        &factor_fit,
        "top-level results character contrast request",
    );
    assert_eq!(character_results, factor_results);
    assert_wald_fit_state_matches(
        &model_frame_fit,
        &factor_fit,
        "top-level model-frame results character contrast request",
    );
    assert_eq!(model_frame_fit.current_model_frame(), Some(&model_frame));
    assert_eq!(model_frame_results, factor_results);
    assert_wald_fit_state_matches(
        &stored_model_frame_fit,
        &factor_fit,
        "top-level stored model-frame results character contrast request",
    );
    assert_eq!(stored_model_frame_results, factor_results);
    assert_wald_fit_state_matches(
        &formula_fit,
        &factor_fit,
        "top-level stored model-frame formula contrast request",
    );
    assert_eq!(formula_results, factor_results);
    assert_wald_fit_state_matches(
        &formula_name_fit,
        &factor_fit,
        "top-level stored model-frame formula named request",
    );
    assert_wald_fit_state_matches(
        &formula_name_fit_only,
        &formula_name_fit,
        "top-level stored model-frame formula named fit-only request",
    );
    assert_result_rows_close(
        &formula_name_results.rows,
        &factor_results.rows,
        "top-level stored model-frame formula named results",
    );
    assert_wald_fit_state_matches(
        &formula_default_fit,
        &factor_fit,
        "top-level stored model-frame formula default request",
    );
    assert_wald_fit_state_matches(
        &formula_default_fit_only,
        &formula_default_fit,
        "top-level stored model-frame formula default fit-only request",
    );
    assert_result_rows_close(
        &formula_default_results.rows,
        &factor_results.rows,
        "top-level stored model-frame formula default results",
    );

    let (_numeric_fit, numeric_results) = builder
        .clone()
        .fit_with_results_contrast_request(
            &counts,
            &design,
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            None::<&[String]>,
        )
        .unwrap();
    assert_eq!(
        numeric_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        numeric_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    assert_eq!(numeric_results.metadata.contrast, Some(vec![0.0, 1.0]));

    let (_list_fit, list_results) = builder
        .fit_with_results_contrast_request(
            &counts,
            &design,
            &ResultsContrast::list(vec!["condition_B_vs_A".into()], Vec::new()),
            None::<&[String]>,
        )
        .unwrap();
    assert_eq!(
        list_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        list_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: condition_B_vs_A effect")
    );
    assert_eq!(list_results.metadata.contrast, Some(vec![0.0, 1.0]));
}

#[test]
fn top_level_results_contrast_request_canonicalizes_cleaned_sample_level_aliases() {
    let counts = native_wald_counts_with_zero_row();
    let design = DesignMatrix::from_row_major(
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
    let builder = glm_mu_native_wald_builder().disable_independent_filtering();
    let levels = [
        "A cell", "A cell", "A cell", "A cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();

    let (_factor_fit, factor_results) = builder
        .clone()
        .fit_with_results_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("cell type", "B cell", "A cell", &levels),
        )
        .unwrap();
    let (_character_fit, character_results) = builder
        .clone()
        .fit_with_results_contrast_request(
            &counts,
            &design,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&levels),
        )
        .unwrap();
    let replacement_options = CooksReplacementOptions::new(f64::MAX);
    let factor_replacement = builder
        .clone()
        .fit_with_test_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &design,
            FactorLevelContrast::new("cell type", "B cell", "A cell", &levels),
            &replacement_options,
        )
        .unwrap();
    let request_replacement = builder
        .fit_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            &design,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&levels),
            &replacement_options,
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(factor_replacement) = factor_replacement else {
        panic!("factor contrast should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(request_replacement) = request_replacement else {
        panic!("character contrast request should route to Wald");
    };

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
    let ambiguous = builder.fit_with_results_contrast_request(
        &counts,
        &design,
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
fn top_level_formula_contrast_request_uses_cleaned_stored_model_frame_factor_name() {
    let counts = native_wald_counts_with_zero_row();
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
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder().model_frame(model_frame);
    let design = builder
        .expanded_formula_design("~ `cell type`")
        .unwrap()
        .standard_design;

    let (explicit_fit, explicit_results) = builder
        .fit_with_results_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("cell type", "B-1", "A 0", &levels),
        )
        .unwrap();
    let (formula_fit, formula_results) = builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &formula_fit,
        &explicit_fit,
        "top-level formula cleaned model-frame factor contrast",
    );
    assert_eq!(formula_results, explicit_results);
    assert_eq!(
        formula_results.metadata.result_name.as_deref(),
        Some("cell.type_B.1_vs_A.0")
    );
}

#[test]
fn top_level_formula_contrast_request_uses_formula_derived_factor_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let builder = glm_mu_native_wald_builder().model_frame(model_frame);
    let formula = "~ factor(condition, levels=c('B','A'))";
    let formula_design = builder
        .expanded_formula_design_with_offsets(formula)
        .unwrap();
    let design = &formula_design.design.standard_design;

    let transformed_factor = "factor(condition, levels = c(\"B\", \"A\"))";
    let (explicit_fit, explicit_results) = builder
        .fit_with_results_factor_level_contrast(
            &counts,
            design,
            FactorLevelContrast::with_reference(transformed_factor, "A", "B", "B", &levels),
        )
        .unwrap();
    let (formula_fit, formula_results) = builder
        .fit_formula_with_results_contrast_request(
            &counts,
            formula,
            &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &formula_fit,
        &explicit_fit,
        "top-level formula-derived factor contrast",
    );
    assert_eq!(formula_results, explicit_results);
    let mut detached_fit = formula_fit.clone();
    detached_fit.model_frame = None;
    let explicit_object_results = detached_fit
        .wald_results_factor_level_contrast(
            &counts,
            FactorLevelContrast::with_reference(transformed_factor, "A", "B", "B", &levels),
            counts.gene_names(),
        )
        .unwrap();
    assert_eq!(explicit_object_results.rows, formula_results.rows);
    assert_eq!(explicit_object_results.metadata, formula_results.metadata);
    let missing_metadata = detached_fit
        .wald_results_contrast_request(
            &counts,
            &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
            None,
        )
        .unwrap_err();
    assert!(
        missing_metadata
            .to_string()
            .contains("use wald_results_factor_level_contrast")
    );
    assert_eq!(
        formula_results.metadata.result_name.as_deref(),
        Some("factor.condition..levels...c..B....A..._A_vs_B")
    );
}

#[test]
fn top_level_formula_contrast_replacement_uses_formula_derived_factor_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let builder = glm_mu_native_wald_builder().model_frame(model_frame);
    let formula = "~ factor(condition, levels=c('B','A'))";
    let options = CooksReplacementOptions::new(f64::MAX);

    let output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            formula,
            &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(output) = output else {
        panic!("formula-derived contrast replacement should route to Wald");
    };

    assert_eq!(output.refit_plan.n_refit, 0);
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("factor.condition..levels...c..B....A..._A_vs_B")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: factor(condition, levels = c(\"B\", \"A\")) A vs B")
    );

    let ordered_formula = "~ ordered(condition, levels=c('B','A'))";
    let ordered_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            ordered_formula,
            &ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(ordered_output) = ordered_output else {
        panic!("ordered formula-derived contrast replacement should route to Wald");
    };

    assert_eq!(ordered_output.refit_plan.n_refit, 0);
    assert_eq!(ordered_output.results.rows, output.results.rows);
    assert_eq!(
        ordered_output.results.metadata.result_name.as_deref(),
        Some("ordered.condition..levels...c..B....A..._A_vs_B")
    );
    assert_eq!(
        ordered_output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: ordered(condition, levels = c(\"B\", \"A\")) A vs B")
    );

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
    let explicit_offset_output = builder
        .clone()
        .normalization_factors(normalization_factors.clone())
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            formula,
            &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let formula_offset_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ factor(condition, levels=c('B','A')) + offset(exposure)",
            &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(explicit_offset_output) = explicit_offset_output else {
        panic!("explicit-offset formula-derived contrast replacement should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(formula_offset_output) = formula_offset_output else {
        panic!("formula-offset formula-derived contrast replacement should route to Wald");
    };
    assert_eq!(
        formula_offset_output.original_fit.normalization_factors,
        Some(normalization_factors.clone())
    );
    assert_wald_fit_state_matches(
        &formula_offset_output.original_fit,
        &explicit_offset_output.original_fit,
        "formula-derived contrast replacement offset original",
    );
    assert_eq!(
        formula_offset_output.original_results,
        explicit_offset_output.original_results
    );
    assert_eq!(
        formula_offset_output.refit_results,
        explicit_offset_output.refit_results
    );
    assert_eq!(
        formula_offset_output.results,
        explicit_offset_output.results
    );

    let explicit_ordered_offset_output = builder
        .clone()
        .normalization_factors(normalization_factors.clone())
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            ordered_formula,
            &ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let formula_ordered_offset_output = builder
        .fit_formula_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ ordered(condition, levels=c('B','A')) + offset(exposure)",
            &ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B"),
            &options,
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(explicit_ordered_offset_output) =
        explicit_ordered_offset_output
    else {
        panic!("explicit-offset ordered formula-derived contrast replacement should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(formula_ordered_offset_output) =
        formula_ordered_offset_output
    else {
        panic!("formula-offset ordered formula-derived contrast replacement should route to Wald");
    };
    assert_wald_fit_state_matches(
        &formula_ordered_offset_output.original_fit,
        &explicit_ordered_offset_output.original_fit,
        "ordered formula-derived contrast replacement offset original",
    );
    assert_eq!(
        formula_ordered_offset_output.results,
        explicit_ordered_offset_output.results
    );
}

#[test]
fn top_level_results_character_contrast_request_requires_sample_levels() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let err = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .fit_with_results_contrast_request(
            &counts,
            &design,
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
fn top_level_formula_model_frame_cleans_factor_result_names() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            0, 0, 10, 10, //
            10, 20, 40, 50,
        ],
    )
    .unwrap();
    let levels = ["A 0", "A 0", "B-1", "B-1"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .size_factors(vec![1.0; 4])
        .disable_independent_filtering()
        .disable_cooks_cutoff()
        .model_frame(model_frame);
    let formula_design = builder.expanded_formula_design("~ condition").unwrap();
    assert_eq!(
        formula_design.standard_design.coefficient_names().unwrap(),
        &["Intercept".to_string(), "condition_B-1_vs_A 0".to_string()]
    );

    let (_fit, results) = builder
        .fit_fixed_dispersion_wald_results_contrast::<String>(
            &counts,
            &formula_design.standard_design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B-1", "A 0"),
            None::<&[String]>,
        )
        .unwrap();

    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B.1_vs_A.0")
    );
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B-1 vs A 0")
    );
    assert!(results.rows[0].log2_fold_change.is_some());
    assert!(results.rows[0].pvalue.is_some());
}
