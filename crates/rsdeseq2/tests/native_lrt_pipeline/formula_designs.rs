use super::*;

#[test]
fn top_level_lrt_formula_helpers_build_full_and_reduced_designs() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
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
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let builder = native_lrt_builder().model_frame(model_frame);

    let (explicit_fit, explicit_results) = builder
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap();
    let (formula_fit, formula_results) = builder
        .fit_lrt_formula_with_results(&counts, "~ condition", "~ 1")
        .unwrap();
    let (formula_named_fit, formula_named_results) = builder
        .fit_lrt_formula_with_results_name(&counts, "~ condition", "~ 1", "condition_B_vs_A")
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
    let formula_offset_builder = builder.clone();
    let (explicit_offset_fit, explicit_offset_results) = explicit_offset_builder
        .fit_lrt_formula_with_results(&counts, "~ condition", "~ 1")
        .unwrap();
    let (formula_offset_fit, formula_offset_results) = formula_offset_builder
        .fit_lrt_formula_with_results(&counts, "~ condition + offset(exposure)", "~ 1")
        .unwrap();
    assert_eq!(
        formula_offset_fit.normalization_factors,
        Some(normalization_factors.clone())
    );
    assert_lrt_fit_state_matches(
        &formula_offset_fit,
        &explicit_offset_fit,
        "formula LRT full offset",
    );
    assert_eq!(formula_offset_results.rows, explicit_offset_results.rows);

    let offset_err = builder
        .fit_lrt_formula_with_results(&counts, "~ condition", "~ offset(exposure)")
        .unwrap_err();
    let DeseqError::UnsupportedFeature { feature } = offset_err else {
        panic!("reduced formula offsets should be rejected as unsupported, got {offset_err:?}");
    };
    assert!(feature.contains("formula offsets"));
    assert_lrt_fit_state_matches(&formula_fit, &explicit_fit, "formula LRT default");
    assert_lrt_fit_state_matches(&formula_named_fit, &explicit_fit, "formula LRT named");
    assert_eq!(formula_results.rows, explicit_results.rows);
    assert_eq!(formula_named_results.rows, explicit_results.rows);
    assert_eq!(
        formula_named_results.metadata.comparison.as_deref(),
        Some("full model versus reduced model")
    );
}

#[test]
fn fixed_lrt_supports_stored_model_frame_dot_power_design() {
    let counts = counts_with_zero_row();
    let condition = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let batch = ["X 0", "Y-1", "X 0", "Y-1", "X 0", "Y-1", "X 0", "Y-1"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: condition,
                levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch name".to_string(),
                sample_levels: batch,
                levels: Some(vec!["X 0".to_string(), "Y-1".to_string()]),
                reference: None,
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose value".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5],
        }],
    };
    let builder = native_lrt_builder().model_frame(model_frame.clone());

    let dot_full = builder.expanded_formula_design("~ .^2").unwrap();
    let explicit_full = builder
        .expanded_formula_design(
            "~ `cell type` + `batch name` + `dose value` + `cell type`:`batch name` + `cell type`:`dose value` + `batch name`:`dose value`",
        )
        .unwrap();
    let reduced = builder.expanded_formula_design("~ 1").unwrap();
    assert_eq!(dot_full, explicit_full);

    let coefficient = resolve_coefficient_index(
        &dot_full.standard_design,
        "cell.type_B.cell_vs_T.cell:batch.name_Y.1_vs_X.0",
    )
    .unwrap();
    let dispersions = vec![0.2; counts.n_genes()];
    let (dot_fit, dot_results) = builder
        .clone()
        .fit_fixed_dispersion_lrt(
            &counts,
            &dot_full.standard_design,
            &reduced.standard_design,
            &dispersions,
            coefficient,
        )
        .unwrap();
    let (explicit_fit, explicit_results) = builder
        .fit_fixed_dispersion_lrt(
            &counts,
            &explicit_full.standard_design,
            &reduced.standard_design,
            &dispersions,
            coefficient,
        )
        .unwrap();

    assert_lrt_fit_state_matches(
        &dot_fit,
        &explicit_fit,
        "stored model-frame dot-power fixed LRT",
    );
    assert_eq!(dot_results, explicit_results);
    assert_eq!(dot_fit.current_model_frame(), Some(&model_frame));
}

#[test]
fn top_level_lrt_formula_contrast_helpers_use_model_frame_metadata() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
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
    let builder = native_lrt_builder().model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let (_factor_fit, factor_results) = builder
        .fit_lrt_with_results_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (_formula_character_fit, formula_character_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ condition",
            "~ 1",
            &ResultsContrast::character("condition", "B", "A"),
        )
        .unwrap();
    let (_formula_list_fit, formula_list_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ condition",
            "~ 1",
            &ResultsContrast::list(vec!["condition_B_vs_A".into()], Vec::new()),
        )
        .unwrap();
    let (_formula_numeric_fit, formula_numeric_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ condition",
            "~ 1",
            &ResultsContrast::numeric(vec![0.0, 1.0]),
        )
        .unwrap();
    let formula_replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ condition",
            "~ 1",
            &ResultsContrast::character("condition", "B", "A"),
            &options,
        )
        .unwrap();
    let explicit_replacement = builder
        .fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &options,
        )
        .unwrap();

    assert_eq!(formula_character_results, factor_results);
    assert_eq!(formula_list_results.rows, factor_results.rows);
    assert_eq!(formula_numeric_results.rows, factor_results.rows);
    assert_eq!(formula_replacement.results, explicit_replacement.results);
    assert_eq!(
        formula_character_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    assert_eq!(
        formula_list_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: condition_B_vs_A effect")
    );
    assert_eq!(
        formula_numeric_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
}

#[test]
fn native_formula_lrt_default_uses_stored_formula_low_count_cooks_gate() {
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
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
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
    let (_formula_fit, formula_results) = builder
        .model_frame(model_frame)
        .fit_lrt_formula_with_results(&counts, "~ condition", "~ 1")
        .unwrap();

    assert_eq!(formula_results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(formula_results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
}

#[test]
fn top_level_lrt_formula_fit_state_retains_model_frame_metadata() {
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
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose value".to_string(),
            values: vec![0.5, 1.5, 0.25, 1.25, 0.75, 1.75, 0.4, 1.4],
        }],
    };
    let builder = native_lrt_builder().model_frame(model_frame.clone());
    let options = CooksReplacementOptions::new(f64::MAX);
    let transformed_formula = "~ relevel(`cell type`, ref='B cell')";

    let (fit, _results) = builder
        .clone()
        .fit_lrt_formula_with_results(&counts, transformed_formula, "~ 1")
        .unwrap();
    let (named_fit, _named_results) = builder
        .clone()
        .fit_lrt_formula_with_results_name(
            &counts,
            transformed_formula,
            "~ 1",
            "relevel.cell.type..ref....B.cell.._T.cell_vs_B.cell",
        )
        .unwrap();
    let fit_only = builder
        .clone()
        .fit_lrt_formula(&counts, transformed_formula, "~ 1")
        .unwrap();
    let named_fit_only = builder
        .clone()
        .fit_lrt_formula_name(
            &counts,
            transformed_formula,
            "~ 1",
            "relevel.cell.type..ref....B.cell.._T.cell_vs_B.cell",
        )
        .unwrap();
    let replacement = builder
        .clone()
        .fit_lrt_formula_with_results_with_cooks_replacement(
            &counts,
            transformed_formula,
            "~ 1",
            &options,
        )
        .unwrap();
    let named_replacement = builder
        .fit_lrt_formula_with_results_name_with_cooks_replacement(
            &counts,
            transformed_formula,
            "~ 1",
            "relevel.cell.type..ref....B.cell.._T.cell_vs_B.cell",
            &options,
        )
        .unwrap();

    let transformed_model_frame = fit.current_model_frame().unwrap();
    assert_eq!(
        transformed_model_frame
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["cell type", "relevel(cell type, ref = \"B cell\")"]
    );
    assert_eq!(
        named_fit.current_model_frame(),
        Some(transformed_model_frame)
    );
    assert_eq!(
        fit_only.current_model_frame(),
        Some(transformed_model_frame)
    );
    assert_eq!(
        named_fit_only.current_model_frame(),
        Some(transformed_model_frame)
    );
    assert_eq!(
        replacement.original_fit.current_model_frame(),
        Some(transformed_model_frame)
    );
    assert_eq!(
        named_replacement.original_fit.current_model_frame(),
        Some(transformed_model_frame)
    );
    assert_eq!(
        replacement
            .refit_fit
            .as_ref()
            .map(DeseqFit::current_model_frame),
        named_replacement
            .refit_fit
            .as_ref()
            .map(DeseqFit::current_model_frame)
    );
    assert_eq!(
        replacement
            .refit_fit
            .as_ref()
            .map(DeseqFit::current_model_frame),
        None
    );

    let numeric_formula = "~ `cell type` + as.numeric(`dose value`)";
    let numeric_fit = builder
        .fit_lrt_formula(&counts, numeric_formula, "~ `cell type`")
        .expect("LRT formula with derived numeric metadata should fit");
    let numeric_frame = numeric_fit.current_model_frame().unwrap();
    assert_eq!(
        numeric_frame
            .numeric_covariates
            .iter()
            .map(|covariate| covariate.name.as_str())
            .collect::<Vec<_>>(),
        vec!["dose value", "dose value_as.numeric"]
    );
    assert_eq!(
        numeric_frame
            .resolved_numeric_covariate_by_alias("dose.value_as.numeric")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "dose value_as.numeric",
            values: &[0.5, 1.5, 0.25, 1.25, 0.75, 1.75, 0.4, 1.4],
        }
    );

    let factor_code_formula = "~ as.double(`cell type`)";
    let factor_code_fit = builder
        .clone()
        .fit_lrt_formula(&counts, factor_code_formula, "~ 1")
        .expect("LRT formula with factor-code numeric metadata should fit");
    let factor_code_replacement = builder
        .fit_lrt_formula_with_results_name_with_cooks_replacement(
            &counts,
            factor_code_formula,
            "~ 1",
            "cell.type_as.double",
            &options,
        )
        .expect("LRT factor-code formula replacement should resolve cleaned coefficient alias");
    let factor_code_frame = factor_code_fit.current_model_frame().unwrap();
    assert_eq!(
        factor_code_frame
            .numeric_covariates
            .iter()
            .map(|covariate| covariate.name.as_str())
            .collect::<Vec<_>>(),
        vec!["dose value", "cell type_as.double"]
    );
    assert_eq!(
        factor_code_frame
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "cell type_as.double",
            values: &[1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
        }
    );
    assert_eq!(factor_code_replacement.refit_plan.n_refit, 0);
    assert_eq!(
        factor_code_replacement.original_fit.current_model_frame(),
        factor_code_fit.current_model_frame()
    );
    assert_eq!(
        factor_code_replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("cell type_as.double")
    );
}
