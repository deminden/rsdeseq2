use super::*;

#[test]
fn builder_try_model_frame_validates_wrapper_metadata() {
    let valid = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: vec![0.0, 1.0],
        }],
    };
    let builder = DeseqBuilder::new().try_model_frame(valid.clone()).unwrap();
    assert_eq!(builder.current_model_frame(), Some(&valid));

    let ambiguous = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "cell-type".to_string(),
            values: vec![0.0, 1.0],
        }],
    };
    assert_eq!(
        DeseqBuilder::new()
            .model_frame(ambiguous.clone())
            .current_model_frame(),
        Some(&ambiguous)
    );
    let error = DeseqBuilder::new()
        .try_model_frame(ambiguous)
        .unwrap_err()
        .to_string();
    assert!(error.contains("resolves ambiguously after R-style cleanup"));
}

#[test]
fn model_frame_contrast_routes_validate_metadata_for_numeric_requests() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![10, 20]).unwrap();
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
    let ambiguous = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "cell-type".to_string(),
            values: vec![0.0, 1.0],
        }],
    };

    let error = DeseqBuilder::new()
        .size_factors(vec![1.0; 2])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1],
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            &ambiguous,
        )
        .unwrap_err()
        .to_string();

    assert!(error.contains("resolves ambiguously after R-style cleanup"));
}

#[test]
fn fixed_dispersion_wald_results_character_contrast_uses_formula_model_frame() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        6,
        vec![
            0, 0, 0, 0, 10, 10, //
            10, 10, 20, 20, 40, 40,
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
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_explicit_fit, explicit_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
        )
        .unwrap();
    let (model_frame_fit, model_frame_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
        )
        .unwrap();
    let builder_with_model_frame = builder.clone().model_frame(model_frame.clone());
    assert_eq!(
        builder_with_model_frame.current_model_frame(),
        Some(&model_frame)
    );
    let (stored_model_frame_fit, stored_model_frame_results) = builder_with_model_frame
        .fit_fixed_dispersion_wald_results_contrast::<String>(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            None,
        )
        .unwrap();

    assert_eq!(model_frame_results, explicit_results);
    assert_eq!(stored_model_frame_results, explicit_results);
    assert_eq!(model_frame_fit.current_model_frame(), Some(&model_frame));
    assert_eq!(
        stored_model_frame_fit.current_model_frame(),
        Some(&model_frame)
    );
    assert_eq!(model_frame_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(model_frame_results.rows[0].pvalue, Some(1.0));
    assert_eq!(
        model_frame_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn fixed_dispersion_wald_model_frame_contrast_uses_declared_factor_levels() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        6,
        vec![
            8, 5, 9, 6, 14, 16, //
            30, 20, 36, 24, 50, 58,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, //
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
    let observed_levels = ["B", "A", "B", "A", "C", "C"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: observed_levels.clone(),
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_explicit_fit, explicit_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::with_reference("condition", "C", "B", "A", &observed_levels),
        )
        .unwrap();
    let (_model_frame_fit, model_frame_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "C", "B"),
            &model_frame,
        )
        .unwrap();

    assert_eq!(model_frame_results, explicit_results);
    assert_eq!(
        model_frame_results.metadata.result_name.as_deref(),
        Some("condition_C_vs_B")
    );
    assert_eq!(
        model_frame_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition C vs B")
    );
}

#[test]
fn fixed_dispersion_wald_model_frame_contrast_uses_unused_declared_reference_with_expanded_design()
{
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            8, 9, 14, 16, //
            30, 36, 50, 58,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            0.0, 1.0, //
            0.0, 1.0,
        ],
        Some(vec!["conditionB".into(), "conditionC".into()]),
    )
    .unwrap();
    let observed_levels = ["B", "B", "C", "C"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: observed_levels.clone(),
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 4])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_explicit_fit, explicit_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::with_reference("condition", "C", "B", "A", &observed_levels),
        )
        .unwrap();
    let (_model_frame_fit, model_frame_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "C", "B"),
            &model_frame,
        )
        .unwrap();

    assert_eq!(model_frame_results, explicit_results);
    assert_eq!(
        model_frame_results.metadata.result_name.as_deref(),
        Some("condition_C_vs_B")
    );
    assert_eq!(
        model_frame_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition C vs B")
    );
}

#[test]
fn fixed_dispersion_wald_model_frame_contrast_accepts_cleaned_factor_name_alias() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            0, 0, 10, 10, //
            10, 20, 40, 50,
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
        Some(vec!["Intercept".into(), "cell.type_B.1_vs_A.0".into()]),
    )
    .unwrap();
    let levels = ["A 0", "A 0", "B-1", "B-1"]
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
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 4])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_explicit_fit, explicit_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::new("cell type", "B-1", "A 0", &levels),
        )
        .unwrap();
    let (_model_frame_fit, model_frame_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            &model_frame,
        )
        .unwrap();
    let (_stored_fit, stored_results) = builder
        .clone()
        .model_frame(model_frame)
        .fit_fixed_dispersion_wald_results_contrast::<String>(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("cell.type", "B-1", "A 0"),
            None,
        )
        .unwrap();
    let (_alias_fit, alias_results) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("cell.type", "B.1", "A.0"),
            &FormulaModelFrame {
                factors: vec![FormulaFactorColumn {
                    name: "cell type".to_string(),
                    sample_levels: levels.clone(),
                    levels: Some(vec!["A 0".to_string(), "B-1".to_string()]),
                    reference: None,
                }],
                numeric_covariates: Vec::new(),
            },
        )
        .unwrap();

    assert_eq!(model_frame_results, explicit_results);
    assert_eq!(stored_results, explicit_results);
    assert_eq!(alias_results, explicit_results);
    assert_eq!(
        model_frame_results.metadata.result_name.as_deref(),
        Some("cell.type_B.1_vs_A.0")
    );
    assert_eq!(
        model_frame_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: cell type B-1 vs A 0")
    );
}

#[test]
fn fixed_dispersion_wald_formula_ordered_contrast_requires_polynomial_selectors() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            8, 9, 14, 16, //
            30, 36, 50, 58,
        ],
    )
    .unwrap();
    let levels = ["A", "A", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 4])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .model_frame(FormulaModelFrame {
            factors: vec![FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: levels,
                levels: Some(vec!["A".to_string(), "B".to_string()]),
                reference: Some("A".to_string()),
            }],
            numeric_covariates: Vec::new(),
        });
    let formula_design = builder
        .expanded_formula_design_with_offsets("~ ordered(condition, levels=c('B','A'))")
        .unwrap();
    let design = &formula_design.design.standard_design;
    let model_frame = &formula_design.model_frame;
    let dispersions = [0.1, 0.1];
    let invalid = ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B");
    let error = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            design,
            &dispersions,
            &invalid,
            model_frame,
        )
        .unwrap_err();
    assert!(error.to_string().contains("ordered factor"));
    let numeric = ResultsContrast::numeric(vec![0.0, 1.0]);
    let named = ResultsContrast::list(vec![design.coefficient_names().unwrap()[1].clone()], vec![]);
    let (_, numerical) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            design,
            &dispersions,
            &numeric,
            model_frame,
        )
        .unwrap();
    let (_, selected) = builder
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame(
            &counts,
            design,
            &dispersions,
            &named,
            model_frame,
        )
        .unwrap();
    assert_eq!(numerical.rows, selected.rows);
}
