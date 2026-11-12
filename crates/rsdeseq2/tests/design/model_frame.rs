use super::*;

#[test]
fn formula_model_frame_resolves_r_cleaned_factor_reference_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec![
                    "T cell".to_string(),
                    "B cell".to_string(),
                    "T cell".to_string(),
                ],
                levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "if".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string(), "A".to_string()],
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "TRUE".to_string(),
                sample_levels: vec!["off".to_string(), "on".to_string(), "off".to_string()],
                levels: None,
                reference: Some("off".to_string()),
            },
        ],
        numeric_covariates: Vec::new(),
    };

    assert_eq!(
        model_frame
            .resolved_factor_reference_by_alias("cell.type")
            .unwrap(),
        Some("T cell")
    );
    assert_eq!(
        model_frame
            .resolved_factor_reference_by_alias("if.")
            .unwrap(),
        Some("A")
    );
    assert_eq!(
        model_frame
            .resolved_factor_reference_by_alias("TRUE.")
            .unwrap(),
        Some("off")
    );
    assert_eq!(
        model_frame
            .resolved_factor_reference_by_alias("missing")
            .unwrap(),
        None
    );
}

#[test]
fn formula_model_frame_factor_reference_aliases_prefer_exact_and_reject_ambiguous() {
    let exact = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "cell.type".to_string(),
                sample_levels: vec!["X".to_string(), "Y".to_string()],
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert_eq!(
        exact
            .resolved_factor_reference_by_alias("cell.type")
            .unwrap(),
        Some("X")
    );

    let ambiguous = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "cell-type".to_string(),
                sample_levels: vec!["X".to_string(), "Y".to_string()],
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert!(
        ambiguous
            .resolved_factor_reference_by_alias("cell.type")
            .is_err()
    );
}

#[test]
fn formula_model_frame_numeric_aliases_prefer_exact_and_reject_ambiguous() {
    let exact = FormulaModelFrame {
        factors: Vec::new(),
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose value".to_string(),
                values: vec![1.0, 2.0],
            },
            FormulaNumericColumn {
                name: "dose.value".to_string(),
                values: vec![3.0, 4.0],
            },
            FormulaNumericColumn {
                name: "if".to_string(),
                values: vec![5.0, 6.0],
            },
        ],
    };
    assert_eq!(
        exact
            .resolved_numeric_covariate_by_alias("dose.value")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "dose.value",
            values: &[3.0, 4.0],
        }
    );
    assert_eq!(
        exact
            .resolved_numeric_covariate_by_alias("if.")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "if",
            values: &[5.0, 6.0],
        }
    );
    assert_eq!(
        exact
            .resolved_numeric_covariate_by_alias("missing")
            .unwrap(),
        None
    );

    let ambiguous = FormulaModelFrame {
        factors: Vec::new(),
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose value".to_string(),
                values: vec![1.0, 2.0],
            },
            FormulaNumericColumn {
                name: "dose-value".to_string(),
                values: vec![3.0, 4.0],
            },
        ],
    };
    assert!(
        ambiguous
            .resolved_numeric_covariate_by_alias("dose.value")
            .is_err()
    );
}

#[test]
fn formula_model_frame_reference_accepts_cleaned_level_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "B cell".to_string(),
                "A cell".to_string(),
            ],
            levels: Some(vec!["A cell".to_string(), "B cell".to_string()]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    assert_eq!(
        model_frame.resolved_factor_reference("cell type").unwrap(),
        Some("A cell")
    );
    assert_eq!(
        model_frame.resolved_factor_references().unwrap(),
        vec![ResolvedFormulaFactorReference {
            factor: "cell type",
            reference: "A cell",
            levels: Some(&["A cell".to_string(), "B cell".to_string()][..]),
        }]
    );

    let design = expanded_formula_design_from_model_frame("~ `cell type`", &model_frame).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "cell type_B cell_vs_A cell".to_string(),
        ]
    );
}

#[test]
fn formula_model_frame_reference_alias_can_resolve_observed_levels_without_declared_levels() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "B cell".to_string(),
                "A cell".to_string(),
            ],
            levels: None,
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    assert_eq!(
        model_frame.resolved_factor_reference("cell type").unwrap(),
        Some("A cell")
    );

    let design = expanded_formula_design_from_model_frame("~ `cell type`", &model_frame).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "cell type_B cell_vs_A cell".to_string(),
        ]
    );
}

#[test]
fn formula_model_frame_reference_rejects_ambiguous_cleaned_level_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "A-cell".to_string(),
                "B cell".to_string(),
            ],
            levels: Some(vec![
                "A cell".to_string(),
                "A-cell".to_string(),
                "B cell".to_string(),
            ]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    assert!(
        model_frame
            .validate()
            .unwrap_err()
            .to_string()
            .contains("reference level 'A.cell' resolves ambiguously")
    );

    let declared_collision = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "A-cell".to_string(),
                "B cell".to_string(),
            ],
            levels: Some(vec![
                "A cell".to_string(),
                "A-cell".to_string(),
                "B cell".to_string(),
            ]),
            reference: Some("B cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        declared_collision
            .validate()
            .unwrap_err()
            .to_string()
            .contains("declared level 'A-cell' resolves ambiguously after R-style cleanup")
    );
}

#[test]
fn expanded_formula_design_from_model_frame_uses_declared_factor_levels() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "B".to_string(),
                "A".to_string(),
                "C".to_string(),
                "B".to_string(),
            ],
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let design = expanded_formula_design_from_model_frame("~ condition", &model_frame).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_C_vs_A".to_string(),
        ]
    );
    assert_eq!(
        design.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "conditionC".to_string(),
        ]
    );
    assert_eq!(
        design.factor_levels,
        vec![vec!["A".to_string(), "B".to_string(), "C".to_string()]]
    );
}

#[test]
fn expanded_formula_design_from_model_frame_supports_quoted_column_names() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "T cell".to_string(),
                "B cell".to_string(),
                "T cell".to_string(),
                "B cell".to_string(),
            ],
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose value".to_string(),
                values: vec![1.0, 2.0, 3.0, 4.0],
            },
            FormulaNumericColumn {
                name: "size/value".to_string(),
                values: vec![0.5, 1.0, 2.0, 4.0],
            },
        ],
    };

    let design = expanded_formula_design_with_offsets_from_model_frame(
        "~ `cell type` * I(`dose value` + 1) + log2(`size/value`) + offset(`dose value`)",
        &model_frame,
    )
    .unwrap();

    assert_eq!(design.offsets, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(
        design.design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "cell type_B cell_vs_T cell".to_string(),
            "dose value_plus_1".to_string(),
            "size/value_log2".to_string(),
            "cell type_B cell_vs_T cell:dose value_plus_1".to_string(),
        ]
    );
    assert_eq!(
        design.design.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 2.0, -1.0, 0.0, //
            1.0, 1.0, 3.0, 0.0, 3.0, //
            1.0, 0.0, 4.0, 1.0, 0.0, //
            1.0, 1.0, 5.0, 2.0, 5.0,
        ]
    );

    let subtracted = expanded_formula_design_from_model_frame(
        "~ `cell type` + `dose value` - `dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        subtracted.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "cell type_B cell_vs_T cell".to_string(),
        ]
    );
}

#[test]
fn expanded_formula_design_from_model_frame_supports_dot_main_effects() {
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec![
                    "B".to_string(),
                    "A".to_string(),
                    "B".to_string(),
                    "A".to_string(),
                ],
                levels: Some(vec!["A".to_string(), "B".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec![
                    "T".to_string(),
                    "T".to_string(),
                    "B".to_string(),
                    "B".to_string(),
                ],
                levels: Some(vec!["T".to_string(), "B".to_string()]),
                reference: None,
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose value".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0],
        }],
    };

    let dot = expanded_formula_design_from_model_frame("~ .", &model_frame).unwrap();
    let explicit = expanded_formula_design_from_model_frame(
        "~ condition + `cell type` + `dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(dot, explicit);
    assert_eq!(
        dot.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "cell type_B_vs_T".to_string(),
            "dose value".to_string(),
        ]
    );

    let removed =
        expanded_formula_design_from_model_frame("~ . - `dose value` - condition", &model_frame)
            .unwrap();
    let expected = expanded_formula_design_from_model_frame("~ `cell type`", &model_frame).unwrap();
    assert_eq!(removed, expected);

    let star_dot =
        expanded_formula_design_from_model_frame("~ condition * .", &model_frame).unwrap();
    let star_explicit = expanded_formula_design_from_model_frame(
        "~ condition + `cell type` + `dose value` + condition:`cell type` + condition:`dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(star_dot, star_explicit);

    let interaction_dot =
        expanded_formula_design_from_model_frame("~ condition:.", &model_frame).unwrap();
    let interaction_explicit = expanded_formula_design_from_model_frame(
        "~ condition + condition:`cell type` + condition:`dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(interaction_dot, interaction_explicit);

    let nested_dot =
        expanded_formula_design_from_model_frame("~ condition / .", &model_frame).unwrap();
    let nested_explicit = expanded_formula_design_from_model_frame(
        "~ condition + condition:`cell type` + condition:`dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(nested_dot, nested_explicit);

    let star_dot_removed =
        expanded_formula_design_from_model_frame("~ condition * . - condition:.", &model_frame)
            .unwrap();
    let star_dot_removed_explicit =
        expanded_formula_design_from_model_frame("~ `cell type` + `dose value`", &model_frame)
            .unwrap();
    assert_eq!(star_dot_removed, star_dot_removed_explicit);

    let nested_dot_removed =
        expanded_formula_design_from_model_frame("~ condition / . - condition:.", &model_frame)
            .unwrap();
    let nested_dot_removed_explicit =
        expanded_formula_design_from_model_frame("~ 1", &model_frame).unwrap();
    assert_eq!(nested_dot_removed, nested_dot_removed_explicit);

    let all_dot_removed =
        expanded_formula_design_from_model_frame("~ condition * . - condition * .", &model_frame)
            .unwrap();
    let intercept_only = expanded_formula_design_from_model_frame("~ 1", &model_frame).unwrap();
    assert_eq!(all_dot_removed, intercept_only);
}
