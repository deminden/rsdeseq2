use super::*;

#[test]
fn model_frame_factor_level_contrast_infers_and_overrides_reference() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "C".to_string(),
            ],
            levels: None,
            reference: Some("B".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0],
        }],
    };

    let inferred_request = ResultsContrast::character("condition", "C", "A");
    let inferred = factor_level_contrast_from_model_frame(&inferred_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(inferred.factor, "condition");
    assert_eq!(inferred.numerator, "C");
    assert_eq!(inferred.denominator, "A");
    assert_eq!(inferred.reference, Some("B"));
    assert_eq!(
        inferred.sample_levels,
        &model_frame.factors[0].sample_levels
    );

    let explicit_request = ResultsContrast::character_with_reference("condition", "C", "B", "A");
    let explicit = factor_level_contrast_from_model_frame(&explicit_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(explicit.reference, Some("A"));

    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            &model_frame
        )
        .unwrap()
        .is_none()
    );
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::list(vec!["condition_C_vs_A".into()], Vec::new()),
            &model_frame,
        )
        .unwrap()
        .is_none()
    );
}

#[test]
fn model_frame_factor_level_contrast_uses_declared_level_reference() {
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

    let inferred_request = ResultsContrast::character("condition", "B", "A");
    let inferred = factor_level_contrast_from_model_frame(&inferred_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(inferred.reference, Some("A"));

    let explicit_request = ResultsContrast::character_with_reference("condition", "C", "B", "B");
    let explicit = factor_level_contrast_from_model_frame(&explicit_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(explicit.reference, Some("B"));
}

#[test]
fn model_frame_factor_level_contrast_allows_unused_declared_reference() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["B".to_string(), "C".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("condition", "C", "B");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "condition");
    assert_eq!(contrast.numerator, "C");
    assert_eq!(contrast.denominator, "B");
    assert_eq!(contrast.reference, Some("A"));
    assert_eq!(
        contrast.sample_levels,
        &model_frame.factors[0].sample_levels
    );

    let missing_reference = ResultsContrast::character_with_reference("condition", "C", "B", "D");
    assert!(factor_level_contrast_from_model_frame(&missing_reference, &model_frame).is_err());
}

#[test]
fn model_frame_factor_level_contrast_accepts_r_cleaned_factor_name_alias() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "T cell".to_string(),
                "B-cell".to_string(),
                "B-cell".to_string(),
                "T cell".to_string(),
            ],
            levels: Some(vec!["T cell".to_string(), "B-cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B-cell", "T cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "cell type");
    assert_eq!(contrast.reference, Some("T cell"));
    assert_eq!(
        contrast.sample_levels,
        &model_frame.factors[0].sample_levels
    );

    let relabeled_model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control factor(A)\", \"treated factor(B)\"))".to_string(),
            sample_levels: vec![
                "control factor(A)".to_string(),
                "treated factor(B)".to_string(),
                "control factor(A)".to_string(),
                "treated factor(B)".to_string(),
            ],
            levels: Some(vec![
                "control factor(A)".to_string(),
                "treated factor(B)".to_string(),
            ]),
            reference: Some("control factor(A)".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let relabeled_request = ResultsContrast::character(
        "factor.condition..levels...c..A....B....labels...c..control.factor.A.....treated.factor.B....",
        "treated factor(B)",
        "control factor(A)",
    );
    let relabeled_contrast =
        factor_level_contrast_from_model_frame(&relabeled_request, &relabeled_model_frame)
            .unwrap()
            .unwrap();
    assert_eq!(
        relabeled_contrast.factor,
        relabeled_model_frame.factors[0].name
    );
    assert_eq!(relabeled_contrast.reference, Some("control factor(A)"));
    assert_eq!(
        relabeled_contrast.sample_levels,
        &relabeled_model_frame.factors[0].sample_levels
    );

    let escaped_model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: r#"factor(condition, levels = c("A", "B"), labels = c("control \"quote\"", "treated \\ path"))"#.to_string(),
            sample_levels: vec![
                "control \"quote\"".to_string(),
                "treated \\ path".to_string(),
                "control \"quote\"".to_string(),
                "treated \\ path".to_string(),
            ],
            levels: Some(vec![
                "control \"quote\"".to_string(),
                "treated \\ path".to_string(),
            ]),
            reference: Some("control \"quote\"".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let escaped_request = ResultsContrast::character(
        r#"factor.condition..levels...c..A....B....labels...c..control...quote......treated....path..."#,
        "treated \\ path",
        "control \"quote\"",
    );
    let escaped_contrast =
        factor_level_contrast_from_model_frame(&escaped_request, &escaped_model_frame)
            .unwrap()
            .unwrap();
    assert_eq!(escaped_contrast.factor, escaped_model_frame.factors[0].name);
    assert_eq!(escaped_contrast.reference, Some("control \"quote\""));
    assert_eq!(
        escaped_contrast.sample_levels,
        &escaped_model_frame.factors[0].sample_levels
    );
}

#[test]
fn model_frame_factor_level_contrast_accepts_r_cleaned_level_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "T cell".to_string(),
                "B cell".to_string(),
                "B cell".to_string(),
                "T cell".to_string(),
            ],
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B.cell", "T.cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "cell type");
    assert_eq!(contrast.numerator, "B cell");
    assert_eq!(contrast.denominator, "T cell");
    assert_eq!(contrast.reference, Some("T cell"));
    assert_eq!(
        contrast.sample_levels,
        &model_frame.factors[0].sample_levels
    );
}

#[test]
fn model_frame_factor_level_contrast_prefers_exact_levels_before_cleaned_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "A.cell".to_string(),
                "B cell".to_string(),
                "B.cell".to_string(),
            ],
            levels: Some(vec![
                "A cell".to_string(),
                "A.cell".to_string(),
                "B cell".to_string(),
                "B.cell".to_string(),
            ]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B.cell", "A.cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.numerator, "B.cell");
    assert_eq!(contrast.denominator, "A.cell");
    assert_eq!(contrast.reference, Some("A.cell"));

    let cleaned_request = ResultsContrast::character("cell.type", "B.cell.", "A.cell.");
    assert!(
        factor_level_contrast_from_model_frame(&cleaned_request, &model_frame)
            .unwrap_err()
            .to_string()
            .contains("does not contain numerator level")
    );
}

#[test]
fn model_frame_factor_level_contrast_resolves_cleaned_reference_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "T cell".to_string(),
                "B cell".to_string(),
                "B cell".to_string(),
                "T cell".to_string(),
            ],
            levels: Some(vec![
                "A cell".to_string(),
                "T cell".to_string(),
                "B cell".to_string(),
            ]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B.cell", "T.cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.reference, Some("A cell"));

    let explicit =
        ResultsContrast::character_with_reference("cell.type", "B.cell", "T.cell", "A.cell");
    let explicit_contrast = factor_level_contrast_from_model_frame(&explicit, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(explicit_contrast.reference, Some("A cell"));
}

#[test]
fn model_frame_factor_level_contrast_resolves_cleaned_observed_reference_alias() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "B cell".to_string(),
                "B cell".to_string(),
                "A cell".to_string(),
            ],
            levels: None,
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B.cell", "A.cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "cell type");
    assert_eq!(contrast.numerator, "B cell");
    assert_eq!(contrast.denominator, "A cell");
    assert_eq!(contrast.reference, Some("A cell"));
}

#[test]
fn model_frame_factor_level_contrast_rejects_ambiguous_level_aliases() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "A cell".to_string(),
                "B cell".to_string(),
                "B-cell".to_string(),
                "A cell".to_string(),
            ],
            levels: Some(vec![
                "A cell".to_string(),
                "B cell".to_string(),
                "B-cell".to_string(),
            ]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("cell.type", "B.cell", "A.cell"),
            &model_frame,
        )
        .unwrap_err()
        .to_string()
        .contains("ambiguously")
    );
}

#[test]
fn model_frame_factor_level_contrast_validates_factor_metadata() {
    let duplicate = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec!["A".to_string(), "C".to_string()],
                levels: None,
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("condition", "B", "A"),
            &duplicate,
        )
        .unwrap_err()
        .to_string()
        .contains("more than once")
    );

    let missing_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("condition", "C", "A"),
            &missing_level,
        )
        .unwrap_err()
        .to_string()
        .contains("numerator level")
    );
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character_with_reference("condition", "B", "A", "C"),
            &missing_level,
        )
        .unwrap_err()
        .to_string()
        .contains("reference level")
    );

    let ambiguous_cleaned = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "cell-type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("cell.type", "B", "A"),
            &ambiguous_cleaned,
        )
        .unwrap_err()
        .to_string()
        .contains("ambiguously")
    );

    let cross_type_duplicate = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "condition".to_string(),
            values: vec![0.0, 1.0],
        }],
    };
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("condition", "B", "A"),
            &cross_type_duplicate,
        )
        .unwrap_err()
        .to_string()
        .contains("cannot be both a factor and numeric covariate")
    );

    let cross_type_alias = FormulaModelFrame {
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
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("cell type", "B", "A"),
            &cross_type_alias,
        )
        .unwrap_err()
        .to_string()
        .contains("resolves ambiguously after R-style cleanup")
    );

    let bad_numeric_shape = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: vec![0.0],
        }],
    };
    assert!(
        factor_level_contrast_from_model_frame(
            &ResultsContrast::character("condition", "B", "A"),
            &bad_numeric_shape,
        )
        .unwrap_err()
        .to_string()
        .contains("formula numeric covariate values")
    );
}
