use super::*;

#[test]
fn expanded_additive_design_with_interactions_builds_factor_pair_terms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let interactions = [ExpandedFactorInteractionSpec {
        left_factor: "condition",
        right_factor: "batch",
    }];

    let expanded =
        expanded_additive_design_with_interactions(&factors, &[], &interactions).unwrap();

    assert_eq!(
        expanded.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "batchX".to_string(),
            "batchY".to_string(),
            "conditionA:batchX".to_string(),
            "conditionA:batchY".to_string(),
            "conditionB:batchX".to_string(),
            "conditionB:batchY".to_string(),
        ]
    );
    assert_eq!(
        expanded.expanded_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        expanded.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
            "condition_B_vs_A:batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        expanded.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 1.0, 1.0,
        ]
    );
    assert_eq!(
        expanded.coefficient_groups,
        vec![vec![0], vec![2], vec![4], vec![8]]
    );
    assert_eq!(expanded.interactions, vec!["condition:batch".to_string()]);
}

#[test]
fn expanded_additive_design_with_all_interactions_builds_numeric_interactions() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let time = [1.0, 1.0, 2.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];
    let factor_numeric = [ExpandedFactorNumericInteractionSpec {
        factor: "condition",
        numeric: "dose",
    }];
    let numeric_interactions = [ExpandedNumericInteractionSpec {
        left_numeric: "dose",
        right_numeric: "time",
    }];

    let expanded = expanded_additive_design_with_all_interactions(
        &factors,
        &numeric,
        &[],
        &factor_numeric,
        &numeric_interactions,
    )
    .unwrap();

    assert_eq!(
        expanded.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "dose".to_string(),
            "time".to_string(),
            "conditionA:dose".to_string(),
            "conditionB:dose".to_string(),
            "dose:time".to_string(),
        ]
    );
    assert_eq!(
        expanded.expanded_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0,
        ]
    );
    assert_eq!(
        expanded.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose".to_string(),
            "time".to_string(),
            "condition_B_vs_A:dose".to_string(),
            "dose:time".to_string(),
        ]
    );
    assert_eq!(
        expanded.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, 2.0, 0.0, 0.0, //
            1.0, 1.0, 1.0, 2.0, 1.0, 2.0,
        ]
    );
    assert_eq!(
        expanded.coefficient_groups,
        vec![vec![0], vec![2], vec![3], vec![4], vec![6], vec![7]]
    );
    assert_eq!(
        expanded.factor_numeric_interactions,
        vec!["condition:dose".to_string()]
    );
    assert_eq!(expanded.numeric_interactions, vec!["dose:time".to_string()]);
}

#[test]
fn expanded_additive_design_with_interactions_validates_specs() {
    let condition = vec!["A".to_string(), "B".to_string()];
    let batch = vec!["X".to_string(), "Y".to_string()];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];

    assert!(
        expanded_additive_design_with_interactions(
            &factors,
            &[],
            &[ExpandedFactorInteractionSpec {
                left_factor: "condition",
                right_factor: "missing",
            }]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design_with_interactions(
            &factors,
            &[],
            &[ExpandedFactorInteractionSpec {
                left_factor: "condition",
                right_factor: "condition",
            }]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design_with_interactions(
            &factors,
            &[],
            &[
                ExpandedFactorInteractionSpec {
                    left_factor: "condition",
                    right_factor: "batch",
                },
                ExpandedFactorInteractionSpec {
                    left_factor: "batch",
                    right_factor: "condition",
                },
            ]
        )
        .is_err()
    );

    let dose = [0.0, 1.0];
    let time = [1.0, 2.0];
    let numeric = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];
    assert!(
        expanded_additive_design_with_all_interactions(
            &factors,
            &numeric,
            &[],
            &[ExpandedFactorNumericInteractionSpec {
                factor: "condition",
                numeric: "missing",
            }],
            &[],
        )
        .is_err()
    );
    assert!(
        expanded_additive_design_with_all_interactions(
            &factors,
            &numeric,
            &[],
            &[],
            &[ExpandedNumericInteractionSpec {
                left_numeric: "dose",
                right_numeric: "dose",
            }],
        )
        .is_err()
    );
    assert!(
        expanded_additive_design_with_all_interactions(
            &factors,
            &numeric,
            &[],
            &[],
            &[
                ExpandedNumericInteractionSpec {
                    left_numeric: "dose",
                    right_numeric: "time",
                },
                ExpandedNumericInteractionSpec {
                    left_numeric: "time",
                    right_numeric: "dose",
                },
            ],
        )
        .is_err()
    );
}

#[test]
fn expanded_additive_design_validates_numeric_covariates() {
    let dose = [0.0, 1.0, 2.0];
    let bad = [0.0, f64::NAN, 2.0];
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "dose",
                values: &dose,
            }]
        )
        .is_ok()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "",
                values: &dose,
            }]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "dose",
                values: &[],
            }]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "dose",
                values: &bad,
            }]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[
                ExpandedNumericSpec {
                    name: "dose",
                    values: &dose,
                },
                ExpandedNumericSpec {
                    name: "dose",
                    values: &dose,
                },
            ]
        )
        .is_err()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "dose",
                values: &[0.0, 1.0],
            }]
        )
        .is_ok()
    );
    assert!(
        expanded_additive_design(
            &[],
            &[ExpandedNumericSpec {
                name: "Intercept",
                values: &dose,
            }]
        )
        .is_err()
    );
}

#[test]
fn expanded_formula_design_builds_supported_pairwise_terms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let time = [1.0, 1.0, 2.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];

    let parsed = expanded_formula_design(
        "~ condition + dose + time + condition:dose + dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    let direct = expanded_additive_design_with_all_interactions(
        &factors,
        &numeric,
        &[],
        &[ExpandedFactorNumericInteractionSpec {
            factor: "condition",
            numeric: "dose",
        }],
        &[ExpandedNumericInteractionSpec {
            left_numeric: "dose",
            right_numeric: "time",
        }],
    )
    .unwrap();

    assert_eq!(parsed, direct);
}

#[test]
fn expanded_formula_design_supports_intercept_only_formula() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let from_factor = expanded_formula_design("~ 1", &factors, &[]).unwrap();
    assert_eq!(
        from_factor.expanded_design.coefficient_names().unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(
        from_factor.standard_design.coefficient_names().unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(from_factor.expanded_design.matrix().as_slice(), &[1.0; 4]);
    assert_eq!(from_factor.standard_design.matrix().as_slice(), &[1.0; 4]);
    assert_eq!(from_factor.coefficient_groups, vec![vec![0]]);

    let from_numeric = expanded_formula_design("~ 1", &[], &numeric).unwrap();
    assert_eq!(from_numeric, from_factor);

    assert!(expanded_formula_design("~ 1", &[], &[]).is_err());
    assert!(expanded_formula_design("~ 0", &factors, &[]).is_err());
}

#[test]
fn expanded_formula_design_from_model_frame_infers_and_overrides_references() {
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
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: vec![
                    "Y".to_string(),
                    "Y".to_string(),
                    "X".to_string(),
                    "X".to_string(),
                ],
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0],
        }],
    };
    let inferred_condition = vec![
        "B".to_string(),
        "A".to_string(),
        "B".to_string(),
        "A".to_string(),
    ];
    let batch = vec![
        "Y".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "X".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &inferred_condition,
            reference: "B",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    assert_eq!(model_frame.n_samples().unwrap(), 4);
    model_frame.validate().unwrap();
    assert_eq!(
        model_frame.resolved_factor_reference("condition").unwrap(),
        Some("B")
    );
    assert_eq!(
        model_frame.resolved_factor_reference("batch").unwrap(),
        Some("X")
    );
    assert_eq!(
        model_frame.resolved_factor_reference("missing").unwrap(),
        None
    );
    assert_eq!(
        model_frame
            .resolved_factor_reference_by_alias("condition")
            .unwrap(),
        Some("B")
    );
    assert_eq!(
        model_frame.resolved_factor_references().unwrap(),
        vec![
            ResolvedFormulaFactorReference {
                factor: "condition",
                reference: "B",
                levels: None,
            },
            ResolvedFormulaFactorReference {
                factor: "batch",
                reference: "X",
                levels: None,
            },
        ]
    );
    assert_eq!(
        model_frame
            .resolved_numeric_covariate("dose")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "dose",
            values: &[0.0, 1.0, 2.0, 3.0],
        }
    );
    assert_eq!(
        model_frame.resolved_numeric_covariate("missing").unwrap(),
        None
    );
    assert_eq!(
        model_frame.resolved_numeric_covariates().unwrap(),
        vec![ResolvedFormulaNumericCovariate {
            name: "dose",
            values: &[0.0, 1.0, 2.0, 3.0],
        }]
    );

    let from_model_frame =
        expanded_formula_design_from_model_frame("~ condition * batch + dose", &model_frame)
            .unwrap();
    let direct = expanded_formula_design("~ condition * batch + dose", &factors, &numeric).unwrap();

    assert_eq!(from_model_frame, direct);
    assert_eq!(
        from_model_frame
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "condition_A_vs_B".to_string(),
            "batch_Y_vs_X".to_string(),
            "dose".to_string(),
            "condition_A_vs_B:batch_Y_vs_X".to_string(),
        ]
    );
}
