use rsdeseq2::design::formula_has_offset_terms;
use rsdeseq2::prelude::*;
use std::assert_matches;

#[test]
fn intercept_only_design_has_named_all_ones_column() {
    let design = DesignMatrix::intercept_only(3).unwrap();

    assert_eq!(design.n_samples(), 3);
    assert_eq!(design.n_coefficients(), 1);
    assert_eq!(
        design.coefficient_names().unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(design.matrix().as_slice(), &[1.0, 1.0, 1.0]);
    assert!(design.is_full_rank().unwrap());
}

#[test]
fn design_matrix_reports_full_rank() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    assert_eq!(design.rank().unwrap(), 2);
    assert!(design.is_full_rank().unwrap());
    design.validate_full_rank("test").unwrap();
}

#[test]
fn design_matrix_resolves_coefficient_names() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    assert_eq!(design.coefficient_index("Intercept").unwrap(), 0);
    assert_eq!(design.coefficient_index("condition_B_vs_A").unwrap(), 1);
    assert!(design.coefficient_index("missing").is_err());

    let unnamed = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();
    assert!(unnamed.coefficient_index("Intercept").is_err());
}

#[test]
fn design_matrix_sample_rows_accept_legacy_and_new_ranges() {
    let design = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "dose".into()]),
    )
    .unwrap();

    assert_eq!(design.sample_rows(1..3).unwrap(), &[1.0, 1.0, 1.0, 2.0]);
    assert_eq!(
        design
            .sample_rows(core::range::Range { start: 0, end: 2 })
            .unwrap(),
        &[1.0, 0.0, 1.0, 1.0]
    );
    assert_matches!(
        design
            .sample_rows(core::range::Range { start: 0, end: 4 })
            .unwrap_err(),
        DeseqError::InvalidDimensions { .. }
    );
}

#[test]
fn design_matrix_detects_dependent_columns() {
    let design = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "duplicate".into()]),
    )
    .unwrap();

    assert_eq!(design.rank().unwrap(), 1);
    assert!(!design.is_full_rank().unwrap());
    let error = design.validate_full_rank("test").unwrap_err();
    assert!(error.to_string().contains("not full rank"));
}

#[test]
fn design_matrix_reports_zero_columns_in_full_rank_error() {
    let design = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0,
        ],
        Some(vec!["Intercept".into(), "empty_group".into()]),
    )
    .unwrap();

    let error = design.validate_full_rank("test").unwrap_err();
    let message = error.to_string();
    assert!(message.contains("not full rank"));
    assert!(message.contains("empty_group"));
}

#[test]
fn original_design_full_rank_error_reports_unnamed_zero_column_index() {
    let design = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0,
        ],
        None,
    )
    .unwrap();

    let error = design.validate_full_rank("test").unwrap_err();
    let message = error.to_string();

    assert!(message.contains("not full rank"));
    assert!(message.contains("zero columns: 1"));
}

#[test]
fn design_rank_validates_tolerance() {
    let design = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();

    assert!(design.rank_with_tolerance(f64::NAN).is_err());
    assert!(design.rank_with_tolerance(-1.0).is_err());
}

#[test]
fn expanded_factor_design_builds_expanded_and_standard_surfaces() {
    let levels = ["A", "B", "A", "C"];
    let expanded = expanded_factor_design("condition", &levels, "A").unwrap();

    assert_eq!(
        expanded.levels,
        vec!["A".to_string(), "B".to_string(), "C".to_string()]
    );
    assert_eq!(
        expanded.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "conditionC".to_string(),
        ]
    );
    assert_eq!(
        expanded.expanded_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        expanded.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_C_vs_A".to_string(),
        ]
    );
    assert_eq!(
        expanded.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0,
        ]
    );
    assert_eq!(expanded.coefficient_groups, vec![vec![0], vec![2], vec![3]]);
}

#[test]
fn expanded_factor_design_validates_inputs() {
    assert!(expanded_factor_design("", &["A", "B"], "A").is_err());
    assert!(expanded_factor_design("condition", &["A", "B"], "").is_err());
    assert!(expanded_factor_design("condition", &["A", ""], "A").is_err());
    assert!(expanded_factor_design("condition", &["B", "C"], "A").is_err());
    let empty: [&str; 0] = [];
    assert!(expanded_factor_design("condition", &empty, "A").is_err());
}

#[test]
fn expanded_additive_factor_design_builds_multiple_factor_surfaces() {
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

    let expanded = expanded_additive_factor_design(&factors).unwrap();

    assert_eq!(
        expanded.factor_levels,
        vec![
            vec!["A".to_string(), "B".to_string()],
            vec!["X".to_string(), "Y".to_string()],
        ]
    );
    assert_eq!(
        expanded.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "batchX".to_string(),
            "batchY".to_string(),
        ]
    );
    assert_eq!(
        expanded.expanded_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, 1.0, //
            1.0, 0.0, 1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        expanded.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        expanded.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ]
    );
    assert_eq!(expanded.coefficient_groups, vec![vec![0], vec![2], vec![4]]);
    assert!(expanded.numeric_covariates.is_empty());
    assert!(expanded.interactions.is_empty());
}

#[test]
fn expanded_additive_factor_design_validates_inputs() {
    let condition = vec!["A".to_string(), "B".to_string()];
    let short_batch = vec!["X".to_string()];
    let duplicate = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
    ];
    assert!(expanded_additive_factor_design(&[]).is_err());
    assert!(expanded_additive_factor_design(&duplicate).is_err());

    let misaligned = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &short_batch,
            reference: "X",
            levels: None,
        },
    ];
    assert!(expanded_additive_factor_design(&misaligned).is_err());
}

#[test]
fn expanded_additive_design_includes_numeric_covariates_unchanged() {
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

    let expanded = expanded_additive_design(&factors, &numeric).unwrap();

    assert_eq!(
        expanded.factor_levels,
        vec![vec!["A".to_string(), "B".to_string()]]
    );
    assert_eq!(expanded.numeric_covariates, vec!["dose".to_string()]);
    assert_eq!(
        expanded.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA".to_string(),
            "conditionB".to_string(),
            "dose".to_string(),
        ]
    );
    assert_eq!(
        expanded.expanded_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 1.0,
        ]
    );
    assert_eq!(
        expanded.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose".to_string(),
        ]
    );
    assert_eq!(
        expanded.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ]
    );
    assert_eq!(expanded.coefficient_groups, vec![vec![0], vec![2], vec![3]]);
}

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

    assert!(expanded_additive_design_with_interactions(
        &factors,
        &[],
        &[ExpandedFactorInteractionSpec {
            left_factor: "condition",
            right_factor: "missing",
        }]
    )
    .is_err());
    assert!(expanded_additive_design_with_interactions(
        &factors,
        &[],
        &[ExpandedFactorInteractionSpec {
            left_factor: "condition",
            right_factor: "condition",
        }]
    )
    .is_err());
    assert!(expanded_additive_design_with_interactions(
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
    .is_err());

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
    assert!(expanded_additive_design_with_all_interactions(
        &factors,
        &numeric,
        &[],
        &[ExpandedFactorNumericInteractionSpec {
            factor: "condition",
            numeric: "missing",
        }],
        &[],
    )
    .is_err());
    assert!(expanded_additive_design_with_all_interactions(
        &factors,
        &numeric,
        &[],
        &[],
        &[ExpandedNumericInteractionSpec {
            left_numeric: "dose",
            right_numeric: "dose",
        }],
    )
    .is_err());
    assert!(expanded_additive_design_with_all_interactions(
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
    .is_err());
}

#[test]
fn expanded_additive_design_validates_numeric_covariates() {
    let dose = [0.0, 1.0, 2.0];
    let bad = [0.0, f64::NAN, 2.0];
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        }]
    )
    .is_ok());
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "",
            values: &dose,
        }]
    )
    .is_err());
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[],
        }]
    )
    .is_err());
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "dose",
            values: &bad,
        }]
    )
    .is_err());
    assert!(expanded_additive_design(
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
    .is_err());
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[0.0, 1.0],
        }]
    )
    .is_ok());
    assert!(expanded_additive_design(
        &[],
        &[ExpandedNumericSpec {
            name: "Intercept",
            values: &dose,
        }]
    )
    .is_err());
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

#[test]
fn expanded_formula_design_supports_relevel_factor_transform() {
    let condition = vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "B".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let condition_levels = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: Some(&condition_levels),
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];

    let design =
        expanded_formula_design("~ relevel(condition, ref=\"B\") + batch", &factors, &[]).unwrap();

    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(condition, ref = \"B\")_A_vs_B".to_string(),
            "relevel(condition, ref = \"B\")_C_vs_B".to_string(),
            "batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        design.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(condition, ref = \"B\")B".to_string(),
            "relevel(condition, ref = \"B\")A".to_string(),
            "relevel(condition, ref = \"B\")C".to_string(),
            "batchX".to_string(),
            "batchY".to_string(),
        ]
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: condition,
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let from_model_frame =
        expanded_formula_design_from_model_frame("~ relevel(x=condition, ref='C')", &model_frame)
            .unwrap();
    assert_eq!(
        from_model_frame
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(condition, ref = \"C\")_A_vs_C".to_string(),
            "relevel(condition, ref = \"C\")_B_vs_C".to_string(),
        ]
    );
}

#[test]
fn expanded_formula_design_supports_factor_identity_transforms() {
    let condition = vec![
        "A".to_string(),
        "B".to_string(),
        "A".to_string(),
        "B".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let condition_levels = vec!["A".to_string(), "B".to_string()];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: Some(&condition_levels),
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];

    let design =
        expanded_formula_design("~ factor(condition) + as.factor(x=batch)", &factors, &[]).unwrap();

    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)_B_vs_A".to_string(),
            "as.factor(batch)_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        design.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)A".to_string(),
            "factor(condition)B".to_string(),
            "as.factor(batch)X".to_string(),
            "as.factor(batch)Y".to_string(),
        ]
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: condition,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("B".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let from_model_frame =
        expanded_formula_design_from_model_frame("~ as.factor(condition)", &model_frame).unwrap();
    assert_eq!(
        from_model_frame
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "as.factor(condition)_A_vs_B".to_string(),
        ]
    );
}

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
    assert!(ambiguous
        .resolved_factor_reference_by_alias("cell.type")
        .is_err());
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
        "~ condition:`cell type` + condition:`dose value`",
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
    let star_dot_removed_explicit = expanded_formula_design_from_model_frame(
        "~ condition + `cell type` + `dose value`",
        &model_frame,
    )
    .unwrap();
    assert_eq!(star_dot_removed, star_dot_removed_explicit);

    let nested_dot_removed =
        expanded_formula_design_from_model_frame("~ condition / . - condition:.", &model_frame)
            .unwrap();
    let nested_dot_removed_explicit =
        expanded_formula_design_from_model_frame("~ condition", &model_frame).unwrap();
    assert_eq!(nested_dot_removed, nested_dot_removed_explicit);

    let all_dot_removed =
        expanded_formula_design_from_model_frame("~ condition * . - condition * .", &model_frame)
            .unwrap();
    let intercept_only = expanded_formula_design_from_model_frame("~ 1", &model_frame).unwrap();
    assert_eq!(all_dot_removed, intercept_only);
}

#[test]
fn expanded_formula_design_simplifies_duplicate_formula_terms() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let simplified = expanded_formula_design(
        "~ condition + condition + dose + dose + condition:batch + batch:condition + condition:dose + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    let explicit = expanded_formula_design(
        "~ condition + dose + condition:batch + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(simplified, explicit);

    let star_simplified =
        expanded_formula_design("~ condition * batch + condition:batch", &factors, &numeric)
            .unwrap();
    let star_explicit = expanded_formula_design("~ condition * batch", &factors, &numeric).unwrap();
    assert_eq!(star_simplified, star_explicit);

    let higher_order = expanded_formula_design(
        "~ condition:batch:dose + dose:condition:batch",
        &factors,
        &numeric,
    )
    .unwrap();
    let higher_order_explicit =
        expanded_formula_design("~ condition:batch:dose", &factors, &numeric).unwrap();
    assert_eq!(higher_order, higher_order_explicit);
}

#[test]
fn expanded_formula_design_from_model_frame_keeps_unused_declared_reference() {
    let condition = vec!["B".to_string(), "C".to_string(), "B".to_string()];
    let levels = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: Some(&levels),
    }];
    let direct = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_C_vs_A".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0]
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: condition,
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let from_model_frame =
        expanded_formula_design_from_model_frame("~ condition", &model_frame).unwrap();
    assert_eq!(from_model_frame, direct);
    assert_eq!(model_frame.n_samples().unwrap(), 3);
    assert_eq!(
        model_frame.resolved_factor_reference("condition").unwrap(),
        Some("A")
    );
    assert_eq!(
        model_frame.resolved_factor_references().unwrap(),
        vec![ResolvedFormulaFactorReference {
            factor: "condition",
            reference: "A",
            levels: Some(&levels),
        }]
    );
}

#[test]
fn expanded_formula_design_from_model_frame_validates_declared_factor_levels() {
    let duplicate_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string(), "A".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &duplicate_level)
            .unwrap_err()
            .to_string()
            .contains("duplicate declared level")
    );

    let missing_sample_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string(), "C".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &missing_sample_level)
            .unwrap_err()
            .to_string()
            .contains("declared levels")
    );

    let missing_reference = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("C".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &missing_reference)
            .unwrap_err()
            .to_string()
            .contains("reference level")
    );
}

#[test]
fn expanded_formula_design_from_model_frame_handles_offsets_and_validation() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string(), "A".to_string()],
            levels: None,
            reference: Some("A".to_string()),
        }],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose".to_string(),
                values: vec![1.0, 2.0, 4.0],
            },
            FormulaNumericColumn {
                name: "offset_log".to_string(),
                values: vec![0.1, 0.2, 0.3],
            },
        ],
    };

    let parsed = expanded_formula_design_with_offsets_from_model_frame(
        "~ condition + log2(dose) + offset(offset_log)",
        &model_frame,
    )
    .unwrap();
    assert_eq!(parsed.offsets, vec![0.1, 0.2, 0.3]);
    assert_eq!(
        parsed.design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_log2".to_string(),
        ]
    );

    let named_transforms = expanded_formula_design_with_offsets_from_model_frame(
        "~ condition + I(x=dose + 1) + log(x=dose, base=2) + scale(x=dose, center=FALSE, scale=FALSE) + poly(x=dose, degree=2) + offset(I(x=offset_log + dose))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(named_transforms.offsets, vec![1.1, 2.2, 4.3],);
    assert_eq!(
        named_transforms
            .design
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_plus_1".to_string(),
            "dose_log_base_2".to_string(),
            "dose_identity".to_string(),
            "poly(dose, 2)1".to_string(),
            "poly(dose, 2)2".to_string(),
        ]
    );
    for (sample, dose) in [1.0_f64, 2.0, 4.0].iter().copied().enumerate() {
        let row = named_transforms
            .design
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[2], dose + 1.0);
        assert!((row[3] - dose.log2()).abs() < 1e-12);
        assert_eq!(row[4], dose);
    }

    let duplicate = FormulaModelFrame {
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
    assert!(expanded_formula_design_from_model_frame("~ condition", &duplicate).is_err());

    let bad_reference = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: Some("C".to_string()),
        }],
        numeric_covariates: vec![],
    };
    assert!(expanded_formula_design_from_model_frame("~ condition", &bad_reference).is_err());
}

#[test]
fn expanded_formula_design_expands_star_terms() {
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

    let parsed = expanded_formula_design("~ condition * batch", &factors, &[]).unwrap();
    let direct = expanded_additive_design_with_all_interactions(
        &factors,
        &[],
        &[ExpandedFactorInteractionSpec {
            left_factor: "condition",
            right_factor: "batch",
        }],
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(parsed, direct);
}

#[test]
fn expanded_formula_design_accepts_interactions_without_main_effects() {
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
    let dose = [0.0, 1.0, 0.0, 1.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let interaction_only = expanded_formula_design("~ condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        interaction_only
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA:batchX".to_string(),
            "conditionA:batchY".to_string(),
            "conditionB:batchX".to_string(),
            "conditionB:batchY".to_string(),
        ]
    );
    assert_eq!(
        interaction_only.standard_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        interaction_only.coefficient_groups,
        vec![vec![0], vec![1], vec![2], vec![3], vec![4]]
    );

    let one_main = expanded_formula_design("~ condition + condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        one_main.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        one_main.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        one_main.coefficient_groups,
        vec![vec![0], vec![2], vec![4], vec![6]]
    );

    let factor_numeric =
        expanded_formula_design("~ condition + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(
        factor_numeric.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:dose".to_string(),
            "conditionB:dose".to_string(),
        ]
    );
    assert_eq!(
        factor_numeric.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0,
        ]
    );
}

#[test]
fn expanded_formula_design_supports_intercept_removal() {
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

    let no_intercept = expanded_formula_design("~ 0 + condition + batch", &factors, &[]).unwrap();
    assert_eq!(
        no_intercept.standard_design.coefficient_names().unwrap(),
        &[
            "conditionA".to_string(),
            "conditionB".to_string(),
            "batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        no_intercept.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 1.0,
        ]
    );
    assert_eq!(
        no_intercept.coefficient_groups,
        vec![vec![0], vec![1], vec![3]]
    );

    let minus_one = expanded_formula_design("~ condition - 1", &factors, &[]).unwrap();
    assert_eq!(
        minus_one.standard_design.coefficient_names().unwrap(),
        &["conditionA".to_string(), "conditionB".to_string()]
    );
    assert_eq!(
        minus_one.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, //
            1.0, 0.0, //
            0.0, 1.0, //
            0.0, 1.0,
        ]
    );

    let restored = expanded_formula_design("~ 0 + 1 + condition", &factors, &[]).unwrap();
    let restored_from_minus_one =
        expanded_formula_design("~ condition - 1 + 1", &factors, &[]).unwrap();
    let with_intercept = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    assert_eq!(restored, with_intercept);
    assert_eq!(restored_from_minus_one, with_intercept);

    let removed_again = expanded_formula_design("~ 0 + 1 + condition - 1", &factors, &[]).unwrap();
    assert_eq!(removed_again, minus_one);

    let interaction_only = expanded_formula_design("~ 0 + condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        interaction_only
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "conditionA:batchX".to_string(),
            "conditionA:batchY".to_string(),
            "conditionB:batchX".to_string(),
            "conditionB:batchY".to_string(),
        ]
    );
    assert_eq!(
        interaction_only.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, //
            0.0, 0.0, 1.0, 0.0, //
            0.0, 0.0, 0.0, 1.0,
        ]
    );
}

#[test]
fn expanded_formula_design_supports_pairwise_nested_terms() {
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
    let dose = [0.0, 1.0, 0.0, 1.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let nested = expanded_formula_design("~ condition / batch", &factors, &[]).unwrap();
    let expanded = expanded_formula_design("~ condition + condition:batch", &factors, &[]).unwrap();
    assert_eq!(nested, expanded);
    assert_eq!(
        nested.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );

    let nested_numeric = expanded_formula_design("~ condition / dose", &factors, &numeric).unwrap();
    let expanded_numeric =
        expanded_formula_design("~ condition + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(nested_numeric, expanded_numeric);

    let no_intercept = expanded_formula_design("~ 0 + condition / batch", &factors, &[]).unwrap();
    assert_eq!(
        no_intercept.standard_design.coefficient_names().unwrap(),
        &[
            "conditionA".to_string(),
            "conditionB".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        no_intercept.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 1.0,
        ]
    );

    let nested_in = expanded_formula_design("~ batch %in% condition", &factors, &[]).unwrap();
    let nested_in_explicit = expanded_formula_design("~ condition:batch", &factors, &[]).unwrap();
    assert_eq!(nested_in, nested_in_explicit);

    let nested_in_with_outer =
        expanded_formula_design("~ condition + batch %in% condition", &factors, &[]).unwrap();
    assert_eq!(nested_in_with_outer, nested);

    let nested_in_numeric =
        expanded_formula_design("~ dose %in% condition", &factors, &numeric).unwrap();
    let nested_in_numeric_explicit =
        expanded_formula_design("~ condition:dose", &factors, &numeric).unwrap();
    assert_eq!(nested_in_numeric, nested_in_numeric_explicit);

    assert!(expanded_formula_design("~ condition %in%", &factors, &[]).is_err());
    assert!(
        expanded_formula_design("~ condition %in% batch %in% dose", &factors, &numeric).is_err()
    );
}

#[test]
fn expanded_formula_design_supports_three_way_terms() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let direct = expanded_formula_design(
        "~ condition + batch + dose + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
            "dose".to_string(),
            "condition_B_vs_A:batch_Y_vs_X:dose".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 2.0, 0.0, //
            1.0, 1.0, 1.0, 3.0, 3.0,
        ]
    );
    assert_eq!(
        direct.higher_order_interactions,
        vec!["condition:batch:dose".to_string()]
    );

    let star = expanded_formula_design("~ condition * batch * dose", &factors, &numeric).unwrap();
    assert!(star
        .standard_design
        .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose")
        .is_ok());
    assert_eq!(
        star.higher_order_interactions,
        vec!["condition:batch:dose".to_string()]
    );

    let nested = expanded_formula_design("~ condition / batch / dose", &factors, &numeric).unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:batch + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested, expanded_nested);
}

#[test]
fn expanded_formula_design_supports_arbitrary_order_terms() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
    let time = [2.0, 3.0, 5.0, 7.0];
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

    let direct = expanded_formula_design(
        "~ condition + batch + dose + time + condition:batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
            "dose".to_string(),
            "time".to_string(),
            "condition_B_vs_A:batch_Y_vs_X:dose:time".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 2.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 3.0, 0.0, //
            1.0, 1.0, 0.0, 2.0, 5.0, 0.0, //
            1.0, 1.0, 1.0, 3.0, 7.0, 21.0,
        ]
    );
    assert_eq!(
        direct.higher_order_interactions,
        vec!["condition:batch:dose:time".to_string()]
    );

    let star =
        expanded_formula_design("~ condition * batch * dose * time", &factors, &numeric).unwrap();
    assert!(star
        .standard_design
        .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose:time")
        .is_ok());
    assert!(star
        .higher_order_interactions
        .iter()
        .any(|interaction| interaction == "condition:batch:dose:time"));

    let nested =
        expanded_formula_design("~ condition / batch / dose / time", &factors, &numeric).unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:batch + condition:batch:dose + condition:batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested, expanded_nested);
}

#[test]
fn expanded_formula_design_supports_term_subtraction() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
    let time = [2.0, 3.0, 5.0, 7.0];
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

    let removed_main =
        expanded_formula_design("~ condition + batch + dose - batch", &factors, &numeric).unwrap();
    let expected_main = expanded_formula_design("~ condition + dose", &factors, &numeric).unwrap();
    assert_eq!(removed_main, expected_main);

    let removed_pairwise =
        expanded_formula_design("~ condition * batch - condition:batch", &factors, &[]).unwrap();
    let expected_pairwise = expanded_formula_design("~ condition + batch", &factors, &[]).unwrap();
    assert_eq!(removed_pairwise, expected_pairwise);

    let removed_higher = expanded_formula_design(
        "~ condition * batch * dose * time - condition:batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert!(removed_higher
        .standard_design
        .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose:time")
        .is_err());
    assert!(!removed_higher
        .higher_order_interactions
        .iter()
        .any(|interaction| interaction == "condition:batch:dose:time"));

    let removed_nested = expanded_formula_design(
        "~ condition / batch / dose - condition / batch",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_nested =
        expanded_formula_design("~ condition:batch:dose", &factors, &numeric).unwrap();
    assert_eq!(removed_nested, expected_nested);
}

#[test]
fn expanded_formula_design_supports_parenthesized_groups() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let grouped_star =
        expanded_formula_design("~ (condition + batch) * dose", &factors, &numeric).unwrap();
    let expanded_star = expanded_formula_design(
        "~ condition + batch + dose + condition:dose + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(grouped_star, expanded_star);

    let grouped_interaction =
        expanded_formula_design("~ (condition + batch):dose", &factors, &numeric).unwrap();
    let expanded_interaction =
        expanded_formula_design("~ condition:dose + batch:dose", &factors, &numeric).unwrap();
    assert_eq!(grouped_interaction, expanded_interaction);

    let grouped_nested =
        expanded_formula_design("~ (condition + batch) / dose", &factors, &numeric).unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:dose + batch + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(grouped_nested, expanded_nested);

    let grouped_nested_in =
        expanded_formula_design("~ (batch + dose) %in% condition", &factors, &numeric).unwrap();
    let expanded_nested_in =
        expanded_formula_design("~ condition:batch + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(grouped_nested_in, expanded_nested_in);

    let removed_group = expanded_formula_design(
        "~ condition + batch + dose - (batch + dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_group, expected_removed);

    let removed_nested_in = expanded_formula_design(
        "~ condition + condition:batch + condition:dose - (batch + dose) %in% condition",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed_nested_in =
        expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_nested_in, expected_removed_nested_in);
}

#[test]
fn expanded_formula_design_supports_nested_parenthesized_groups() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
    let time = [2.0, 3.0, 5.0, 7.0];
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

    let nested_star =
        expanded_formula_design("~ (condition + (batch + dose)) * time", &factors, &numeric)
            .unwrap();
    let expanded_star = expanded_formula_design(
        "~ condition + batch + dose + time + condition:time + batch:time + dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested_star, expanded_star);

    let nested_interaction =
        expanded_formula_design("~ condition:(batch + (dose + time))", &factors, &numeric).unwrap();
    let expanded_interaction = expanded_formula_design(
        "~ condition:batch + condition:dose + condition:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested_interaction, expanded_interaction);

    let nested_nested =
        expanded_formula_design("~ condition / (batch + (dose + time))", &factors, &numeric)
            .unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:batch + condition:dose + condition:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested_nested, expanded_nested);

    let removed_nested_group = expanded_formula_design(
        "~ condition + batch + dose + time - (batch + (dose + time))",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_nested_group, expected_removed);
}

#[test]
fn expanded_formula_design_supports_parenthesized_power_terms() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
    let time = [2.0, 3.0, 5.0, 7.0];
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

    let pairwise =
        expanded_formula_design("~ (condition + batch + dose)^2", &factors, &numeric).unwrap();
    let explicit_pairwise = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(pairwise, explicit_pairwise);

    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: condition.clone(),
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: batch.clone(),
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose".to_string(),
                values: dose.to_vec(),
            },
            FormulaNumericColumn {
                name: "time".to_string(),
                values: time.to_vec(),
            },
        ],
    };
    let model_frame_pairwise =
        expanded_formula_design_from_model_frame("~ (condition + batch + dose)^2", &model_frame)
            .unwrap();
    assert_eq!(model_frame_pairwise, explicit_pairwise);

    let higher_order =
        expanded_formula_design("~ (condition + batch + dose + time)^3", &factors, &numeric)
            .unwrap();
    let explicit_higher_order = expanded_formula_design(
        "~ condition + batch + dose + time + condition:batch + condition:dose + condition:time + batch:dose + batch:time + dose:time + condition:batch:dose + condition:batch:time + condition:dose:time + batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(higher_order, explicit_higher_order);

    let first_order =
        expanded_formula_design("~ (condition + batch + dose)^1", &factors, &numeric).unwrap();
    let explicit_first_order =
        expanded_formula_design("~ condition + batch + dose", &factors, &numeric).unwrap();
    assert_eq!(first_order, explicit_first_order);

    assert!(expanded_formula_design("~ (condition + batch)^0", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ (condition + batch)^2.5", &factors, &numeric).is_err());
}

#[test]
fn expanded_formula_design_supports_dot_power_terms() {
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
    let dose = [0.0, 1.0, 2.0, 3.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let dot_pairwise = expanded_formula_design("~ .^2", &factors, &numeric).unwrap();
    let explicit_pairwise = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(dot_pairwise, explicit_pairwise);

    let dot_parenthesized = expanded_formula_design("~ (.)^2", &factors, &numeric).unwrap();
    assert_eq!(dot_parenthesized, explicit_pairwise);

    let removed_pairwise =
        expanded_formula_design("~ .^2 - .^2 + condition", &factors, &numeric).unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_pairwise, expected_removed);

    let dot_higher = expanded_formula_design("~ .^3", &factors, &numeric).unwrap();
    let explicit_higher = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(dot_higher, explicit_higher);

    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: condition.clone(),
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: batch.clone(),
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: dose.to_vec(),
        }],
    };
    let model_frame_dot = expanded_formula_design_from_model_frame("~ .^2", &model_frame).unwrap();
    assert_eq!(model_frame_dot, explicit_pairwise);

    assert!(expanded_formula_design("~ .^0", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ .^two", &factors, &numeric).is_err());
}

#[test]
fn expanded_formula_design_supports_numeric_power_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
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

    let design =
        expanded_formula_design("~ condition + dose + I(dose^2)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose".to_string(),
            "dose_pow_2".to_string(),
        ]
    );
    assert_eq!(
        design.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, //
            1.0, 1.0, 2.0, 4.0, //
            1.0, 1.0, 3.0, 9.0,
        ]
    );
    assert_eq!(
        design.numeric_covariates,
        vec!["dose".to_string(), "dose_pow_2".to_string()]
    );

    let interaction =
        expanded_formula_design("~ condition * I(dose^2)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_pow_2")
        .is_ok());

    let removed =
        expanded_formula_design("~ condition + I(dose^2) - I(dose^2)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_numeric_identity_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design = expanded_formula_design("~ condition + I(dose)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_identity".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        assert_eq!(
            design.standard_design.matrix().row(sample).unwrap()[2],
            value
        );
    }
    assert_eq!(design.numeric_covariates, vec!["dose_identity".to_string()]);

    let interaction = expanded_formula_design("~ condition * I(dose)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_identity")
        .is_ok());

    let removed =
        expanded_formula_design("~ condition + I(dose) - I(dose)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_signed_numeric_identity_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design = expanded_formula_design("~ I(+dose) + I(-dose)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_neg".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value);
        assert_eq!(row[2], -value);
    }

    let interaction =
        expanded_formula_design("~ condition * I(-dose)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_neg")
        .is_ok());

    let removed =
        expanded_formula_design("~ condition + I(-dose) - I(-dose)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_numeric_scalar_arithmetic_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design = expanded_formula_design(
        "~ I(dose + 1) + I(dose - 1) + I(dose * 2) + I(dose / 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_1".to_string(),
            "dose_minus_1".to_string(),
            "dose_times_2".to_string(),
            "dose_div_2".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value + 1.0);
        assert_eq!(row[2], value - 1.0);
        assert_eq!(row[3], value * 2.0);
        assert_eq!(row[4], value / 2.0);
    }

    let interaction =
        expanded_formula_design("~ condition * I(dose + 1)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_plus_1")
        .is_ok());

    let removed = expanded_formula_design(
        "~ condition + I(dose + 1) - I(dose + 1)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_scalar_left_numeric_arithmetic_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
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

    let design = expanded_formula_design(
        "~ I(1 + dose) + I(1 - dose) + I(2 * dose) + I(2 / dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_1".to_string(),
            "dose_rminus_1".to_string(),
            "dose_times_2".to_string(),
            "dose_rdiv_2".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], 1.0 + value);
        assert_eq!(row[2], 1.0 - value);
        assert_eq!(row[3], 2.0 * value);
        assert_eq!(row[4], 2.0 / value);
    }

    let interaction =
        expanded_formula_design("~ condition * I(1 + dose)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_plus_1")
        .is_ok());
}

#[test]
fn expanded_formula_design_supports_numeric_binary_arithmetic_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
    let time = [1.0_f64, 2.0, 3.0, 4.0];
    let dose_over_time = [1.0_f64, 1.0, 4.0 / 3.0, 2.0];
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
        ExpandedNumericSpec {
            name: "dose/time",
            values: &dose_over_time,
        },
    ];

    let design = expanded_formula_design(
        "~ I(dose + time) + I(dose - time) + I(dose * time) + I(dose / time)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_time".to_string(),
            "dose_minus_time".to_string(),
            "dose_times_time".to_string(),
            "dose_div_time".to_string(),
        ]
    );
    for sample in 0..dose.len() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], dose[sample] + time[sample]);
        assert_eq!(row[2], dose[sample] - time[sample]);
        assert_eq!(row[3], dose[sample] * time[sample]);
        assert_eq!(row[4], dose[sample] / time[sample]);
    }

    let interaction =
        expanded_formula_design("~ condition * I(dose + time)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_plus_time")
        .is_ok());

    let named_x = expanded_formula_design(
        "~ I(x=dose) + I(x=dose + 1) + I(x=`dose/time`)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_x.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_plus_1".to_string(),
            "dose/time_identity".to_string(),
        ]
    );
    for (sample, dose_value) in dose.iter().copied().enumerate() {
        let row = named_x.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], dose_value);
        assert_eq!(row[2], dose_value + 1.0);
        assert_eq!(row[3], numeric[2].values[sample]);
    }

    let removed = expanded_formula_design(
        "~ condition + I(dose + time) - I(dose + time)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_numeric_function_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 4.0, 9.0, 16.0];
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

    let design = expanded_formula_design(
        "~ condition + log(dose) + log2(dose) + log10(dose) + log1p(dose) + sqrt(dose) + scale(dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_log10".to_string(),
            "dose_log1p".to_string(),
            "dose_sqrt".to_string(),
            "dose_scale".to_string(),
        ]
    );
    let mean = dose.iter().sum::<f64>() / dose.len() as f64;
    let sd = (dose
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / (dose.len() as f64 - 1.0))
        .sqrt();
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[2], value.ln());
        assert_eq!(row[3], value.log2());
        assert_eq!(row[4], value.log10());
        assert_eq!(row[5], value.ln_1p());
        assert_eq!(row[6], value.sqrt());
        assert!((row[7] - ((value - mean) / sd)).abs() < 1e-12);
    }
    assert_eq!(
        design.numeric_covariates,
        vec![
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_log10".to_string(),
            "dose_log1p".to_string(),
            "dose_sqrt".to_string(),
            "dose_scale".to_string(),
        ]
    );

    let interaction =
        expanded_formula_design("~ condition * log2(dose)", &factors, &numeric).unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_log2")
        .is_ok());

    let named_function_args = expanded_formula_design(
        "~ log(x=dose) + log(dose, base=2) + log(x=dose, base=10) + log(dose, 4) + log2(x=dose) + sqrt(x=dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_function_args
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_log".to_string(),
            "dose_log_base_2".to_string(),
            "dose_log_base_10".to_string(),
            "dose_log_base_4".to_string(),
            "dose_log2".to_string(),
            "dose_sqrt".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = named_function_args
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value.ln());
        assert!((row[2] - value.log2()).abs() < 1e-12);
        assert!((row[3] - value.log10()).abs() < 1e-12);
        assert!((row[4] - value.log(4.0)).abs() < 1e-12);
        assert_eq!(row[5], value.log2());
        assert_eq!(row[6], value.sqrt());
    }

    let removed =
        expanded_formula_design("~ condition + sqrt(dose) - sqrt(dose)", &factors, &numeric)
            .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);

    let scaled_interaction =
        expanded_formula_design("~ condition * scale(dose)", &factors, &numeric).unwrap();
    assert!(scaled_interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_scale")
        .is_ok());

    let scaled_options = expanded_formula_design(
        "~ scale(dose, center=FALSE) + scale(dose, scale=FALSE) + scale(dose, center=FALSE, scale=FALSE) + scale(dose, FALSE, FALSE) + scale(dose, F, scale=T)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        scaled_options.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_scale_uncentered".to_string(),
            "dose_center".to_string(),
            "dose_identity".to_string(),
        ]
    );
    let rms =
        (dose.iter().map(|value| value * value).sum::<f64>() / (dose.len() as f64 - 1.0)).sqrt();
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = scaled_options.standard_design.matrix().row(sample).unwrap();
        assert!((row[1] - value / rms).abs() < 1e-12);
        assert!((row[2] - (value - mean)).abs() < 1e-12);
        assert_eq!(row[3], value);
    }

    let positional_scaled_options = expanded_formula_design(
        "~ scale(dose, FALSE, FALSE) + scale(dose, F, T)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        positional_scaled_options
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_scale_uncentered".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = positional_scaled_options
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value);
        assert!((row[2] - value / rms).abs() < 1e-12);
    }

    let named_x_scaled_options = expanded_formula_design(
        "~ scale(x=dose) + scale(x=dose, center=FALSE, scale=FALSE) + scale(center=FALSE, x=dose, scale=TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_x_scaled_options
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_scale".to_string(),
            "dose_identity".to_string(),
            "dose_scale_uncentered".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = named_x_scaled_options
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert!((row[1] - ((value - mean) / sd)).abs() < 1e-12);
        assert_eq!(row[2], value);
        assert!((row[3] - value / rms).abs() < 1e-12);
    }

    let scaled_constants = expanded_formula_design(
        "~ scale(dose, center=5, scale=FALSE) + scale(dose, center=FALSE, scale=2) + scale(dose, center=5, scale=2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        scaled_constants
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_centered".to_string(),
            "dose_scaled".to_string(),
            "dose_centered_scaled".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = scaled_constants
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value - 5.0);
        assert_eq!(row[2], value / 2.0);
        assert_eq!(row[3], (value - 5.0) / 2.0);
    }

    let parenthesized_arguments = expanded_formula_design(
        "~ scale(dose, center=(5), scale=FALSE) + scale(dose, center=FALSE, scale=(2))",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        parenthesized_arguments
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_centered".to_string(),
            "dose_scaled".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = parenthesized_arguments
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value - 5.0);
        assert_eq!(row[2], value / 2.0);
    }
}

#[test]
fn expanded_formula_design_transform_parser_handles_nested_parentheses_and_quotes() {
    let numeric_with_parens = [1.0_f64, 2.0, 4.0, 8.0];
    let factors = [];
    let numeric = [ExpandedNumericSpec {
        name: "dose(value)",
        values: &numeric_with_parens,
    }];

    let design = expanded_formula_design(
        "~ log2(`dose(value)`) + I(`dose(value)` + 1)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose(value)_log2".to_string(),
            "dose(value)_plus_1".to_string(),
        ]
    );
    assert_eq!(
        resolve_coefficient_index(&design.standard_design, "dose.value._log2").unwrap(),
        1
    );
    assert_eq!(
        resolve_coefficient_index(&design.standard_design, "dose.value._plus_1").unwrap(),
        2
    );
    for (sample, value) in numeric_with_parens.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value.log2());
        assert_eq!(row[2], value + 1.0);
    }

    let literal_log = [1.0_f64, 3.0, 5.0, 7.0];
    let literal_poly = [2.0_f64, 4.0, 6.0, 8.0];
    let literal_identity = [0.5_f64, 1.5, 2.5, 3.5];
    let equals_name = [1.0_f64, 2.0, 3.0, 4.0];
    let literal_numeric = [
        ExpandedNumericSpec {
            name: "log2(dose)",
            values: &literal_log,
        },
        ExpandedNumericSpec {
            name: "poly(dose, 2)",
            values: &literal_poly,
        },
        ExpandedNumericSpec {
            name: "I(dose + 1)",
            values: &literal_identity,
        },
        ExpandedNumericSpec {
            name: "dose=a",
            values: &equals_name,
        },
    ];
    let literal = expanded_formula_design(
        "~ `log2(dose)` + `poly(dose, 2)` + `I(dose + 1)` + log2(`dose=a`) + scale(`dose=a`, FALSE, FALSE)",
        &factors,
        &literal_numeric,
    )
    .unwrap();
    assert_eq!(
        literal.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "log2(dose)".to_string(),
            "poly(dose, 2)".to_string(),
            "I(dose + 1)".to_string(),
            "dose=a_log2".to_string(),
            "dose=a_identity".to_string(),
        ]
    );
    for sample in 0..literal_log.len() {
        let row = literal.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], literal_log[sample]);
        assert_eq!(row[2], literal_poly[sample]);
        assert_eq!(row[3], literal_identity[sample]);
        assert_eq!(row[4], equals_name[sample].log2());
        assert_eq!(row[5], equals_name[sample]);
    }
}

#[test]
fn expanded_formula_design_supports_raw_polynomial_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design =
        expanded_formula_design("~ condition + poly(dose, 3, raw=TRUE)", &factors, &numeric)
            .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
            "dose_poly_3".to_string(),
        ]
    );
    assert_eq!(
        design.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 2.0, 4.0, 8.0, //
            1.0, 1.0, 3.0, 9.0, 27.0,
        ]
    );
    assert_eq!(
        design.numeric_covariates,
        vec![
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
            "dose_poly_3".to_string(),
        ]
    );

    let interaction = expanded_formula_design(
        "~ condition:poly(dose, degree=2, raw = TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        interaction.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA:dose_poly_1".to_string(),
            "conditionB:dose_poly_1".to_string(),
            "conditionA:dose_poly_2".to_string(),
            "conditionB:dose_poly_2".to_string(),
        ]
    );

    let treatment_interaction = expanded_formula_design(
        "~ condition + poly(dose, degree=2, raw = TRUE) + condition:poly(dose, degree=2, raw = TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert!(treatment_interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_poly_1")
        .is_ok());
    assert!(treatment_interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:dose_poly_2")
        .is_ok());

    let named_order =
        expanded_formula_design("~ poly(dose, raw = TRUE, degree = 2)", &factors, &numeric)
            .unwrap();
    assert_eq!(
        named_order.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
        ]
    );

    let named_x =
        expanded_formula_design("~ poly(x=dose, degree=2, raw=TRUE)", &factors, &numeric).unwrap();
    assert_eq!(named_x, named_order);

    let removed = expanded_formula_design(
        "~ condition + poly(dose, 2, raw=TRUE) - poly(dose, 2, raw=TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_orthogonal_polynomial_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design =
        expanded_formula_design("~ condition + poly(dose, 3)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "poly(dose, 3)1".to_string(),
            "poly(dose, 3)2".to_string(),
            "poly(dose, 3)3".to_string(),
        ]
    );
    let expected = [
        -0.6708203932499369,
        0.5,
        -0.22360679774997888,
        -0.22360679774997896,
        -0.5,
        0.6708203932499369,
        0.22360679774997896,
        -0.5,
        -0.6708203932499369,
        0.6708203932499369,
        0.5,
        0.22360679774997896,
    ];
    for (row_idx, row) in design
        .standard_design
        .matrix()
        .as_slice()
        .chunks_exact(5)
        .enumerate()
    {
        for (observed, expected) in row[2..].iter().zip(&expected[row_idx * 3..row_idx * 3 + 3]) {
            assert!((observed - expected).abs() < 1e-12);
        }
    }

    let raw_false =
        expanded_formula_design("~ poly(dose, degree=2, raw=FALSE)", &factors, &numeric).unwrap();
    let default = expanded_formula_design("~ poly(dose, degree=2)", &factors, &numeric).unwrap();
    assert_eq!(raw_false, default);

    let named_x = expanded_formula_design("~ poly(x=dose, degree=2)", &factors, &numeric).unwrap();
    assert_eq!(named_x, default);
    let reordered_named =
        expanded_formula_design("~ poly(degree=2, x=dose)", &factors, &numeric).unwrap();
    assert_eq!(reordered_named, default);
    let simple =
        expanded_formula_design("~ poly(dose, 2, simple=TRUE)", &factors, &numeric).unwrap();
    assert_eq!(simple, default);

    let interaction = expanded_formula_design(
        "~ condition + poly(dose, 2) + condition:poly(dose, 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:poly(dose, 2)1")
        .is_ok());
    assert!(interaction
        .standard_design
        .coefficient_index("condition_B_vs_A:poly(dose, 2)2")
        .is_ok());

    let removed = expanded_formula_design(
        "~ condition + poly(dose, 2) - poly(dose, 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected_removed);

    let repeated_dose = [0.0_f64, 0.0, 1.0, 2.0];
    let repeated_numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &repeated_dose,
    }];
    let repeated = expanded_formula_design("~ poly(dose, 2)", &[], &repeated_numeric).unwrap();
    let repeated_expected = [
        -0.45226701686664544,
        0.21320071635561041,
        -0.45226701686664544,
        0.21320071635561041,
        0.15075567228888181,
        -0.8528028654224417,
        0.753778361444409,
        0.4264014327112208,
    ];
    for (observed, expected) in repeated
        .standard_design
        .matrix()
        .as_slice()
        .chunks_exact(3)
        .flat_map(|row| row[1..].iter())
        .zip(repeated_expected)
    {
        assert!((observed - expected).abs() < 1e-12);
    }
}

#[test]
fn expanded_formula_design_with_offsets_extracts_formula_offsets() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
    let exposure = [0.5, 1.0, 1.5, 2.0];
    let batch_offset = [0.1, 0.2, 0.3, 0.4];
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
            name: "exposure",
            values: &exposure,
        },
        ExpandedNumericSpec {
            name: "batch_offset",
            values: &batch_offset,
        },
    ];

    let with_offsets = expanded_formula_design_with_offsets(
        "~ condition + dose + offset(exposure) + offset(batch_offset)",
        &factors,
        &numeric,
    )
    .unwrap();
    let without_offsets =
        expanded_formula_design("~ condition + dose", &factors, &numeric).unwrap();
    assert_eq!(with_offsets.design, without_offsets);
    assert_eq!(with_offsets.offsets, vec![0.6, 1.2, 1.8, 2.4]);

    let intercept_only =
        expanded_formula_design_with_offsets("~ offset(exposure)", &[], &numeric).unwrap();
    assert_eq!(
        intercept_only
            .design
            .standard_design
            .coefficient_names()
            .unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(intercept_only.offsets, exposure);
}

#[test]
fn expanded_formula_design_with_offsets_supports_numeric_transform_offsets() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
    let exposure = [0.5_f64, 1.0, 1.5, 2.0];
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
            name: "exposure value",
            values: &exposure,
        },
    ];

    let transformed = expanded_formula_design_with_offsets(
        "~ condition + offset(log2(dose)) + offset(I(x=dose + `exposure value`))",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_design = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(transformed.design, expected_design);
    assert_eq!(
        transformed.offsets,
        dose.iter()
            .zip(exposure.iter())
            .map(|(dose, exposure)| dose.log2() + dose + exposure)
            .collect::<Vec<_>>()
    );

    let quoted =
        expanded_formula_design_with_offsets("~ offset(`exposure value`)", &[], &numeric).unwrap();
    assert_eq!(quoted.offsets, exposure);

    let scaled = expanded_formula_design_with_offsets(
        "~ offset(scale(dose, center=FALSE, scale=FALSE))",
        &[],
        &numeric,
    )
    .unwrap();
    assert_eq!(scaled.offsets, dose);

    let named_transform_offsets = expanded_formula_design_with_offsets(
        "~ offset(log(x=dose, base=2)) + offset(scale(x=`exposure value`, center=FALSE, scale=FALSE))",
        &[],
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_transform_offsets.offsets,
        dose.iter()
            .zip(exposure.iter())
            .map(|(dose, exposure)| dose.log2() + exposure)
            .collect::<Vec<_>>()
    );

    assert!(expanded_formula_design_with_offsets(
        "~ offset(poly(dose, 2, raw=TRUE))",
        &[],
        &numeric
    )
    .is_err());
}

#[test]
fn formula_offset_detection_is_syntax_level() {
    assert!(!formula_has_offset_terms("~ condition + offset_like").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(exposure)").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(log2(exposure))").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(I(exposure + 1))").unwrap());
    assert!(formula_has_offset_terms("~ condition - offset(exposure)").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset( exposure )").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(exposure").is_err());
    assert!(formula_has_offset_terms("condition + offset(exposure)").is_err());
}

#[test]
fn expanded_formula_design_validates_unsupported_terms() {
    let condition = vec!["A".to_string(), "B".to_string()];
    let batch = vec!["X".to_string(), "Y".to_string()];
    let dose = [0.0, 1.0];
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
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];
    let numeric_with_zero_time = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &[0.0, 1.0],
        },
    ];

    assert!(expanded_formula_design("condition + dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + missing", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ condition + condition:condition", &factors, &numeric).is_err()
    );
    assert!(expanded_formula_design("~ condition - missing", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition -", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ (condition + batch", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + (batch - dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose^x)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose^32)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose + missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(y=dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose / 0)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(1 / dose)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ I(dose / time)", &factors, &numeric_with_zero_time).is_err()
    );
    assert!(expanded_formula_design("~ log(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose + 1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=-2)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(x=dose, x=dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log2(dose, base=2)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ sqrt(dose, extra=1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design(
        "~ log1p(dose)",
        &factors,
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[-1.0, 0.0],
        }]
    )
    .is_err());
    assert!(expanded_formula_design("~ log2(dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=NaN)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=0)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=-1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design(
        "~ scale(dose, center=FALSE, center=TRUE)",
        &factors,
        &numeric
    )
    .is_err());
    assert!(
        expanded_formula_design("~ scale(dose, scale=FALSE, scale=TRUE)", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ scale(x=dose, x=dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design(
        "~ scale(dose)",
        &factors,
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[1.0, 1.0],
        }]
    )
    .is_err());
    assert!(expanded_formula_design(
        "~ poly(dose, 2)",
        &factors,
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[0.0, 0.0],
        }]
    )
    .is_err());
    assert!(expanded_formula_design("~ poly(dose, x, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 32, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(missing, 2, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, coefs=1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, simple=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, raw=maybe)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ poly(x=dose, x=dose, degree=2)", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design("~ poly(dose, degree=2, degree=3)", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design("~ poly(dose, 2, raw=TRUE, raw=FALSE)", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ factor(dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ factor(missing)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ factor(condition, levels=c('A','B'))", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ factor(x=condition, x=batch)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ as.factor(condition, extra=1)", &factors, &numeric).is_err()
    );
    assert!(expanded_formula_design("~ relevel(dose, ref='A')", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ relevel(condition)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ relevel(condition, ref='Z')", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design_with_offsets("~ condition - offset(dose)", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design_with_offsets(
        "~ condition + offset(missing)",
        &factors,
        &numeric
    )
    .is_err());
    assert!(expanded_formula_design_with_offsets(
        "~ condition + offset(dose + 1)",
        &factors,
        &numeric
    )
    .is_err());
    assert!(
        expanded_formula_design("~ condition * batch * dose * missing", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ condition:batch:dose:missing", &factors, &numeric).is_err());
}
