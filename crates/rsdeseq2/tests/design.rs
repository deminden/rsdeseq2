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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
    ];
    assert!(expanded_additive_factor_design(&[]).is_err());
    assert!(expanded_additive_factor_design(&duplicate).is_err());

    let misaligned = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &short_batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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

    let removed_group = expanded_formula_design(
        "~ condition + batch + dose - (batch + dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_group, expected_removed);
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
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
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
    }];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let design = expanded_formula_design(
        "~ condition + log(dose) + log2(dose) + log10(dose) + sqrt(dose) + scale(dose)",
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
        assert_eq!(row[5], value.sqrt());
        assert!((row[6] - ((value - mean) / sd)).abs() < 1e-12);
    }
    assert_eq!(
        design.numeric_covariates,
        vec![
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_log10".to_string(),
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
        "~ scale(dose, center=FALSE) + scale(dose, scale=FALSE) + scale(dose, center=FALSE, scale=FALSE)",
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
fn expanded_formula_design_validates_unsupported_terms() {
    let condition = vec!["A".to_string(), "B".to_string()];
    let batch = vec!["X".to_string(), "Y".to_string()];
    let dose = [0.0, 1.0];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
        },
    ];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    assert!(expanded_formula_design("condition + dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + missing", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + condition", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ condition + condition:condition", &factors, &numeric).is_err()
    );
    assert!(expanded_formula_design("~ condition - missing", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition -", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ (condition + batch", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + (batch - dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose^x)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose^32)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose + 1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log2(dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=NaN)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=0)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=-1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design(
        "~ scale(dose)",
        &factors,
        &[ExpandedNumericSpec {
            name: "dose",
            values: &[1.0, 1.0],
        }]
    )
    .is_err());
    assert!(expanded_formula_design("~ poly(dose, 2)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, x, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 32, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(missing, 2, raw=TRUE)", &factors, &numeric).is_err());
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
