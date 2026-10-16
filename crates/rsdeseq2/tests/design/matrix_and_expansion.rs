use super::*;

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
