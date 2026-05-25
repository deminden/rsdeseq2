use rsdeseq2::prelude::*;

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
