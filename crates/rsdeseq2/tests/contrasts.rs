use rsdeseq2::prelude::*;

fn named_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "batch_Y_vs_X".into(),
        ]),
    )
    .unwrap()
}

#[test]
fn numeric_contrast_resolves_after_validation() {
    let design = named_design();
    let contrast = ContrastSpec::Numeric(vec![0.0, 1.0, -1.0]);
    assert_eq!(
        resolve_contrast(&design, &contrast).unwrap(),
        vec![0.0, 1.0, -1.0]
    );

    assert!(resolve_contrast(&design, &ContrastSpec::Numeric(vec![1.0])).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::Numeric(vec![0.0, 0.0, 0.0])).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::Numeric(vec![0.0, f64::NAN, 0.0])).is_err());
}

#[test]
fn coefficient_name_contrast_resolves_to_unit_vector() {
    let design = named_design();
    let contrast = ContrastSpec::coefficient_name("condition_B_vs_A");
    assert_eq!(
        resolve_contrast(&design, &contrast).unwrap(),
        vec![0.0, 1.0, 0.0]
    );

    assert!(resolve_contrast(&design, &ContrastSpec::coefficient_name("missing")).is_err());

    let unnamed = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();
    assert!(resolve_contrast(&unnamed, &ContrastSpec::coefficient_name("Intercept")).is_err());
}

#[test]
fn contrast_specs_expose_stable_result_metadata_labels() {
    let coefficient = ContrastSpec::coefficient_name("condition_B_vs_A");
    assert_eq!(coefficient.result_name(), "condition_B_vs_A");
    assert_eq!(coefficient.comparison(), "coefficient condition_B_vs_A");

    let factor = ContrastSpec::factor_level("condition", "B", "A");
    assert_eq!(factor.result_name(), "condition_B_vs_A");
    assert_eq!(
        factor.comparison(),
        "factor-level contrast: condition B vs A"
    );

    let list = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        0.5,
        -2.0,
    );
    assert_eq!(list.result_name(), "contrast");
    assert_eq!(
        list.comparison(),
        "coefficient list contrast: condition_B_vs_A at 0.5 vs batch_Y_vs_X at -2"
    );
}

#[test]
fn list_contrast_resolves_like_deseq2_name_lists() {
    let design = named_design();
    let contrast = ContrastSpec::list(vec!["condition_B_vs_A".into()], vec!["batch_Y_vs_X".into()]);
    assert_eq!(
        resolve_contrast(&design, &contrast).unwrap(),
        vec![0.0, 1.0, -1.0]
    );

    let weighted = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        0.5,
        -2.0,
    );
    assert_eq!(
        resolve_contrast(&design, &weighted).unwrap(),
        vec![0.0, 0.5, -2.0]
    );

    let numerator_only = ContrastSpec::list(vec!["condition_B_vs_A".into()], Vec::new());
    assert_eq!(
        resolve_contrast(&design, &numerator_only).unwrap(),
        vec![0.0, 1.0, 0.0]
    );

    let denominator_only =
        ContrastSpec::list_with_values(Vec::new(), vec!["batch_Y_vs_X".into()], 0.5, -0.5);
    assert_eq!(
        resolve_contrast(&design, &denominator_only).unwrap(),
        vec![0.0, 0.0, -0.5]
    );

    let overlap = ContrastSpec::list(
        vec!["condition_B_vs_A".into()],
        vec!["condition_B_vs_A".into()],
    );
    assert!(resolve_contrast(&design, &overlap).is_err());

    let empty = ContrastSpec::list(Vec::new(), Vec::new());
    assert!(resolve_contrast(&design, &empty).is_err());
}

#[test]
fn factor_level_contrast_resolves_standard_reference_shapes() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "condition_C_vs_A".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "B", "A")).unwrap(),
        vec![0.0, 1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "A", "B")).unwrap(),
        vec![0.0, -1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level_with_reference("condition", "C", "B", "A")
        )
        .unwrap(),
        vec![0.0, -1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_resolves_expanded_shapes() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec!["conditionA".into(), "conditionB".into()]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "B", "A")).unwrap(),
        vec![-1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_uses_r_like_make_names_for_candidates() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec!["Intercept".into(), "condition_B.1_vs_A.1".into()]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("condition", "B-1", "A-1")
        )
        .unwrap(),
        vec![0.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_validates_inputs() {
    let design = named_design();
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "A", "A")).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "C", "B")).is_err());
}

#[test]
fn contrast_all_zero_numeric_matches_deseq2_expanded_shape() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        6,
        vec![
            0, 0, 0, 0, 50, 60, //
            10, 12, 0, 0, 50, 60, //
            0, 0, 0, 0, 0, 0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 1.0,
        ],
        Some(vec![
            "conditionA".into(),
            "conditionB".into(),
            "conditionC".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        contrast_all_zero_numeric(&counts, &design, &[-1.0, 1.0, 0.0]).unwrap(),
        vec![true, false, true]
    );
    assert_eq!(
        contrast_all_zero_numeric(&counts, &design, &[0.0, 1.0, 0.0]).unwrap(),
        vec![false, false, false]
    );
}

#[test]
fn contrast_all_zero_numeric_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![0, 0]).unwrap();
    let design = DesignMatrix::from_row_major(2, 2, vec![1.0, 0.0, 0.0, 1.0], None).unwrap();
    let wrong_sample_design = DesignMatrix::from_row_major(1, 2, vec![1.0, 0.0], None).unwrap();

    assert!(contrast_all_zero_numeric(&counts, &design, &[1.0]).is_err());
    assert!(contrast_all_zero_numeric(&counts, &design, &[0.0, 0.0]).is_err());
    assert!(contrast_all_zero_numeric(&counts, &design, &[1.0, f64::NAN]).is_err());
    assert!(contrast_all_zero_numeric(&counts, &wrong_sample_design, &[1.0, -1.0]).is_err());
}

#[test]
fn contrast_all_zero_numeric_rejects_overflowed_design_score() {
    let counts = CountMatrix::from_row_major_u32(1, 1, vec![0]).unwrap();
    let design = DesignMatrix::from_row_major(1, 2, vec![f64::MAX, f64::MAX], None).unwrap();

    let err = contrast_all_zero_numeric(&counts, &design, &[1.0, -1.0]).unwrap_err();

    assert!(err.to_string().contains("contrastAllZero design score"));
}

#[test]
fn contrast_all_zero_factor_levels_matches_deseq2_character_shape() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        6,
        vec![
            0, 0, 0, 0, 50, 60, //
            10, 12, 0, 0, 50, 60, //
            0, 0, 0, 0, 0, 0, //
            0, 0, 0, 4, 50, 60,
        ],
    )
    .unwrap();
    let levels = vec!["A", "A", "B", "B", "C", "C"];

    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "B", "A").unwrap(),
        vec![true, false, true, false]
    );
    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "C", "A").unwrap(),
        vec![false, false, true, false]
    );
}

#[test]
fn original_zero_zero_d_vs_b_contrast_shape_is_preserved() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            100, 110, 0, 0, 100, 110, 0, 0, //
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
    )
    .unwrap();
    let levels = vec!["A", "A", "B", "B", "C", "C", "D", "D"];
    let design = DesignMatrix::from_row_major(
        8,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, 1.0, //
            1.0, 0.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "condition_C_vs_A".into(),
            "condition_D_vs_A".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "D", "B").unwrap(),
        vec![true, true]
    );
    assert_eq!(
        contrast_all_zero_numeric(&counts, &design, &[0.0, -1.0, 0.0, 1.0]).unwrap(),
        vec![true, true]
    );
    assert_eq!(
        contrast_all_zero_numeric(&counts, &design, &[0.0, 0.0, 0.0, 1.0]).unwrap(),
        vec![false, false]
    );
}

#[test]
fn contrast_all_zero_factor_levels_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![0, 0]).unwrap();
    assert!(contrast_all_zero_factor_levels(&counts, &["A"], "B", "A").is_err());
    assert!(contrast_all_zero_factor_levels(&counts, &["A", "B"], "A", "A").is_err());
    assert!(contrast_all_zero_factor_levels(&counts, &["A", "B"], "", "A").is_err());
}
