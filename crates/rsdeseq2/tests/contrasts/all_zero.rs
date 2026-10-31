use super::*;

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
fn contrast_all_zero_numeric_keeps_large_cancelling_design_score_finite() {
    let counts = CountMatrix::from_row_major_u32(1, 1, vec![0]).unwrap();
    let design = DesignMatrix::from_row_major(1, 2, vec![f64::MAX, -f64::MAX], None).unwrap();

    let flags = contrast_all_zero_numeric(&counts, &design, &[1.0, -1.0]).unwrap();

    assert_eq!(flags, vec![true]);
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
fn contrast_all_zero_factor_levels_resolves_cleaned_level_aliases() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            0, 0, 5, 6, //
            0, 3, 0, 0,
        ],
    )
    .unwrap();
    let levels = vec!["A cell", "A cell", "B-cell", "B-cell"];

    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "B.cell", "A.cell").unwrap(),
        vec![false, false]
    );
    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "B-cell", "A.cell").unwrap(),
        vec![false, false]
    );

    let same_level = contrast_all_zero_factor_levels(&counts, &levels, "A.cell", "A cell")
        .unwrap_err()
        .to_string();
    assert!(same_level.contains("resolve to the same level"));
}

#[test]
fn contrast_all_zero_factor_levels_rejects_ambiguous_cleaned_level_aliases() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![0, 0, 0, 0]).unwrap();
    let levels = vec!["A cell", "A-cell", "B", "B"];

    let err = contrast_all_zero_factor_levels(&counts, &levels, "A.cell", "B")
        .unwrap_err()
        .to_string();
    assert!(err.contains("resolves ambiguously"));
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
        resolve_contrast(
            &design,
            &ContrastSpec::list(
                vec!["condition_D_vs_A".into()],
                vec!["condition_B_vs_A".into()],
            )
        )
        .unwrap(),
        vec![0.0, -1.0, 0.0, 1.0]
    );
    assert_eq!(
        contrast_all_zero_factor_levels(&counts, &levels, "D", "A").unwrap(),
        vec![false, true]
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
    assert!(contrast_all_zero_factor_levels(&counts, &["A", "A"], "B", "A").is_err());
    assert!(contrast_all_zero_factor_levels(&counts, &["B", "B"], "B", "A").is_err());
}
