use super::*;

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
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "C", "B")).unwrap(),
        vec![0.0, -1.0, 1.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "B", "C")).unwrap(),
        vec![0.0, 1.0, -1.0]
    );
}

#[test]
fn original_results_condition_factor_contrast_shapes_are_preserved() {
    let design = DesignMatrix::from_row_major(
        12,
        4,
        vec![1.0; 48],
        Some(vec![
            "Intercept".into(),
            "group_2_vs_1".into(),
            "condition_2_vs_1".into(),
            "condition_3_vs_1".into(),
        ]),
    )
    .unwrap();

    // Mirrors passable contrast expectations from DESeq2's test_results.R:
    // condition 1 vs 3 = -condition_3_vs_1,
    // condition 1 vs 2 = -condition_2_vs_1,
    // condition 2 vs 3 = condition_2_vs_1 - condition_3_vs_1.
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "1", "3")).unwrap(),
        vec![0.0, 0.0, 0.0, -1.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "1", "2")).unwrap(),
        vec![0.0, 0.0, -1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "2", "3")).unwrap(),
        vec![0.0, 0.0, 1.0, -1.0]
    );

    // Original DESeq2 test_results.R also exercises list contrasts:
    // contrast=list("condition_3_vs_1", "condition_2_vs_1")
    // and listValues=c(.5, -.5).
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::list(
                vec!["condition_3_vs_1".into()],
                vec!["condition_2_vs_1".into()],
            )
        )
        .unwrap(),
        vec![0.0, 0.0, -1.0, 1.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::list_with_values(
                vec!["condition_3_vs_1".into()],
                vec!["condition_2_vs_1".into()],
                0.5,
                -0.5,
            )
        )
        .unwrap(),
        vec![0.0, 0.0, -0.5, 0.5]
    );
}

#[test]
fn original_results_invalid_contrast_shapes_are_rejected() {
    let design = DesignMatrix::from_row_major(
        12,
        4,
        vec![1.0; 48],
        Some(vec![
            "Intercept".into(),
            "group_2_vs_1".into(),
            "condition_2_vs_1".into(),
            "condition_3_vs_1".into(),
        ]),
    )
    .unwrap();

    // Passable primitive counterparts of DESeq2 test_results.R error checks:
    // missing factor/coefficient names, same numerator/denominator levels,
    // duplicated list entries, empty lists, and all-zero numeric contrasts.
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("foo", "lo", "hi")).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "4", "1")).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "1", "1")).is_err());
    assert!(
        resolve_contrast(
            &design,
            &ContrastSpec::list(vec!["condition_2_vs_1".into()], vec!["foo".into()])
        )
        .is_err()
    );
    assert!(
        resolve_contrast(
            &design,
            &ContrastSpec::list(
                vec!["condition_2_vs_1".into()],
                vec!["condition_2_vs_1".into()],
            )
        )
        .is_err()
    );
    assert!(resolve_contrast(&design, &ContrastSpec::list(Vec::new(), Vec::new())).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::Numeric(vec![0.0, 0.0, 0.0, 0.0])).is_err());
}

#[test]
fn results_contrast_character_preserves_character_all_zero_semantics() {
    let design = DesignMatrix::from_row_major(
        12,
        4,
        vec![1.0; 48],
        Some(vec![
            "Intercept".into(),
            "group_2_vs_1".into(),
            "condition_2_vs_1".into(),
            "condition_3_vs_1".into(),
        ]),
    )
    .unwrap();

    let resolved =
        resolve_results_contrast(&design, &ResultsContrast::character("condition", "1", "3"))
            .unwrap();

    assert_eq!(resolved.numeric, vec![0.0, 0.0, 0.0, -1.0]);
    assert_eq!(resolved.result_name, "condition_X1_vs_X3");
    assert_eq!(
        resolved.comparison,
        "factor-level contrast: condition 1 vs 3"
    );
    assert_eq!(
        resolved.all_zero,
        ResultsContrastAllZero::Character {
            factor: "condition".into(),
            numerator: "1".into(),
            denominator: "3".into(),
        }
    );

    let formula_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec![
            "Intercept".into(),
            "factor(condition, levels = c(\"B\", \"A\"))_A_vs_B".into(),
        ]),
    )
    .unwrap();
    let formula_resolved = resolve_results_contrast(
        &formula_design,
        &ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B"),
    )
    .unwrap();
    assert_eq!(formula_resolved.numeric, vec![0.0, 1.0]);
    assert_eq!(
        formula_resolved.result_name,
        "factor.condition..levels...c..B....A..._A_vs_B"
    );
    assert_eq!(
        formula_resolved.comparison,
        "factor-level contrast: factor.condition..levels...c..B....A... A vs B"
    );
    assert_eq!(
        formula_resolved.all_zero,
        ResultsContrastAllZero::Character {
            factor: "factor.condition..levels...c..B....A...".into(),
            numerator: "A".into(),
            denominator: "B".into(),
        }
    );
}

#[test]
fn results_contrast_list_and_numeric_use_numeric_all_zero_semantics() {
    let design = named_design();

    let list = resolve_results_contrast(
        &design,
        &ResultsContrast::list_with_values(
            vec!["condition_B_vs_A".into()],
            vec!["batch_Y_vs_X".into()],
            0.5,
            -0.5,
        ),
    )
    .unwrap();
    assert_eq!(list.numeric, vec![0.0, 0.5, -0.5]);
    assert_eq!(list.result_name, "contrast");
    assert_eq!(
        list.comparison,
        "coefficient list contrast: 0.5 condition_B_vs_A vs 0.5 batch_Y_vs_X"
    );
    assert_eq!(list.all_zero, ResultsContrastAllZero::Numeric);

    let numeric =
        resolve_results_contrast(&design, &ResultsContrast::numeric(vec![0.0, 1.0, -1.0])).unwrap();
    assert_eq!(numeric.numeric, vec![0.0, 1.0, -1.0]);
    assert_eq!(numeric.result_name, "contrast");
    assert_eq!(numeric.comparison, "primitive numeric contrast");
    assert_eq!(numeric.all_zero, ResultsContrastAllZero::Numeric);
}
