use super::*;

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
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level_with_reference("condition", "B", "A", "C")
        )
        .unwrap(),
        vec![-1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_resolves_formula_factor_transform_names() {
    let design = DesignMatrix::from_row_major(
        4,
        11,
        vec![1.0; 44],
        Some(vec![
            "Intercept".into(),
            "factor(condition, levels = c(\"B\", \"A\"))_A_vs_B".into(),
            "ordered(condition, levels = c(\"B\", \"A\"))_A_vs_B".into(),
            "relevel(ordered(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_C_vs_B"
                .into(),
            "droplevels(condition)_A_vs_B".into(),
            "relevel(droplevels(condition), ref = \"A\")_B_vs_A".into(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))_treated_vs_control".into(),
            "ordered(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))_treated_vs_control".into(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control group\", \"treated-group\"))_treated-group_vs_control group".into(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control factor(A)\", \"treated factor(B)\"))_treated factor(B)_vs_control factor(A)".into(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control, baseline\", \"treated, response\"))_treated, response_vs_control, baseline".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("factor.condition..levels...c..B....A...", "A", "B")
        )
        .unwrap(),
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("ordered.condition..levels...c..B....A...", "A", "B")
        )
        .unwrap(),
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "relevel.ordered.condition..levels...c..C....B....A.....ref....B..",
                "C",
                "B",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("droplevels.condition.", "A", "B")
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("relevel.droplevels.condition...ref....A..", "B", "A")
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "factor.condition..levels...c..A....B....labels...c..control....treated...",
                "treated",
                "control",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "ordered.condition..levels...c..A....B....labels...c..control....treated...",
                "treated",
                "control",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "factor.condition..levels...c..A....B....labels...c..control.group....treated.group...",
                "treated-group",
                "control group",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "factor.condition..levels...c..A....B....labels...c..control.factor.A.....treated.factor.B....",
                "treated factor(B)",
                "control factor(A)",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "factor.condition..levels...c..A....B....labels...c..control..baseline....treated..response...",
                "treated, response",
                "control, baseline",
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_resolves_no_intercept_formula_factor_labels() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            0.0, 1.0, //
            0.0, 1.0,
        ],
        Some(vec![
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control group\", \"treated-group\"))control group".into(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control group\", \"treated-group\"))treated-group".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level(
                "factor.condition..levels...c..A....B....labels...c..control.group....treated.group...",
                "treated-group",
                "control group",
            )
        )
        .unwrap(),
        vec![-1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_infers_shared_reference_with_r_like_names() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "condition_B.1_vs_A.1".into(),
            "condition_C.1_vs_A.1".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("condition", "C-1", "B-1")
        )
        .unwrap(),
        vec![0.0, -1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_infers_shared_reference_without_factor_level_separator() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "conditionB_vs_A".into(),
            "conditionC_vs_A".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "C", "B")).unwrap(),
        vec![0.0, -1.0, 1.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "B", "A")).unwrap(),
        vec![0.0, 1.0, 0.0]
    );
}

#[test]
fn factor_level_contrast_rejects_ambiguous_shared_reference_names() {
    let design = DesignMatrix::from_row_major(
        4,
        4,
        vec![1.0; 16],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "conditionB_vs_A".into(),
            "condition_C_vs_A".into(),
        ]),
    )
    .unwrap();

    let err =
        resolve_contrast(&design, &ContrastSpec::factor_level("condition", "C", "B")).unwrap_err();
    assert!(err.to_string().contains("resolves ambiguously"));
}

#[test]
fn factor_level_contrast_rejects_ambiguous_candidate_name_aliases() {
    let treatment_design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "condition_B.1_vs_A.1".into(),
            "conditionB.1_vs_A.1".into(),
        ]),
    )
    .unwrap();

    let treatment_err = resolve_contrast(
        &treatment_design,
        &ContrastSpec::factor_level("condition", "B-1", "A-1"),
    )
    .unwrap_err();
    assert!(treatment_err.to_string().contains("ambiguously"));

    let expanded_design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "conditionA.1".into(),
            "conditionB.1".into(),
            "conditionB.1".into(),
        ]),
    )
    .unwrap();

    let expanded_err = resolve_contrast(
        &expanded_design,
        &ContrastSpec::factor_level_with_reference("condition", "B-1", "A-1", "A-1"),
    )
    .unwrap_err();
    assert!(expanded_err.to_string().contains("ambiguously"));
}

#[test]
fn factor_level_contrast_rejects_ambiguous_formula_factor_transform_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "factor(condition, levels = c(\"B\", \"A\"))_A_vs_B".into(),
            "factor.condition..levels...c..B....A..._A_vs_B".into(),
        ]),
    )
    .unwrap();

    let err = resolve_contrast(
        &design,
        &ContrastSpec::factor_level("factor.condition..levels...c..B....A...", "A", "B"),
    )
    .unwrap_err();
    assert!(err.to_string().contains("ambiguously"));
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
fn factor_level_contrast_uses_r_reserved_word_make_names_for_candidates() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "condition_if._vs_TRUE.".into(),
            "condition_function._vs_TRUE.".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("condition", "if", "TRUE")
        )
        .unwrap(),
        vec![0.0, 1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::factor_level("condition", "function", "if")
        )
        .unwrap(),
        vec![0.0, -1.0, 1.0]
    );
}

#[test]
fn factor_level_contrast_validates_inputs() {
    let design = named_design();
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "A", "A")).is_err());
    assert!(resolve_contrast(&design, &ContrastSpec::factor_level("condition", "C", "B")).is_err());
}
