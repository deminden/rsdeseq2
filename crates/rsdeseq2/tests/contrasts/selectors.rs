use super::*;

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
fn coefficient_name_contrast_resolves_r_cleaned_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "(Intercept)".into(),
            "if.".into(),
            "condition.B.1".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::coefficient_name("Intercept")).unwrap(),
        vec![1.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::coefficient_name("if")).unwrap(),
        vec![0.0, 1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::coefficient_name("condition-B 1")).unwrap(),
        vec![0.0, 0.0, 1.0]
    );
}

#[test]
fn coefficient_name_contrast_resolves_component_cleaned_factor_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("cell.type_B.cell_vs_T.cell")
        )
        .unwrap(),
        vec![0.0, 1.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "cell.type_B.cell_vs_T.cell").unwrap(),
        1
    );
}

#[test]
fn coefficient_name_contrast_resolves_component_cleaned_interaction_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell:dose value_plus_1".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("cell.type_B.cell_vs_T.cell:dose.value_plus_1")
        )
        .unwrap(),
        vec![0.0, 1.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "cell.type_B.cell_vs_T.cell:dose.value_plus_1").unwrap(),
        1
    );
}

#[test]
fn coefficient_name_contrast_resolves_cleaned_formula_transform_interaction_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        6,
        vec![1.0; 24],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A:poly(dose, 2)1".into(),
            "relevel(condition, ref = \"B\")_A_vs_B".into(),
            "factor(condition, levels = c(\"B\", \"A\"))_A_vs_B".into(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_C_vs_B"
                .into(),
            "droplevels(condition)_A_vs_B".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A:poly.dose..2.1")
        )
        .unwrap(),
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "condition_B_vs_A:poly.dose..2.1").unwrap(),
        1
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("relevel.condition..ref....B.._A_vs_B")
        )
        .unwrap(),
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "relevel.condition..ref....B.._A_vs_B").unwrap(),
        2
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("factor.condition..levels...c..B....A..._A_vs_B")
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name(
                "relevel.factor.condition..levels...c..C....B....A.....ref....B.._C_vs_B"
            )
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    );
    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("droplevels.condition._A_vs_B")
        )
        .unwrap(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    );
}

#[test]
fn coefficient_name_contrast_rejects_ambiguous_component_cleaned_interaction_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell:dose value_plus_1".into(),
            "cell-type_B cell_vs_T cell:dose value_plus_1".into(),
        ]),
    )
    .unwrap();

    let err = resolve_contrast(
        &design,
        &ContrastSpec::coefficient_name("cell.type_B.cell_vs_T.cell:dose.value_plus_1"),
    )
    .unwrap_err();
    assert!(err.to_string().contains("ambiguously"));
    assert!(
        resolve_coefficient_index(&design, "cell.type_B.cell_vs_T.cell:dose.value_plus_1")
            .unwrap_err()
            .to_string()
            .contains("ambiguously")
    );
}

#[test]
fn coefficient_name_contrast_matches_requested_cleaned_alias_against_stored_names() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec!["Intercept".into(), "dose(value)_log2".into()]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::coefficient_name("dose.value._log2")).unwrap(),
        vec![0.0, 1.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "dose.value._log2").unwrap(),
        1
    );
}

#[test]
fn coefficient_name_contrast_prefers_exact_name_over_cleaned_alias() {
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1.0; 8], Some(vec!["a-b".into(), "a.b".into()]))
            .unwrap();

    assert_eq!(
        resolve_contrast(&design, &ContrastSpec::coefficient_name("a-b")).unwrap(),
        vec![1.0, 0.0]
    );
}

#[test]
fn coefficient_name_contrast_rejects_ambiguous_r_cleaned_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0; 8],
        Some(vec![".Intercept.".into(), "Intercept".into()]),
    )
    .unwrap();

    assert!(resolve_contrast(&design, &ContrastSpec::coefficient_name("(Intercept)")).is_err());
    assert!(resolve_coefficient_index(&design, "(Intercept)").is_err());
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
    let cleaned_factor = ContrastSpec::factor_level("if", "B-1", "A 0");
    assert_eq!(cleaned_factor.result_name(), "if._B.1_vs_A.0");
    assert_eq!(
        cleaned_factor.comparison(),
        "factor-level contrast: if B-1 vs A 0"
    );
    let formula_factor =
        ContrastSpec::factor_level("factor(condition, levels = c(\"B\", \"A\"))", "A", "B");
    assert_eq!(
        formula_factor.result_name(),
        "factor.condition..levels...c..B....A..._A_vs_B"
    );
    assert_eq!(
        formula_factor.comparison(),
        "factor-level contrast: factor(condition, levels = c(\"B\", \"A\")) A vs B"
    );
    let relabeled_formula_factor = ContrastSpec::factor_level(
        "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))",
        "treated",
        "control",
    );
    assert_eq!(
        relabeled_formula_factor.result_name(),
        "factor.condition..levels...c..A....B....labels...c..control....treated..._treated_vs_control"
    );
    assert_eq!(
        relabeled_formula_factor.comparison(),
        "factor-level contrast: factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\")) treated vs control"
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
        "coefficient list contrast: 0.5 condition_B_vs_A vs 2 batch_Y_vs_X"
    );

    let positive_only = ContrastSpec::list(vec!["condition_B_vs_A".into()], Vec::new());
    assert_eq!(
        positive_only.comparison(),
        "coefficient list contrast: condition_B_vs_A effect"
    );

    let negative_only =
        ContrastSpec::list_with_values(Vec::new(), vec!["batch_Y_vs_X".into()], 0.5, -0.5);
    assert_eq!(
        negative_only.comparison(),
        "coefficient list contrast: -0.5 batch_Y_vs_X effect"
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
fn list_contrast_resolves_r_cleaned_coefficient_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "(Intercept)".into(),
            "if.".into(),
            "condition.B.1".into(),
        ]),
    )
    .unwrap();
    let contrast = ContrastSpec::list(vec!["if".into()], vec!["condition-B 1".into()]);

    assert_eq!(
        resolve_contrast(&design, &contrast).unwrap(),
        vec![0.0, 1.0, -1.0]
    );

    let duplicated = ContrastSpec::list(vec!["if".into(), "if.".into()], Vec::new());
    assert_eq!(
        resolve_contrast(&design, &duplicated).unwrap(),
        vec![0.0, 1.0, 0.0]
    );

    let overlap = ContrastSpec::list(vec!["Intercept".into()], vec!["(Intercept)".into()]);
    assert!(resolve_contrast(&design, &overlap).is_err());
}

#[test]
fn list_contrast_resolves_component_cleaned_interaction_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell:dose value_plus_1".into(),
            "batch name_Y-1_vs_X 0".into(),
        ]),
    )
    .unwrap();
    let contrast = ContrastSpec::list(
        vec!["cell.type_B.cell_vs_T.cell:dose.value_plus_1".into()],
        vec!["batch.name_Y.1_vs_X.0".into()],
    );

    assert_eq!(
        resolve_contrast(&design, &contrast).unwrap(),
        vec![0.0, 1.0, -1.0]
    );
}

#[test]
fn list_contrast_rejects_ambiguous_component_cleaned_interaction_aliases() {
    let design = DesignMatrix::from_row_major(
        4,
        3,
        vec![1.0; 12],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell:dose value_plus_1".into(),
            "cell-type_B cell_vs_T cell:dose value_plus_1".into(),
        ]),
    )
    .unwrap();
    let contrast = ContrastSpec::list(
        vec!["cell.type_B.cell_vs_T.cell:dose.value_plus_1".into()],
        Vec::new(),
    );

    let err = resolve_contrast(&design, &contrast).unwrap_err();
    assert!(err.to_string().contains("ambiguously"));
}

#[test]
fn list_contrast_requires_deseq2_style_signed_list_values() {
    let design = named_design();
    let positive_zero = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        0.0,
        -1.0,
    );
    let positive_negative = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        -1.0,
        -1.0,
    );
    let negative_zero = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        1.0,
        0.0,
    );
    let negative_positive = ContrastSpec::list_with_values(
        vec!["condition_B_vs_A".into()],
        vec!["batch_Y_vs_X".into()],
        1.0,
        1.0,
    );

    assert!(resolve_contrast(&design, &positive_zero).is_err());
    assert!(resolve_contrast(&design, &positive_negative).is_err());
    assert!(resolve_contrast(&design, &negative_zero).is_err());
    assert!(resolve_contrast(&design, &negative_positive).is_err());
}
