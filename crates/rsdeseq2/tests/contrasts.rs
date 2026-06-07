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
        2,
        vec![1.0; 8],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A:poly(dose, 2)1".into(),
        ]),
    )
    .unwrap();

    assert_eq!(
        resolve_contrast(
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A:poly.dose..2.1")
        )
        .unwrap(),
        vec![0.0, 1.0]
    );
    assert_eq!(
        resolve_coefficient_index(&design, "condition_B_vs_A:poly.dose..2.1").unwrap(),
        1
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
    assert!(resolve_contrast(
        &design,
        &ContrastSpec::list(vec!["condition_2_vs_1".into()], vec!["foo".into()])
    )
    .is_err());
    assert!(resolve_contrast(
        &design,
        &ContrastSpec::list(
            vec!["condition_2_vs_1".into()],
            vec!["condition_2_vs_1".into()],
        )
    )
    .is_err());
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
            numerator: "1".into(),
            denominator: "3".into(),
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

#[test]
fn model_frame_factor_level_contrast_infers_and_overrides_reference() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "C".to_string(),
            ],
            levels: None,
            reference: Some("B".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0],
        }],
    };

    let inferred_request = ResultsContrast::character("condition", "C", "A");
    let inferred = factor_level_contrast_from_model_frame(&inferred_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(inferred.factor, "condition");
    assert_eq!(inferred.numerator, "C");
    assert_eq!(inferred.denominator, "A");
    assert_eq!(inferred.reference, Some("B"));
    assert_eq!(
        inferred.sample_levels,
        &model_frame.factors[0].sample_levels
    );

    let explicit_request = ResultsContrast::character_with_reference("condition", "C", "B", "A");
    let explicit = factor_level_contrast_from_model_frame(&explicit_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(explicit.reference, Some("A"));

    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::numeric(vec![0.0, 1.0]),
        &model_frame
    )
    .unwrap()
    .is_none());
    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::list(vec!["condition_C_vs_A".into()], Vec::new()),
        &model_frame,
    )
    .unwrap()
    .is_none());
}

#[test]
fn model_frame_factor_level_contrast_uses_declared_level_reference() {
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

    let inferred_request = ResultsContrast::character("condition", "B", "A");
    let inferred = factor_level_contrast_from_model_frame(&inferred_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(inferred.reference, Some("A"));

    let explicit_request = ResultsContrast::character_with_reference("condition", "C", "B", "B");
    let explicit = factor_level_contrast_from_model_frame(&explicit_request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(explicit.reference, Some("B"));
}

#[test]
fn model_frame_factor_level_contrast_allows_unused_declared_reference() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["B".to_string(), "C".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("condition", "C", "B");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "condition");
    assert_eq!(contrast.numerator, "C");
    assert_eq!(contrast.denominator, "B");
    assert_eq!(contrast.reference, Some("A"));
    assert_eq!(
        contrast.sample_levels,
        &model_frame.factors[0].sample_levels
    );

    let missing_reference = ResultsContrast::character_with_reference("condition", "C", "B", "D");
    assert!(factor_level_contrast_from_model_frame(&missing_reference, &model_frame).is_err());
}

#[test]
fn model_frame_factor_level_contrast_accepts_r_cleaned_factor_name_alias() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec![
                "T cell".to_string(),
                "B-cell".to_string(),
                "B-cell".to_string(),
                "T cell".to_string(),
            ],
            levels: Some(vec!["T cell".to_string(), "B-cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };

    let request = ResultsContrast::character("cell.type", "B-cell", "T cell");
    let contrast = factor_level_contrast_from_model_frame(&request, &model_frame)
        .unwrap()
        .unwrap();
    assert_eq!(contrast.factor, "cell type");
    assert_eq!(contrast.reference, Some("T cell"));
    assert_eq!(
        contrast.sample_levels,
        &model_frame.factors[0].sample_levels
    );
}

#[test]
fn model_frame_factor_level_contrast_validates_factor_metadata() {
    let duplicate = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: vec!["A".to_string(), "C".to_string()],
                levels: None,
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::character("condition", "B", "A"),
        &duplicate,
    )
    .unwrap_err()
    .to_string()
    .contains("more than once"));

    let missing_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::character("condition", "C", "A"),
        &missing_level,
    )
    .unwrap_err()
    .to_string()
    .contains("numerator level"));
    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::character_with_reference("condition", "B", "A", "C"),
        &missing_level,
    )
    .unwrap_err()
    .to_string()
    .contains("reference level"));

    let ambiguous_cleaned = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "cell-type".to_string(),
                sample_levels: vec!["A".to_string(), "B".to_string()],
                levels: None,
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    assert!(factor_level_contrast_from_model_frame(
        &ResultsContrast::character("cell.type", "B", "A"),
        &ambiguous_cleaned,
    )
    .unwrap_err()
    .to_string()
    .contains("ambiguously"));
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
