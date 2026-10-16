use super::*;

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

    let interaction =
        expanded_formula_design("~ relevel(condition, ref=\"B\"):batch", &factors, &[]).unwrap();
    assert_eq!(
        interaction.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(condition, ref = \"B\")B:batchX".to_string(),
            "relevel(condition, ref = \"B\")B:batchY".to_string(),
            "relevel(condition, ref = \"B\")A:batchX".to_string(),
            "relevel(condition, ref = \"B\")A:batchY".to_string(),
            "relevel(condition, ref = \"B\")C:batchX".to_string(),
            "relevel(condition, ref = \"B\")C:batchY".to_string(),
        ]
    );

    let nested =
        expanded_formula_design("~ relevel(factor(condition), ref=\"B\")", &factors, &[]).unwrap();
    assert_eq!(
        nested.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(factor(condition), ref = \"B\")_A_vs_B".to_string(),
            "relevel(factor(condition), ref = \"B\")_C_vs_B".to_string(),
        ]
    );

    let cell_levels = vec![
        "T cell".to_string(),
        "B-cell".to_string(),
        "NK cell".to_string(),
    ];
    let cell_type = vec![
        "T cell".to_string(),
        "B-cell".to_string(),
        "NK cell".to_string(),
        "B-cell".to_string(),
    ];
    let cell_factors = [ExpandedFactorSpec {
        factor: "cell type",
        sample_levels: &cell_type,
        reference: "T cell",
        levels: Some(&cell_levels),
    }];
    let cleaned_ref =
        expanded_formula_design("~ relevel(`cell type`, ref=\"B.cell\")", &cell_factors, &[])
            .unwrap();
    assert_eq!(
        cleaned_ref.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(cell type, ref = \"B-cell\")_T cell_vs_B-cell".to_string(),
            "relevel(cell type, ref = \"B-cell\")_NK cell_vs_B-cell".to_string(),
        ]
    );

    let nested_named = expanded_formula_design(
        "~ relevel(x=as.factor(condition), ref=\"C\")",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        nested_named.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(as.factor(condition), ref = \"C\")_A_vs_C".to_string(),
            "relevel(as.factor(condition), ref = \"C\")_B_vs_C".to_string(),
        ]
    );

    let nested_ordered_named = expanded_formula_design(
        "~ relevel(x=as.ordered(condition), ref=\"C\")",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        nested_ordered_named
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(as.ordered(condition), ref = \"C\")_A_vs_C".to_string(),
            "relevel(as.ordered(condition), ref = \"C\")_B_vs_C".to_string(),
        ]
    );

    let nested_as_factor =
        expanded_formula_design("~ factor(relevel(condition, ref=\"B\"))", &factors, &[]).unwrap();
    assert_eq!(
        nested_as_factor
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "factor(relevel(condition, ref = \"B\"))_A_vs_B".to_string(),
            "factor(relevel(condition, ref = \"B\"))_C_vs_B".to_string(),
        ]
    );

    let nested_ordered_relevel =
        expanded_formula_design("~ ordered(relevel(condition, ref=\"B\"))", &factors, &[]).unwrap();
    assert_eq!(
        nested_ordered_relevel
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(relevel(condition, ref = \"B\"))_A_vs_B".to_string(),
            "ordered(relevel(condition, ref = \"B\"))_C_vs_B".to_string(),
        ]
    );

    let nested_levels = expanded_formula_design(
        "~ relevel(factor(condition, levels=c('C','B','A')), ref='B')",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        nested_levels.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_C_vs_B"
                .to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_A_vs_B"
                .to_string(),
        ]
    );

    let nested_ordered_levels = expanded_formula_design(
        "~ relevel(ordered(condition, levels=c('C','B','A')), ref='B')",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        nested_ordered_levels
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(ordered(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_C_vs_B"
                .to_string(),
            "relevel(ordered(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_A_vs_B"
                .to_string(),
        ]
    );

    let spaced_nested = expanded_formula_design(
        "~ relevel ( factor ( condition, levels = c('C', 'B', 'A') ), ref = 'B' )",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        spaced_nested.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_C_vs_B"
                .to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"B\")_A_vs_B"
                .to_string(),
        ]
    );

    let nested_named_levels = expanded_formula_design(
        "~ relevel(x=factor(levels=c('C','B','A'), x=condition), ref='C')",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        nested_named_levels
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"C\")_B_vs_C"
                .to_string(),
            "relevel(factor(condition, levels = c(\"C\", \"B\", \"A\")), ref = \"C\")_A_vs_C"
                .to_string(),
        ]
    );

    let equals_condition = vec![
        "A=B".to_string(),
        "C".to_string(),
        "A=B".to_string(),
        "C".to_string(),
    ];
    let equals_levels = vec!["A=B".to_string(), "C".to_string()];
    let equals_factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &equals_condition,
        reference: "A=B",
        levels: Some(&equals_levels),
    }];
    let positional_equals_reference = expanded_formula_design(
        "~ relevel(factor(condition, levels=c('A=B','C')), 'A=B')",
        &equals_factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        positional_equals_reference
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(factor(condition, levels = c(\"A=B\", \"C\")), ref = \"A=B\")_C_vs_A=B"
                .to_string(),
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

    let ordered =
        expanded_formula_design("~ ordered(condition) + as.ordered(x=batch)", &factors, &[])
            .unwrap();
    assert_eq!(
        ordered.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition)_B_vs_A".to_string(),
            "as.ordered(batch)_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        ordered.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition)A".to_string(),
            "ordered(condition)B".to_string(),
            "as.ordered(batch)X".to_string(),
            "as.ordered(batch)Y".to_string(),
        ]
    );

    let interaction =
        expanded_formula_design("~ factor(condition):as.factor(batch)", &factors, &[]).unwrap();
    assert_eq!(
        interaction.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)A:as.factor(batch)X".to_string(),
            "factor(condition)A:as.factor(batch)Y".to_string(),
            "factor(condition)B:as.factor(batch)X".to_string(),
            "factor(condition)B:as.factor(batch)Y".to_string(),
        ]
    );
    assert_eq!(
        interaction.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)A:as.factor(batch)X".to_string(),
            "factor(condition)A:as.factor(batch)Y".to_string(),
            "factor(condition)B:as.factor(batch)X".to_string(),
            "factor(condition)B:as.factor(batch)Y".to_string(),
        ]
    );

    let spaced = expanded_formula_design(
        "~ factor ( condition ) + as.factor ( x = batch )",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        spaced.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)_B_vs_A".to_string(),
            "as.factor(batch)_Y_vs_X".to_string(),
        ]
    );

    let reordered =
        expanded_formula_design("~ factor(condition, levels=c('B','A'))", &factors, &[]).unwrap();
    assert_eq!(
        reordered.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"B\", \"A\"))_A_vs_B".to_string(),
        ]
    );
    assert_eq!(
        reordered.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"B\", \"A\"))B".to_string(),
            "factor(condition, levels = c(\"B\", \"A\"))A".to_string(),
        ]
    );

    let reordered_named = expanded_formula_design(
        "~ factor(levels=c(\"B\", \"A\"), x=condition)",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(reordered_named, reordered);

    let relabeled = expanded_formula_design(
        "~ factor(condition, levels=c('A','B'), labels=c('control','treated'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))_treated_vs_control".to_string(),
        ]
    );
    assert_eq!(
        relabeled.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))control".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))treated".to_string(),
        ]
    );

    let relabeled_positional = expanded_formula_design(
        "~ factor(condition, c('A','B'), c('control','treated'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(relabeled_positional, relabeled);

    let ordered_reordered =
        expanded_formula_design("~ ordered(condition, levels=c('B','A'))", &factors, &[]).unwrap();
    assert_eq!(
        ordered_reordered
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition, levels = c(\"B\", \"A\"))_A_vs_B".to_string(),
        ]
    );
    assert_eq!(
        ordered_reordered
            .expanded_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition, levels = c(\"B\", \"A\"))B".to_string(),
            "ordered(condition, levels = c(\"B\", \"A\"))A".to_string(),
        ]
    );

    let ordered_relabeled = expanded_formula_design(
        "~ ordered(condition, levels=c('A','B'), labels=c('control','treated'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        ordered_relabeled
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))_treated_vs_control".to_string(),
        ]
    );
    assert_eq!(
        ordered_relabeled
            .expanded_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "ordered(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))control".to_string(),
            "ordered(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))treated".to_string(),
        ]
    );
    let ordered_relabeled_named = expanded_formula_design(
        "~ ordered(x=condition, labels=c('control','treated'), levels=c('A','B'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(ordered_relabeled_named, ordered_relabeled);

    let relabeled_non_syntactic = expanded_formula_design(
        "~ factor(condition, levels=c('A','B'), labels=c('control group','treated-group'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled_non_syntactic
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control group\", \"treated-group\"))_treated-group_vs_control group".to_string(),
        ]
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: condition.clone(),
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

    let relabeled_with_offsets = expanded_formula_design_with_offsets(
        "~ factor(condition, levels=c('A','B'), labels=c('control','treated'))",
        &factors,
        &[],
    )
    .unwrap();
    let relabeled_factor = relabeled_with_offsets
        .model_frame
        .factors
        .iter()
        .find(|factor| factor.name.starts_with("factor(condition, levels"))
        .unwrap();
    assert_eq!(
        relabeled_factor.sample_levels,
        vec![
            "control".to_string(),
            "treated".to_string(),
            "control".to_string(),
            "treated".to_string()
        ]
    );
    assert_eq!(
        relabeled_factor.levels.as_deref(),
        Some(&["control".to_string(), "treated".to_string()][..])
    );
    assert_eq!(relabeled_factor.reference.as_deref(), Some("control"));

    let relabeled_non_syntactic_with_offsets = expanded_formula_design_with_offsets(
        "~ factor(condition, levels=c('A','B'), labels=c('control group','treated-group'))",
        &factors,
        &[],
    )
    .unwrap();
    let relabeled_non_syntactic_factor = relabeled_non_syntactic_with_offsets
        .model_frame
        .factors
        .iter()
        .find(|factor| factor.name.contains("control group"))
        .unwrap();
    assert_eq!(
        relabeled_non_syntactic_factor.levels.as_deref(),
        Some(&["control group".to_string(), "treated-group".to_string()][..])
    );
    assert_eq!(
        relabeled_non_syntactic_factor.reference.as_deref(),
        Some("control group")
    );

    let relabeled_transform_text = expanded_formula_design(
        "~ factor(condition, levels=c('A','B'), labels=c('control factor(A)','treated factor(B)'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled_transform_text
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control factor(A)\", \"treated factor(B)\"))_treated factor(B)_vs_control factor(A)".to_string(),
        ]
    );
    let relabeled_transform_text_with_offsets = expanded_formula_design_with_offsets(
        "~ factor(condition, levels=c('A','B'), labels=c('control factor(A)','treated factor(B)'))",
        &factors,
        &[],
    )
    .unwrap();
    let relabeled_transform_text_factor = relabeled_transform_text_with_offsets
        .model_frame
        .factors
        .iter()
        .find(|factor| factor.name.contains("control factor(A)"))
        .unwrap();
    assert_eq!(
        relabeled_transform_text_factor.levels.as_deref(),
        Some(
            &[
                "control factor(A)".to_string(),
                "treated factor(B)".to_string()
            ][..]
        )
    );
    assert_eq!(
        relabeled_transform_text_factor.reference.as_deref(),
        Some("control factor(A)")
    );

    let relabeled_comma_text = expanded_formula_design(
        "~ factor(condition, levels=c('A','B'), labels=c('control, baseline','treated, response'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled_comma_text
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control, baseline\", \"treated, response\"))_treated, response_vs_control, baseline".to_string(),
        ]
    );
    let relabeled_comma_text_with_offsets = expanded_formula_design_with_offsets(
        "~ factor(condition, levels=c('A','B'), labels=c('control, baseline','treated, response'))",
        &factors,
        &[],
    )
    .unwrap();
    let relabeled_comma_text_factor = relabeled_comma_text_with_offsets
        .model_frame
        .factors
        .iter()
        .find(|factor| factor.name.contains("control, baseline"))
        .unwrap();
    assert_eq!(
        relabeled_comma_text_factor.levels.as_deref(),
        Some(
            &[
                "control, baseline".to_string(),
                "treated, response".to_string()
            ][..]
        )
    );
    assert_eq!(
        relabeled_comma_text_factor.reference.as_deref(),
        Some("control, baseline")
    );

    let relabeled_escaped_text = expanded_formula_design(
        r#"~ factor(condition, levels=c('A','B'), labels=c("control \"quote\"","treated \\ path"))"#,
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled_escaped_text
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            r#"factor(condition, levels = c("A", "B"), labels = c("control \"quote\"", "treated \\ path"))_treated \ path_vs_control "quote""#.to_string(),
        ]
    );
    let relabeled_escaped_text_with_offsets = expanded_formula_design_with_offsets(
        r#"~ factor(condition, levels=c('A','B'), labels=c("control \"quote\"","treated \\ path"))"#,
        &factors,
        &[],
    )
    .unwrap();
    let relabeled_escaped_text_factor = relabeled_escaped_text_with_offsets
        .model_frame
        .factors
        .iter()
        .find(|factor| factor.name.contains(r#"control \"quote\""#))
        .unwrap();
    assert_eq!(
        relabeled_escaped_text_factor.levels.as_deref(),
        Some(
            &[
                "control \"quote\"".to_string(),
                "treated \\ path".to_string()
            ][..]
        )
    );
    assert_eq!(
        relabeled_escaped_text_factor.reference.as_deref(),
        Some("control \"quote\"")
    );
}

#[test]
fn expanded_formula_design_supports_droplevels_factor_transform() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "B".to_string(),
                "A".to_string(),
                "B".to_string(),
                "A".to_string(),
            ],
            levels: Some(vec!["Z".to_string(), "B".to_string(), "A".to_string()]),
            reference: Some("Z".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };

    let unchanged = expanded_formula_design_from_model_frame("~ condition", &model_frame).unwrap();
    assert_eq!(
        unchanged.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_Z".to_string(),
            "condition_A_vs_Z".to_string(),
        ]
    );

    let dropped =
        expanded_formula_design_from_model_frame("~ droplevels(condition)", &model_frame).unwrap();
    assert_eq!(
        dropped.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(condition)_A_vs_B".to_string(),
        ]
    );
    assert_eq!(
        dropped.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(condition)B".to_string(),
            "droplevels(condition)A".to_string(),
        ]
    );

    let relevel_dropped = expanded_formula_design_from_model_frame(
        "~ relevel(droplevels(condition), ref='A')",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        relevel_dropped.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "relevel(droplevels(condition), ref = \"A\")_B_vs_A".to_string(),
        ]
    );

    let dropped_reordered = expanded_formula_design_from_model_frame(
        "~ droplevels(factor(condition, levels=c('A','B','Z')))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        dropped_reordered
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(factor(condition, levels = c(\"A\", \"B\", \"Z\")))_B_vs_A".to_string(),
        ]
    );
    assert_eq!(
        dropped_reordered
            .expanded_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(factor(condition, levels = c(\"A\", \"B\", \"Z\")))A".to_string(),
            "droplevels(factor(condition, levels = c(\"A\", \"B\", \"Z\")))B".to_string(),
        ]
    );

    let dropped_ordered = expanded_formula_design_from_model_frame(
        "~ droplevels(ordered(condition, levels=c('A','B','Z')))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        dropped_ordered.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(ordered(condition, levels = c(\"A\", \"B\", \"Z\")))_B_vs_A".to_string(),
        ]
    );
    assert_eq!(
        dropped_ordered.expanded_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(ordered(condition, levels = c(\"A\", \"B\", \"Z\")))A".to_string(),
            "droplevels(ordered(condition, levels = c(\"A\", \"B\", \"Z\")))B".to_string(),
        ]
    );

    let dropped_releveled = expanded_formula_design_from_model_frame(
        "~ droplevels(relevel(condition, ref='A'))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        dropped_releveled
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(relevel(condition, ref = \"A\"))_B_vs_A".to_string(),
        ]
    );
    let dropped_releveled_with_offsets = expanded_formula_design_with_offsets_from_model_frame(
        "~ droplevels(relevel(condition, ref='A'))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(
        dropped_releveled_with_offsets
            .model_frame
            .factors
            .iter()
            .find(|factor| factor.name == "droplevels(relevel(condition, ref = \"A\"))")
            .and_then(|factor| factor.reference.as_deref()),
        Some("A")
    );

    let observed_reference = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec![
                "B".to_string(),
                "A".to_string(),
                "B".to_string(),
                "A".to_string(),
            ],
            levels: Some(vec!["Z".to_string(), "B".to_string(), "A".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let dropped_observed_reference =
        expanded_formula_design_from_model_frame("~ droplevels(condition)", &observed_reference)
            .unwrap();
    assert_eq!(
        dropped_observed_reference
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "droplevels(condition)_B_vs_A".to_string(),
        ]
    );
}
