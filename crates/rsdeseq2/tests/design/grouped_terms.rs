use super::*;

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
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
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
    assert!(
        removed_higher
            .standard_design
            .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose:time")
            .is_err()
    );
    assert!(
        !removed_higher
            .higher_order_interactions
            .iter()
            .any(|interaction| interaction == "condition:batch:dose:time")
    );

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
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
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

    let grouped_nested_in =
        expanded_formula_design("~ (batch + dose) %in% condition", &factors, &numeric).unwrap();
    let expanded_nested_in =
        expanded_formula_design("~ condition:batch + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(grouped_nested_in, expanded_nested_in);

    let signed_grouped_star =
        expanded_formula_design("~ condition * (batch + dose - dose)", &factors, &numeric).unwrap();
    let expanded_signed_star =
        expanded_formula_design("~ condition + batch + condition:batch", &factors, &numeric)
            .unwrap();
    assert_eq!(signed_grouped_star, expanded_signed_star);

    let signed_grouped_interaction =
        expanded_formula_design("~ condition:(batch + dose - dose)", &factors, &numeric).unwrap();
    let expanded_signed_interaction =
        expanded_formula_design("~ condition:batch", &factors, &numeric).unwrap();
    assert_eq!(signed_grouped_interaction, expanded_signed_interaction);

    let signed_grouped_nested =
        expanded_formula_design("~ condition / (batch + dose - dose)", &factors, &numeric).unwrap();
    let expanded_signed_nested =
        expanded_formula_design("~ condition + condition:batch", &factors, &numeric).unwrap();
    assert_eq!(signed_grouped_nested, expanded_signed_nested);

    let signed_grouped_nested_in = expanded_formula_design(
        "~ (condition + batch - batch) %in% dose",
        &factors,
        &numeric,
    )
    .unwrap();
    let expanded_signed_nested_in =
        expanded_formula_design("~ dose:condition", &factors, &numeric).unwrap();
    assert_eq!(signed_grouped_nested_in, expanded_signed_nested_in);

    let removed_group = expanded_formula_design(
        "~ condition + batch + dose - (batch + dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_group, expected_removed);

    let signed_group =
        expanded_formula_design("~ condition + (batch + dose - dose)", &factors, &numeric).unwrap();
    let expected_signed_group =
        expanded_formula_design("~ condition + batch", &factors, &numeric).unwrap();
    assert_eq!(signed_group, expected_signed_group);

    let inverted_signed_group =
        expanded_formula_design("~ condition + batch - (batch - dose)", &factors, &numeric)
            .unwrap();
    let expected_inverted_signed_group =
        expanded_formula_design("~ condition + dose", &factors, &numeric).unwrap();
    assert_eq!(inverted_signed_group, expected_inverted_signed_group);

    let removed_nested_in = expanded_formula_design(
        "~ condition + condition:batch + condition:dose - (batch + dose) %in% condition",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed_nested_in =
        expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_nested_in, expected_removed_nested_in);
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
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
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

    let nested_signed_group = expanded_formula_design(
        "~ condition + batch + (dose + (time - dose))",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_nested_signed_group =
        expanded_formula_design("~ condition + batch + time", &factors, &numeric).unwrap();
    assert_eq!(nested_signed_group, expected_nested_signed_group);
}

#[test]
fn expanded_formula_design_supports_parenthesized_power_terms() {
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
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
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

    let pairwise =
        expanded_formula_design("~ (condition + batch + dose)^2", &factors, &numeric).unwrap();
    let explicit_pairwise = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(pairwise, explicit_pairwise);

    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: condition.clone(),
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: batch.clone(),
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose".to_string(),
                values: dose.to_vec(),
            },
            FormulaNumericColumn {
                name: "time".to_string(),
                values: time.to_vec(),
            },
        ],
    };
    let model_frame_pairwise =
        expanded_formula_design_from_model_frame("~ (condition + batch + dose)^2", &model_frame)
            .unwrap();
    assert_eq!(model_frame_pairwise, explicit_pairwise);

    let higher_order =
        expanded_formula_design("~ (condition + batch + dose + time)^3", &factors, &numeric)
            .unwrap();
    let explicit_higher_order = expanded_formula_design(
        "~ condition + batch + dose + time + condition:batch + condition:dose + condition:time + batch:dose + batch:time + dose:time + condition:batch:dose + condition:batch:time + condition:dose:time + batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(higher_order, explicit_higher_order);

    let first_order =
        expanded_formula_design("~ (condition + batch + dose)^1", &factors, &numeric).unwrap();
    let explicit_first_order =
        expanded_formula_design("~ condition + batch + dose", &factors, &numeric).unwrap();
    assert_eq!(first_order, explicit_first_order);

    let signed_base =
        expanded_formula_design("~ (condition + batch + dose - dose)^2", &factors, &numeric)
            .unwrap();
    let explicit_signed_base =
        expanded_formula_design("~ condition + batch + condition:batch", &factors, &numeric)
            .unwrap();
    assert_eq!(signed_base, explicit_signed_base);

    let signed_nested_base = expanded_formula_design(
        "~ (condition + (batch + dose) - dose)^2",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(signed_nested_base, explicit_signed_base);

    let signed_absent_removal =
        expanded_formula_design("~ (condition + batch - missing)^2", &factors, &numeric).unwrap();
    assert_eq!(signed_absent_removal, explicit_signed_base);

    assert!(expanded_formula_design("~ (condition + batch)^0", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ (condition + batch)^2.5", &factors, &numeric).is_err());
}

#[test]
fn expanded_formula_design_supports_dot_power_terms() {
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
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let dot_pairwise = expanded_formula_design("~ .^2", &factors, &numeric).unwrap();
    let explicit_pairwise = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(dot_pairwise, explicit_pairwise);

    let dot_parenthesized = expanded_formula_design("~ (.)^2", &factors, &numeric).unwrap();
    assert_eq!(dot_parenthesized, explicit_pairwise);

    let removed_pairwise =
        expanded_formula_design("~ .^2 - .^2 + condition", &factors, &numeric).unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed_pairwise, expected_removed);

    let dot_higher = expanded_formula_design("~ .^3", &factors, &numeric).unwrap();
    let explicit_higher = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose + batch:dose + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(dot_higher, explicit_higher);

    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: condition.clone(),
                levels: None,
                reference: Some("A".to_string()),
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: batch.clone(),
                levels: None,
                reference: Some("X".to_string()),
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose".to_string(),
            values: dose.to_vec(),
        }],
    };
    let model_frame_dot = expanded_formula_design_from_model_frame("~ .^2", &model_frame).unwrap();
    assert_eq!(model_frame_dot, explicit_pairwise);

    assert!(expanded_formula_design("~ .^0", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ .^two", &factors, &numeric).is_err());
}
