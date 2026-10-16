use super::*;

#[test]
fn expanded_formula_design_expands_star_terms() {
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

    let parsed = expanded_formula_design("~ condition * batch", &factors, &[]).unwrap();
    let direct = expanded_additive_design_with_all_interactions(
        &factors,
        &[],
        &[ExpandedFactorInteractionSpec {
            left_factor: "condition",
            right_factor: "batch",
        }],
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(parsed, direct);
}

#[test]
fn expanded_formula_design_accepts_interactions_without_main_effects() {
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
    let dose = [0.0, 1.0, 0.0, 1.0];
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

    let interaction_only = expanded_formula_design("~ condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        interaction_only
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA:batchX".to_string(),
            "conditionA:batchY".to_string(),
            "conditionB:batchX".to_string(),
            "conditionB:batchY".to_string(),
        ]
    );
    assert_eq!(
        interaction_only.standard_design.matrix().as_slice(),
        &[
            1.0, 1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        interaction_only.coefficient_groups,
        vec![vec![0], vec![1], vec![2], vec![3], vec![4]]
    );

    let one_main = expanded_formula_design("~ condition + condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        one_main.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        one_main.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0,
        ]
    );
    assert_eq!(
        one_main.coefficient_groups,
        vec![vec![0], vec![2], vec![4], vec![6]]
    );

    let factor_numeric =
        expanded_formula_design("~ condition + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(
        factor_numeric.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:dose".to_string(),
            "conditionB:dose".to_string(),
        ]
    );
    assert_eq!(
        factor_numeric.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, 1.0,
        ]
    );
}

#[test]
fn expanded_formula_design_supports_intercept_removal() {
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

    let no_intercept = expanded_formula_design("~ 0 + condition + batch", &factors, &[]).unwrap();
    assert_eq!(
        no_intercept.standard_design.coefficient_names().unwrap(),
        &[
            "conditionA".to_string(),
            "conditionB".to_string(),
            "batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        no_intercept.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            0.0, 1.0, 0.0, //
            0.0, 1.0, 1.0,
        ]
    );
    assert_eq!(
        no_intercept.coefficient_groups,
        vec![vec![0], vec![1], vec![3]]
    );

    let minus_one = expanded_formula_design("~ condition - 1", &factors, &[]).unwrap();
    assert_eq!(
        minus_one.standard_design.coefficient_names().unwrap(),
        &["conditionA".to_string(), "conditionB".to_string()]
    );
    assert_eq!(
        minus_one.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, //
            1.0, 0.0, //
            0.0, 1.0, //
            0.0, 1.0,
        ]
    );

    let relabeled_no_intercept = expanded_formula_design(
        "~ 0 + factor(condition, levels=c('A','B'), labels=c('control','treated'))",
        &factors,
        &[],
    )
    .unwrap();
    assert_eq!(
        relabeled_no_intercept
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))control".to_string(),
            "factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\"))treated".to_string(),
        ]
    );
    assert_eq!(
        relabeled_no_intercept.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, //
            1.0, 0.0, //
            0.0, 1.0, //
            0.0, 1.0,
        ]
    );
    assert_eq!(
        relabeled_no_intercept.coefficient_groups,
        vec![vec![0], vec![1]]
    );

    let restored = expanded_formula_design("~ 0 + 1 + condition", &factors, &[]).unwrap();
    let restored_from_minus_one =
        expanded_formula_design("~ condition - 1 + 1", &factors, &[]).unwrap();
    let with_intercept = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    assert_eq!(restored, with_intercept);
    assert_eq!(restored_from_minus_one, with_intercept);

    let plus_zero = expanded_formula_design("~ condition + 0", &factors, &[]).unwrap();
    assert_eq!(plus_zero, minus_one);

    let minus_zero = expanded_formula_design("~ condition - 0", &factors, &[]).unwrap();
    assert_eq!(minus_zero, with_intercept);

    let zero_minus_zero = expanded_formula_design("~ 0 - 0 + condition", &factors, &[]).unwrap();
    assert_eq!(zero_minus_zero, with_intercept);

    let plus_zero_plus_one = expanded_formula_design("~ condition + 0 + 1", &factors, &[]).unwrap();
    assert_eq!(plus_zero_plus_one, with_intercept);

    let removed_again = expanded_formula_design("~ 0 + 1 + condition - 1", &factors, &[]).unwrap();
    assert_eq!(removed_again, minus_one);

    let plus_zero_minus_one =
        expanded_formula_design("~ condition + 0 - 1", &factors, &[]).unwrap();
    assert_eq!(plus_zero_minus_one, minus_one);

    let minus_one_minus_zero =
        expanded_formula_design("~ condition - 1 - 0", &factors, &[]).unwrap();
    assert_eq!(minus_one_minus_zero, with_intercept);

    let minus_negative_one = expanded_formula_design("~ condition - -1", &factors, &[]).unwrap();
    assert_eq!(minus_negative_one, with_intercept);

    let zero_minus_negative_one =
        expanded_formula_design("~ 0 - -1 + condition", &factors, &[]).unwrap();
    assert_eq!(zero_minus_negative_one, with_intercept);

    let plus_zero_minus_negative_one =
        expanded_formula_design("~ condition + 0 - -1", &factors, &[]).unwrap();
    assert_eq!(plus_zero_minus_negative_one, with_intercept);

    let interaction_only = expanded_formula_design("~ 0 + condition:batch", &factors, &[]).unwrap();
    assert_eq!(
        interaction_only
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "conditionA:batchX".to_string(),
            "conditionA:batchY".to_string(),
            "conditionB:batchX".to_string(),
            "conditionB:batchY".to_string(),
        ]
    );
    assert_eq!(
        interaction_only.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, //
            0.0, 0.0, 1.0, 0.0, //
            0.0, 0.0, 0.0, 1.0,
        ]
    );
}

#[test]
fn expanded_formula_design_supports_leading_unary_plus() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
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

    let leading_plus = expanded_formula_design("~ + condition + dose", &factors, &numeric).unwrap();
    let repeated_leading_plus =
        expanded_formula_design("~ + + condition + dose", &factors, &numeric).unwrap();
    let explicit = expanded_formula_design("~ condition + dose", &factors, &numeric).unwrap();
    assert_eq!(leading_plus, explicit);
    assert_eq!(repeated_leading_plus, explicit);

    let leading_plus_intercept =
        expanded_formula_design("~ + 1 + condition", &factors, &numeric).unwrap();
    let explicit_intercept =
        expanded_formula_design("~ 1 + condition", &factors, &numeric).unwrap();
    assert_eq!(leading_plus_intercept, explicit_intercept);

    let repeated_plus =
        expanded_formula_design("~ condition + + dose", &factors, &numeric).unwrap();
    assert_eq!(repeated_plus, explicit);

    let plus_minus =
        expanded_formula_design("~ condition + dose + - dose", &factors, &numeric).unwrap();
    let minus_plus =
        expanded_formula_design("~ condition + dose - + dose", &factors, &numeric).unwrap();
    let repeated_minus =
        expanded_formula_design("~ condition + dose - - dose", &factors, &numeric).unwrap();
    let condition_only = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(plus_minus, condition_only);
    assert_eq!(minus_plus, condition_only);
    assert_eq!(repeated_minus, condition_only);

    let plus_minus_intercept =
        expanded_formula_design("~ condition + - 1", &factors, &numeric).unwrap();
    let minus_plus_intercept =
        expanded_formula_design("~ condition - + 1", &factors, &numeric).unwrap();
    let no_intercept = expanded_formula_design("~ condition - 1", &factors, &numeric).unwrap();
    assert_eq!(plus_minus_intercept, no_intercept);
    assert_eq!(minus_plus_intercept, no_intercept);

    let minus_negative_one_intercept =
        expanded_formula_design("~ condition - - 1", &factors, &numeric).unwrap();
    assert_eq!(minus_negative_one_intercept, condition_only);

    assert!(expanded_formula_design("~ condition +", &factors, &numeric).is_err());
}

#[test]
fn expanded_formula_design_supports_pairwise_nested_terms() {
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
    let dose = [0.0, 1.0, 0.0, 1.0];
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

    let nested = expanded_formula_design("~ condition / batch", &factors, &[]).unwrap();
    let expanded = expanded_formula_design("~ condition + condition:batch", &factors, &[]).unwrap();
    assert_eq!(nested, expanded);
    assert_eq!(
        nested.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );

    let nested_numeric = expanded_formula_design("~ condition / dose", &factors, &numeric).unwrap();
    let expanded_numeric =
        expanded_formula_design("~ condition + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(nested_numeric, expanded_numeric);

    let no_intercept = expanded_formula_design("~ 0 + condition / batch", &factors, &[]).unwrap();
    assert_eq!(
        no_intercept.standard_design.coefficient_names().unwrap(),
        &[
            "conditionA".to_string(),
            "conditionB".to_string(),
            "conditionA:batch_Y_vs_X".to_string(),
            "conditionB:batch_Y_vs_X".to_string(),
        ]
    );
    assert_eq!(
        no_intercept.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 1.0,
        ]
    );

    let nested_in = expanded_formula_design("~ batch %in% condition", &factors, &[]).unwrap();
    let nested_in_explicit = expanded_formula_design("~ condition:batch", &factors, &[]).unwrap();
    assert_eq!(nested_in, nested_in_explicit);

    let nested_in_with_outer =
        expanded_formula_design("~ condition + batch %in% condition", &factors, &[]).unwrap();
    assert_eq!(nested_in_with_outer, nested);

    let nested_in_numeric =
        expanded_formula_design("~ dose %in% condition", &factors, &numeric).unwrap();
    let nested_in_numeric_explicit =
        expanded_formula_design("~ condition:dose", &factors, &numeric).unwrap();
    assert_eq!(nested_in_numeric, nested_in_numeric_explicit);

    assert!(expanded_formula_design("~ condition %in%", &factors, &[]).is_err());
    assert!(
        expanded_formula_design("~ condition %in% batch %in% dose", &factors, &numeric).is_err()
    );
}

#[test]
fn expanded_formula_design_supports_three_way_terms() {
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

    let direct = expanded_formula_design(
        "~ condition + batch + dose + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
            "dose".to_string(),
            "condition_B_vs_A:batch_Y_vs_X:dose".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, 2.0, 0.0, //
            1.0, 1.0, 1.0, 3.0, 3.0,
        ]
    );
    assert_eq!(
        direct.higher_order_interactions,
        vec!["condition:batch:dose".to_string()]
    );

    let star = expanded_formula_design("~ condition * batch * dose", &factors, &numeric).unwrap();
    assert!(
        star.standard_design
            .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose")
            .is_ok()
    );
    assert_eq!(
        star.higher_order_interactions,
        vec!["condition:batch:dose".to_string()]
    );

    let nested = expanded_formula_design("~ condition / batch / dose", &factors, &numeric).unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:batch + condition:batch:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested, expanded_nested);
}

#[test]
fn expanded_formula_design_supports_arbitrary_order_terms() {
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

    let direct = expanded_formula_design(
        "~ condition + batch + dose + time + condition:batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "batch_Y_vs_X".to_string(),
            "dose".to_string(),
            "time".to_string(),
            "condition_B_vs_A:batch_Y_vs_X:dose:time".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 2.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 3.0, 0.0, //
            1.0, 1.0, 0.0, 2.0, 5.0, 0.0, //
            1.0, 1.0, 1.0, 3.0, 7.0, 21.0,
        ]
    );
    assert_eq!(
        direct.higher_order_interactions,
        vec!["condition:batch:dose:time".to_string()]
    );

    let star =
        expanded_formula_design("~ condition * batch * dose * time", &factors, &numeric).unwrap();
    assert!(
        star.standard_design
            .coefficient_index("condition_B_vs_A:batch_Y_vs_X:dose:time")
            .is_ok()
    );
    assert!(
        star.higher_order_interactions
            .iter()
            .any(|interaction| interaction == "condition:batch:dose:time")
    );

    let nested =
        expanded_formula_design("~ condition / batch / dose / time", &factors, &numeric).unwrap();
    let expanded_nested = expanded_formula_design(
        "~ condition + condition:batch + condition:batch:dose + condition:batch:dose:time",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested, expanded_nested);
}
