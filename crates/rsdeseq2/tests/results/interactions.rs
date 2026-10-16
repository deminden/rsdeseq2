use super::*;

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_numeric_covariates() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let direct_design = expanded_additive_design(&factors, &numeric_covariates).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 2).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        2,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("dose")
    );
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_factor_interactions() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
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
    let interactions = [ExpandedFactorInteractionSpec {
        left_factor: "condition",
        right_factor: "batch",
    }];
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let direct_design =
        expanded_additive_design_with_interactions(&factors, &[], &interactions).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 3).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &interactions,
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        3,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A:batch_Y_vs_X")
    );
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_numeric_interactions() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let time = [1.0, 1.0, 2.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];
    let factor_numeric = [ExpandedFactorNumericInteractionSpec {
        factor: "condition",
        numeric: "dose",
    }];
    let numeric_interactions = [ExpandedNumericInteractionSpec {
        left_numeric: "dose",
        right_numeric: "time",
    }];
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let direct_design = expanded_additive_design_with_all_interactions(
        &factors,
        &numeric_covariates,
        &[],
        &factor_numeric,
        &numeric_interactions,
    )
    .unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 4).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            interactions: &[],
            factor_numeric_interactions: &factor_numeric,
            numeric_interactions: &numeric_interactions,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        4,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A:dose")
    );
}
