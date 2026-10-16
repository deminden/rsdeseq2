use super::*;

#[test]
fn fit_expanded_additive_beta_prior_wald_results_builds_design_and_matches_direct_workflow() {
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
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
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
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 1).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options: options.clone(),
        },
        1,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let contrast = fit_expanded_additive_beta_prior_wald_contrast_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
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
        &[0.0, 1.0, -1.0],
    )
    .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &additive.fit,
        &[0.0, 1.0, -1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast.fit, additive.fit);
    assert_eq!(contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_additive_beta_prior_wald_replacement_builds_design_and_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
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
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
    let direct = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: &counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &direct_design.expanded_design,
                standard_design: &direct_design.standard_design,
                coefficient_groups: &direct_design.coefficient_groups,
            },
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: options.clone(),
        },
        1,
        &replacement_options,
    )
    .unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: options.clone(),
        },
        1,
        &replacement_options,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.replacement, direct);

    let contrast = fit_expanded_additive_beta_prior_wald_contrast_results_with_cooks_replacement(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options,
        },
        &[0.0, 1.0, -1.0],
        &replacement_options,
    )
    .unwrap();
    assert!(contrast.replacement.refit.is_some());
    assert_eq!(
        contrast.replacement.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_normalization_factors_and_weights() {
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
    let normalization_factors = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 0.9, 1.0, 0.8, //
            1.0, 1.0, 0.95, 0.9,
        ],
    )
    .unwrap();
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
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design: direct_design_input,
        normalization_factors: &normalization_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        direct_input,
        2,
    )
    .unwrap();

    let additive =
        fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
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

    let contrast =
        fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: BetaPriorRefitOptions {
                    fit_options: IrlsOptions::default(),
                    variance_options: BetaPriorVarianceOptions {
                        method: BetaPriorVarianceMethod::Quantile,
                        upper_quantile: 0.5,
                        ..BetaPriorVarianceOptions::default()
                    },
                },
            },
            &[0.0, 1.0, -1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &additive.fit,
        &[0.0, 1.0, -1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast.fit, additive.fit);
    assert_eq!(contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_additive_beta_prior_wald_normalization_factor_replacement_matches_direct_workflow()
{
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
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
    let normalization_factors = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
    let direct =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &direct_design.expanded_design,
                    standard_design: &direct_design.standard_design,
                    coefficient_groups: &direct_design.coefficient_groups,
                },
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options: options.clone(),
            },
            1,
            &replacement_options,
        )
        .unwrap();

    let additive =
        fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options: options.clone(),
            },
            1,
            &replacement_options,
        )
        .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.replacement, direct);

    let contrast =
        fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options,
            },
            &[0.0, 1.0, -1.0],
            &replacement_options,
        )
        .unwrap();
    assert!(contrast.replacement.refit.is_some());
    assert_eq!(
        contrast.replacement.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}
