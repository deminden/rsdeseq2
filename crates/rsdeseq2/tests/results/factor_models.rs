use super::*;

#[test]
fn fit_expanded_factor_beta_prior_wald_results_builds_design_and_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
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

    let factor = fit_expanded_factor_beta_prior_wald_results(
        ExpandedFactorBetaPriorWaldResultsInput {
            counts: &counts,
            factor: "condition",
            sample_levels: &sample_levels,
            reference: "A",
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

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.fit, direct.fit);
    assert_eq!(factor.results, direct.results);
    assert_eq!(
        factor.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let factor_contrast = fit_expanded_factor_beta_prior_wald_contrast_results(
        ExpandedFactorBetaPriorWaldResultsInput {
            counts: &counts,
            factor: "condition",
            sample_levels: &sample_levels,
            reference: "A",
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        &[0.0, 1.0],
    )
    .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &factor.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(factor_contrast.fit, factor.fit);
    assert_eq!(factor_contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_factor_beta_prior_wald_replacement_builds_design_and_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
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

    let factor = fit_expanded_factor_beta_prior_wald_results_with_cooks_replacement(
        ExpandedFactorBetaPriorWaldResultsInput {
            counts: &counts,
            factor: "condition",
            sample_levels: &sample_levels,
            reference: "A",
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

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.replacement, direct);

    let factor_contrast =
        fit_expanded_factor_beta_prior_wald_contrast_results_with_cooks_replacement(
            ExpandedFactorBetaPriorWaldResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
                size_factors: &size_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options,
            },
            &[0.0, 1.0],
            &replacement_options,
        )
        .unwrap();
    assert!(factor_contrast.replacement.refit.is_some());
    assert_eq!(
        factor_contrast
            .replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_factor_beta_prior_wald_results_accepts_normalization_factors_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
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
        1,
    )
    .unwrap();

    let factor =
        fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options,
            },
            1,
        )
        .unwrap();

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.fit, direct.fit);
    assert_eq!(factor.results, direct.results);

    let factor_contrast =
        fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
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
            &[0.0, 1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &factor.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(factor_contrast.fit, factor.fit);
    assert_eq!(factor_contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_factor_beta_prior_wald_normalization_factor_replacement_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
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

    let factor =
        fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
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

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.replacement, direct);

    let factor_contrast =
        fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options,
            },
            &[0.0, 1.0],
            &replacement_options,
        )
        .unwrap();
    assert!(factor_contrast.replacement.refit.is_some());
    assert_eq!(
        factor_contrast
            .replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("contrast")
    );
}
