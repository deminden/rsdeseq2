use super::*;

#[test]
fn fit_expanded_beta_prior_wald_results_with_cooks_replacement_refits_marked_rows() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
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
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let groups = [vec![0], vec![1, 2]];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let design = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded_design,
        standard_design: &standard_design,
        coefficient_groups: &groups,
    };
    let input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options,
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };

    let output = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        input.clone(),
        1,
        &replacement_options,
    )
    .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(output.refit.is_some());
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        &[10, 12, 41, 24, 30, 33, 40, 54]
    );
    for gene in output.refit_plan.refit_rows.iter().copied() {
        let refit = output.refit.as_ref().unwrap();
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit.results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].base_mean,
            output.refit_plan.replaced_base_mean[gene]
        );
    }

    let contrast_output = fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
        input,
        &[0.0, 1.0],
        &replacement_options,
    )
    .unwrap();
    assert!(contrast_output.refit.is_some());
    assert_eq!(
        contrast_output.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_beta_prior_wald_results_with_cooks_replacement_skips_unmarked_rows() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
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
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let groups = [vec![0], vec![1, 2]];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: ExpandedModelBetaPriorDesignInput {
            expanded_design: &expanded_design,
            standard_design: &standard_design,
            coefficient_groups: &groups,
        },
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: None,
        options,
    };

    let output = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        input,
        1,
        &CooksReplacementOptions::new(f64::MAX),
    )
    .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit.is_none());
    assert_eq!(output.refit_plan.replacement.replaced_counts, counts);
    assert_eq!(output.results.rows[0].base_mean, base_mean[0]);
}

#[test]
fn fit_expanded_beta_prior_wald_normalization_factor_cooks_replacement_refits_marked_rows() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
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
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
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
    let groups = [vec![0], vec![1, 2]];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design: ExpandedModelBetaPriorDesignInput {
            expanded_design: &expanded_design,
            standard_design: &standard_design,
            coefficient_groups: &groups,
        },
        normalization_factors: &normalization_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: None,
        options,
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };

    let output =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            input.clone(),
            1,
            &replacement_options,
        )
        .unwrap();

    assert!(output.refit_plan.should_refit);
    assert!(output.refit.is_some());
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        &[10, 12, 41, 24, 30, 33, 40, 54]
    );

    let contrast_output =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            input,
            &[0.0, 1.0],
            &replacement_options,
        )
        .unwrap();
    assert!(contrast_output.refit.is_some());
    assert_eq!(
        contrast_output.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_beta_prior_wald_results_accepts_offsets_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
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
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
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
    let groups = [vec![0], vec![1, 2]];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let design = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded_design,
        standard_design: &standard_design,
        coefficient_groups: &groups,
    };
    let input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design,
        normalization_factors: &normalization_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options,
    };

    let coefficient_workflow =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            input.clone(),
            1,
        )
        .unwrap();
    let direct_fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
            },
            &dispersions,
            &base_mean,
            &disp_fit,
            input.options.clone(),
        )
        .unwrap();
    let direct_results = build_wald_results_from_expanded_beta_prior_fit(
        &base_mean,
        &direct_fit,
        1,
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(coefficient_workflow.fit, direct_fit);
    assert_eq!(coefficient_workflow.results, direct_results);

    let contrast_workflow =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            input,
            &[0.0, 1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &coefficient_workflow.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast_workflow.fit, coefficient_workflow.fit);
    assert_eq!(contrast_workflow.results, direct_contrast_results);
}

#[test]
fn fit_expanded_beta_prior_wald_cooks_replacement_carries_weights_into_refit() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
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
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 0.8, 0.25, 0.9, //
            0.7, 1.0, 0.95, 0.85,
        ],
    )
    .unwrap();
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let groups = [vec![0], vec![1, 2]];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: ExpandedModelBetaPriorDesignInput {
            expanded_design: &expanded_design,
            standard_design: &standard_design,
            coefficient_groups: &groups,
        },
        size_factors: &size_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: None,
        options,
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };

    let output = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        input.clone(),
        1,
        &replacement_options,
    )
    .unwrap();
    assert!(output.refit.is_some());

    let direct_refit = fit_expanded_beta_prior_wald_results(
        ExpandedBetaPriorWaldResultsInput {
            counts: &output.refit_plan.replacement.replaced_counts,
            design: input.design,
            size_factors: &size_factors,
            weights: Some(&weights),
            dispersions: &dispersions,
            base_mean: &output.refit_plan.replaced_base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: input.options,
        },
        1,
    )
    .unwrap();

    assert_eq!(output.refit.as_ref().unwrap(), &direct_refit);
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            direct_refit.results.rows[gene].log2_fold_change
        );
    }
}
