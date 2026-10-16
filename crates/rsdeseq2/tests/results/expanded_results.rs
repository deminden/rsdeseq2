use super::*;

#[test]
fn build_wald_results_from_expanded_model_fit_reports_collapsed_coefficients() {
    let expanded_design = DesignMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_A".into(),
            "condition_B".into(),
            "condition_C".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, true],
        beta: RowMajorMatrix::from_row_major(
            2,
            4,
            vec![
                4.0, 1.0, 3.0, 5.0, //
                6.0, -2.0, 2.0, 4.0,
            ],
        )
        .unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 4, vec![1.0; 8]).unwrap(),
        beta_optim_start: RowMajorMatrix::from_elem(2, 4, f64::NAN).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                16,
                vec![
                    4.0, 1.0, 2.0, 3.0, //
                    1.0, 9.0, 4.0, 5.0, //
                    2.0, 4.0, 16.0, 6.0, //
                    3.0, 5.0, 6.0, 25.0, //
                    1.0, 0.0, 0.0, 0.0, //
                    0.0, 4.0, 1.0, 1.0, //
                    0.0, 1.0, 9.0, 1.0, //
                    0.0, 1.0, 1.0, 16.0,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        beta_optim_iter: vec![f64::NAN; 2],
        beta_optim_start_objective: vec![f64::NAN; 2],
        beta_optim_objective: vec![f64::NAN; 2],
        beta_optim_gradient_norm: vec![f64::NAN; 2],
        model_matrix: expanded_design,
        n_terms: 4,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_wald_results_from_expanded_model_fit(
        &[10.0, 20.0],
        &expanded_fit,
        &standard_design,
        &[vec![0], vec![1, 2]],
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_relative_eq!(
        results.rows[0].lfc_se.unwrap(),
        8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].stat.unwrap(),
        2.0 / 8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn build_wald_contrast_results_from_expanded_model_fit_reports_collapsed_contrast() {
    let expanded_design = DesignMatrix::from_row_major(
        2,
        3,
        vec![
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_A".into(),
            "condition_B".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, true],
        beta: RowMajorMatrix::from_row_major(
            2,
            3,
            vec![
                4.0, 1.0, 3.0, //
                6.0, -2.0, 2.0,
            ],
        )
        .unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 3, vec![1.0; 6]).unwrap(),
        beta_optim_start: RowMajorMatrix::from_elem(2, 3, f64::NAN).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                9,
                vec![
                    4.0, 1.0, 2.0, //
                    1.0, 9.0, 4.0, //
                    2.0, 4.0, 16.0, //
                    1.0, 0.0, 0.0, //
                    0.0, 4.0, 1.0, //
                    0.0, 1.0, 9.0,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        beta_optim_iter: vec![f64::NAN; 2],
        beta_optim_start_objective: vec![f64::NAN; 2],
        beta_optim_objective: vec![f64::NAN; 2],
        beta_optim_gradient_norm: vec![f64::NAN; 2],
        model_matrix: expanded_design,
        n_terms: 3,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_wald_contrast_results_from_expanded_model_fit(
        &[10.0, 20.0],
        &expanded_fit,
        &standard_design,
        &[vec![0], vec![1, 2]],
        &[0.0, 1.0],
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_relative_eq!(
        results.rows[0].lfc_se.unwrap(),
        8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].stat.unwrap(),
        2.0 / 8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(results.metadata.result_name.as_deref(), Some("contrast"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
}

#[test]
fn build_wald_results_from_expanded_beta_prior_fit_uses_collapsed_prior_fit() {
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let prior_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, false],
        beta: RowMajorMatrix::from_row_major(2, 2, vec![4.0, 2.0, 6.0, -1.0]).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 2, vec![0.5, 0.25, 0.75, 0.5]).unwrap(),
        beta_optim_start: RowMajorMatrix::from_elem(2, 2, f64::NAN).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                4,
                vec![
                    0.25, 0.0, //
                    0.0, 0.0625, //
                    0.5625, 0.0, //
                    0.0, 0.25,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        beta_optim_iter: vec![f64::NAN; 2],
        beta_optim_start_objective: vec![f64::NAN; 2],
        beta_optim_objective: vec![f64::NAN; 2],
        beta_optim_gradient_norm: vec![f64::NAN; 2],
        model_matrix: standard_design,
        n_terms: 2,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let fit = ExpandedModelBetaPriorGlmFit {
        expanded_mle_fit: prior_fit.clone(),
        expanded_prior_fit: prior_fit.clone(),
        prior_fit: prior_fit.clone(),
        beta_prior_variance: vec![1e6, 1.0],
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let coefficient_results = build_wald_results_from_expanded_beta_prior_fit(
        &[10.0, 20.0],
        &fit,
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();
    let direct_results = build_wald_results(
        &[10.0, 20.0],
        &prior_fit,
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(coefficient_results, direct_results);
    assert_eq!(
        coefficient_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &[10.0, 20.0],
        &fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();
    let direct_contrast = wald_test_contrast(&prior_fit, &[0.0, 1.0]).unwrap();
    let direct_contrast_results = build_wald_contrast_results(
        &[10.0, 20.0],
        &prior_fit,
        &direct_contrast,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(contrast_results, direct_contrast_results);
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_beta_prior_wald_results_runs_fit_and_result_workflow() {
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

    let coefficient_workflow = fit_expanded_beta_prior_wald_results(input.clone(), 1).unwrap();
    let direct_fit = fit_expanded_glms_with_estimated_beta_prior_variance(
        &counts,
        design,
        &size_factors,
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
    assert_eq!(
        coefficient_workflow.results.metadata.result_name.as_deref(),
        Some("condition_or_batch")
    );

    let contrast_workflow =
        fit_expanded_beta_prior_wald_contrast_results(input, &[0.0, 1.0]).unwrap();
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
    assert_eq!(
        contrast_workflow.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}
