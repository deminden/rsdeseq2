use super::*;

#[test]
fn expanded_model_beta_average_collapses_grouped_coefficients() {
    let expanded = RowMajorMatrix::from_row_major(
        2,
        5,
        vec![
            4.0, 1.0, 3.0, 10.0, 14.0, //
            6.0, -2.0, 2.0, 5.0, 9.0,
        ],
    )
    .unwrap();

    let collapsed =
        average_expanded_model_coefficients(&expanded, &[vec![0], vec![1, 2], vec![3, 4]]).unwrap();

    assert_eq!(collapsed.n_rows(), 2);
    assert_eq!(collapsed.n_cols(), 3);
    assert_relative_eq!(collapsed.row(0).unwrap()[0], 4.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.row(0).unwrap()[1], 2.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.row(0).unwrap()[2], 12.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.row(1).unwrap()[0], 6.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.row(1).unwrap()[1], 0.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.row(1).unwrap()[2], 7.0, epsilon = 1e-12);
}

#[test]
fn expanded_model_covariance_average_propagates_grouped_uncertainty() {
    let covariance = RowMajorMatrix::from_row_major(
        1,
        16,
        vec![
            4.0, 1.0, 2.0, 3.0, //
            1.0, 9.0, 4.0, 5.0, //
            2.0, 4.0, 16.0, 6.0, //
            3.0, 5.0, 6.0, 25.0,
        ],
    )
    .unwrap();

    let collapsed =
        average_expanded_model_covariances(&covariance, 4, &[vec![0, 1], vec![2, 3]]).unwrap();

    assert_eq!(collapsed.n_rows(), 1);
    assert_eq!(collapsed.n_cols(), 4);
    let row = collapsed.row(0).unwrap();
    assert_relative_eq!(row[0], (4.0 + 1.0 + 1.0 + 9.0) / 4.0, epsilon = 1e-12);
    assert_relative_eq!(row[1], (2.0 + 3.0 + 4.0 + 5.0) / 4.0, epsilon = 1e-12);
    assert_relative_eq!(row[2], (2.0 + 4.0 + 3.0 + 5.0) / 4.0, epsilon = 1e-12);
    assert_relative_eq!(row[3], (16.0 + 6.0 + 6.0 + 25.0) / 4.0, epsilon = 1e-12);
}

#[test]
fn expanded_model_fit_collapse_returns_standard_glm_surface() {
    let expanded_design = DesignMatrix::from_row_major(
        3,
        4,
        vec![
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 0.0, 1.0,
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
        3,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 0.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, false],
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
        mu: RowMajorMatrix::from_row_major(2, 3, vec![10.0, 20.0, 30.0, 5.0, 6.0, 7.0]).unwrap(),
        beta_iter: vec![5, 6],
        beta_optim_iter: vec![f64::NAN; 2],
        beta_optim_start_objective: vec![f64::NAN; 2],
        beta_optim_objective: vec![f64::NAN; 2],
        beta_optim_gradient_norm: vec![f64::NAN; 2],
        model_matrix: expanded_design,
        n_terms: 4,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 3, vec![0.1; 6]).unwrap(),
    };

    let collapsed =
        collapse_expanded_model_fit(&expanded_fit, &standard_design, &[vec![0], vec![1, 2]])
            .unwrap();

    assert_eq!(collapsed.model_matrix, standard_design);
    assert_eq!(collapsed.n_terms, 2);
    assert_eq!(collapsed.log_like, expanded_fit.log_like);
    assert_eq!(collapsed.beta_converged, expanded_fit.beta_converged);
    assert_eq!(collapsed.mu, expanded_fit.mu);
    assert_eq!(collapsed.hat_diagonal, expanded_fit.hat_diagonal);
    assert_relative_eq!(collapsed.beta.row(0).unwrap()[0], 4.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.beta.row(0).unwrap()[1], 2.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.beta.row(1).unwrap()[0], 6.0, epsilon = 1e-12);
    assert_relative_eq!(collapsed.beta.row(1).unwrap()[1], 0.0, epsilon = 1e-12);

    let covariance = collapsed.beta_covariance.as_ref().unwrap();
    assert_eq!(covariance.n_cols(), 4);
    assert_relative_eq!(covariance.row(0).unwrap()[0], 4.0, epsilon = 1e-12);
    assert_relative_eq!(covariance.row(0).unwrap()[3], 8.25, epsilon = 1e-12);
    assert_relative_eq!(collapsed.beta_se.row(0).unwrap()[0], 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        collapsed.beta_se.row(0).unwrap()[1],
        8.25_f64.sqrt(),
        epsilon = 1e-12
    );
}

#[test]
fn expanded_model_group_contrast_uses_average_weights() {
    let contrast = expanded_model_group_contrast(5, &[1, 2], &[3, 4]).unwrap();

    assert_eq!(contrast, vec![0.0, 0.5, 0.5, -0.5, -0.5]);

    let expanded = RowMajorMatrix::from_row_major(1, 5, vec![4.0, 1.0, 3.0, 10.0, 14.0]).unwrap();
    let beta = expanded.row(0).unwrap();
    let effect = beta
        .iter()
        .zip(contrast.iter())
        .map(|(beta, weight)| beta * weight)
        .sum::<f64>();
    assert_relative_eq!(effect, 2.0 - 12.0, epsilon = 1e-12);
}

#[test]
fn expanded_model_group_helpers_validate_inputs() {
    let expanded = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap();

    assert!(average_expanded_model_coefficients(&expanded, &[]).is_err());
    assert!(average_expanded_model_coefficients(&expanded, &[vec![0, 2]]).is_err());
    assert!(average_expanded_model_coefficients(&expanded, &[vec![0, 0]]).is_err());

    let nonfinite = RowMajorMatrix::from_row_major(1, 2, vec![1.0, f64::INFINITY]).unwrap();
    assert!(average_expanded_model_coefficients(&nonfinite, &[vec![0, 1]]).is_err());

    let covariance = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
    assert!(average_expanded_model_covariances(&covariance, 0, &[vec![0]]).is_err());
    assert!(average_expanded_model_covariances(&covariance, 3, &[vec![0]]).is_err());
    assert!(average_expanded_model_covariances(&covariance, 2, &[vec![0, 2]]).is_err());
    let nonfinite_covariance =
        RowMajorMatrix::from_row_major(1, 4, vec![1.0, f64::NAN, 0.0, 1.0]).unwrap();
    assert!(average_expanded_model_covariances(&nonfinite_covariance, 2, &[vec![0, 1]]).is_err());
    let expanded_design = DesignMatrix::from_row_major(1, 2, vec![1.0, 0.0], None).unwrap();
    let standard_design = DesignMatrix::from_row_major(1, 1, vec![1.0], None).unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![0.0],
        beta_converged: vec![true],
        beta: expanded.clone(),
        beta_se: expanded.clone(),
        beta_optim_start: RowMajorMatrix::from_elem(1, 2, f64::NAN).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(1, 1, vec![1.0]).unwrap(),
        beta_iter: vec![1],
        beta_optim_iter: vec![f64::NAN],
        beta_optim_start_objective: vec![f64::NAN],
        beta_optim_objective: vec![f64::NAN],
        beta_optim_gradient_norm: vec![f64::NAN],
        model_matrix: expanded_design,
        n_terms: 2,
        hat_diagonal: RowMajorMatrix::from_row_major(1, 1, vec![0.0]).unwrap(),
    };
    assert!(collapse_expanded_model_fit(&expanded_fit, &standard_design, &[vec![0]]).is_err());

    assert!(expanded_model_group_contrast(2, &[], &[1]).is_err());
    assert!(expanded_model_group_contrast(2, &[0], &[0]).is_err());
    assert!(expanded_model_group_contrast(2, &[0, 0], &[1]).is_err());
    assert!(expanded_model_group_contrast(2, &[0], &[2]).is_err());
}

#[test]
fn expanded_beta_prior_refit_runs_expanded_fit_and_collapses_standard_surface() {
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
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let groups = [vec![0], vec![1, 2]];

    let wrapped = fit_expanded_glms_with_estimated_beta_prior_variance(
        &counts,
        ExpandedModelBetaPriorDesignInput {
            expanded_design: &expanded_design,
            standard_design: &standard_design,
            coefficient_groups: &groups,
        },
        &size_factors,
        &dispersions,
        &base_mean,
        &disp_fit,
        options.clone(),
    )
    .unwrap();

    let expanded_mle = fit_irls(
        &counts,
        &expanded_design,
        &size_factors,
        &dispersions,
        options.fit_options.clone(),
    )
    .unwrap();
    let beta_prior_variance = estimate_beta_prior_variance(
        &expanded_mle.beta,
        &base_mean,
        &disp_fit,
        expanded_design.coefficient_names(),
        options.variance_options,
    )
    .unwrap();
    let expanded_prior = fit_glms_with_beta_prior_variance(
        &counts,
        &expanded_design,
        &size_factors,
        &dispersions,
        &beta_prior_variance,
        options.fit_options,
    )
    .unwrap();
    let collapsed =
        collapse_expanded_model_fit(&expanded_prior, &standard_design, &groups).unwrap();

    assert_eq!(wrapped.expanded_mle_fit, expanded_mle);
    assert_eq!(wrapped.expanded_prior_fit, expanded_prior);
    assert_eq!(wrapped.prior_fit, collapsed);
    assert_eq!(wrapped.beta_prior_variance, beta_prior_variance);
    assert_eq!(wrapped.prior_fit.model_matrix, standard_design);
    assert_eq!(wrapped.prior_fit.n_terms, 2);
    for gene in 0..counts.n_genes() {
        let expanded_beta = wrapped.expanded_prior_fit.beta.row(gene).unwrap();
        assert_relative_eq!(
            wrapped.prior_fit.beta.row(gene).unwrap()[1],
            (expanded_beta[1] + expanded_beta[2]) / 2.0,
            epsilon = 1e-12
        );
    }
}
