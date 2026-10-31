use super::*;

#[test]
fn fit_with_dispersion_wraps_fixed_dispersion_irls() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let options = no_ridge_options();

    let wrapped = fit_with_dispersion(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        options.clone(),
    )
    .unwrap();
    let direct =
        fit_fixed_dispersion_irls(&counts, &design, &[1.0, 1.0, 1.0, 1.0], &[0.05], options)
            .unwrap();

    assert_eq!(wrapped, direct);
}

#[test]
fn irls_optim_fallback_refits_nonconverged_rows_by_default_and_can_be_disabled() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let without_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: false,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let with_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(!without_optim.beta_converged[0]);
    assert!(with_optim.beta_converged[0]);
    assert_relative_eq!(
        with_optim.beta.as_slice()[0],
        10.0_f64.log2(),
        epsilon = 1e-4
    );
    assert_relative_eq!(
        with_optim.beta.as_slice()[1],
        2.0_f64.log2(),
        epsilon = 1e-4
    );
    assert!(with_optim.log_like[0] >= without_optim.log_like[0] - 1e-8);
}

#[test]
fn irls_optim_fallback_preserves_pre_optim_hat_diagonal() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let pre_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: false,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let with_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: true,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(with_optim.beta_converged[0]);
    for (actual, expected) in with_optim
        .hat_diagonal
        .as_slice()
        .iter()
        .zip(pre_optim.hat_diagonal.as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_optim_fallback_log_likelihood_uses_min_mu_floored_means() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 1, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let fit = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: true,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let stored_mu = fit.mu.row(0).unwrap();
    let inference_mu = stored_mu
        .iter()
        .copied()
        .map(|value| value.max(IrlsOptions::default().min_mu))
        .collect::<Vec<_>>();
    let expected_log_like =
        nbinom_log_likelihood(counts.row(0).unwrap(), &inference_mu, 0.05).unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.log_like[0], expected_log_like, epsilon = 1e-12);
}

#[test]
fn irls_qr_solver_matches_normal_equations_for_two_group_design() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_qr_options(),
    )
    .unwrap();

    assert_eq!(normal.beta_converged, qr.beta_converged);
    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.beta_se.as_slice().iter().zip(normal.beta_se.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn original_qr_results_match_non_qr_results_on_stable_rows() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_qr_options(),
    )
    .unwrap();

    let normal_wald = wald_test_coefficient(&normal, 1).unwrap();
    let qr_wald = wald_test_coefficient(&qr, 1).unwrap();

    for gene in 0..counts.n_genes() {
        assert_relative_eq!(
            qr.beta.row(gene).unwrap()[1],
            normal.beta.row(gene).unwrap()[1],
            epsilon = 1e-8,
            max_relative = 1e-8
        );
        assert_relative_eq!(
            qr.beta_se.row(gene).unwrap()[1],
            normal.beta_se.row(gene).unwrap()[1],
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
    for (actual, expected) in qr_wald.stat.iter().zip(normal_wald.stat.iter()) {
        assert_relative_eq!(
            actual.unwrap(),
            expected.unwrap(),
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
    for (actual, expected) in qr_wald.pvalue.iter().zip(normal_wald.pvalue.iter()) {
        assert_relative_eq!(
            actual.unwrap(),
            expected.unwrap(),
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
}

#[test]
fn irls_qr_solver_handles_default_ridge() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let fit = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            solver: IrlsSolver::Qr,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-5);
    assert!(
        fit.hat_diagonal
            .as_slice()
            .iter()
            .all(|value| value.is_finite())
    );
}

#[test]
fn irls_vector_ridge_matches_scalar_ridge_when_values_match() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let scalar = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.25,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let vector = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.25, 0.25]),
    )
    .unwrap();

    for (actual, expected) in vector.beta.as_slice().iter().zip(scalar.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-10, max_relative = 1e-10);
    }
    for (actual, expected) in vector
        .beta_covariance
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .zip(scalar.beta_covariance.as_ref().unwrap().as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-10, max_relative = 1e-10);
    }
}

#[test]
fn irls_vector_ridge_can_penalize_one_coefficient() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let unpenalized = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        no_ridge_options(),
    )
    .unwrap();
    let penalized = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();

    assert!(penalized.beta_converged[0]);
    assert!(penalized.beta.as_slice()[1].abs() < unpenalized.beta.as_slice()[1].abs());
    assert_relative_eq!(
        penalized.beta.as_slice()[0],
        15.0_f64.log2(),
        epsilon = 2e-1
    );
}

#[test]
fn irls_qr_solver_matches_normal_equations_with_vector_ridge() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            solver: IrlsSolver::Qr,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();

    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_all_one_weights_match_unweighted_fit() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(2, 4, vec![1.0; 8]).unwrap();

    let unweighted = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let weighted = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();

    assert_eq!(weighted.beta_converged, unweighted.beta_converged);
    for (actual, expected) in weighted
        .beta
        .as_slice()
        .iter()
        .zip(unweighted.beta.as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in weighted.mu.as_slice().iter().zip(unweighted.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in weighted.log_like.iter().zip(unweighted.log_like.iter()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_zero_weight_sample_does_not_influence_fit() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 1_000, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 1.0, 1.0]).unwrap();

    let fit = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[0], 10.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.mu.as_slice()[1], 10.0, epsilon = 1e-7);
}

#[test]
fn irls_qr_solver_matches_weighted_normal_equations() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 1_000, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 1.0, 1.0]).unwrap();

    let normal = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_qr_options(),
    )
    .unwrap();

    assert_eq!(qr.beta_converged, normal.beta_converged);
    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(
        fit_fixed_dispersion_irls(&counts, &design, &[1.0], &[0.1], no_ridge_options()).is_err()
    );
    assert!(
        fit_fixed_dispersion_irls(&counts, &design, &[1.0, 1.0, 1.0], &[], no_ridge_options())
            .is_err()
    );

    let bad_design = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();
    assert!(
        fit_fixed_dispersion_irls(
            &counts,
            &bad_design,
            &[1.0, 1.0, 1.0],
            &[0.1],
            no_ridge_options()
        )
        .is_err()
    );

    let bad_weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0, -1.0, 1.0]).unwrap();
    assert!(
        fit_fixed_dispersion_irls_with_weights(
            &counts,
            &design,
            &[1.0, 1.0, 1.0],
            &[0.1],
            Some(&bad_weights),
            no_ridge_options()
        )
        .is_err()
    );

    assert!(
        fit_fixed_dispersion_irls(
            &counts,
            &design,
            &[1.0, 1.0, 1.0],
            &[0.1],
            no_ridge_options().ridge_lambda_by_coefficient(vec![0.0, 0.0])
        )
        .is_err()
    );
    assert!(
        fit_fixed_dispersion_irls(
            &counts,
            &design,
            &[1.0, 1.0, 1.0],
            &[0.1],
            no_ridge_options().ridge_lambda_by_coefficient(vec![-1.0])
        )
        .is_err()
    );
}

#[test]
fn irls_rejects_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![0, 0, 0]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();
    let err = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("all zero"));
}

#[test]
fn ordered_batches_report_the_lowest_failed_gene_at_one_and_four_workers() {
    for early in [3, 2051] {
        let mut values = vec![10_u32; 4097 * 8];
        values[early * 8..(early + 1) * 8].fill(0);
        values[3000 * 8..3001 * 8].fill(0);
        let counts = CountMatrix::from_row_major_u32(4097, 8, values).unwrap();
        let design = DesignMatrix::from_row_major(
            8,
            2,
            (0..8).flat_map(|i| [1.0, (i / 4) as f64]).collect(),
            None,
        )
        .unwrap();
        for workers in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            let error = pool
                .install(|| {
                    rsdeseq2::glm::fit_fixed_dispersion_irls(
                        &counts,
                        &design,
                        &[1.0; 8],
                        &vec![0.1; 4097],
                        IrlsOptions::default(),
                    )
                })
                .unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains(&format!("gene {early} is all zero")),
                "{error}"
            );
        }
    }
}
