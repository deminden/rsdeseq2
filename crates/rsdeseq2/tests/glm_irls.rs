use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn no_ridge_options() -> IrlsOptions {
    IrlsOptions {
        ridge_lambda: 0.0,
        ..IrlsOptions::default()
    }
}

fn no_ridge_qr_options() -> IrlsOptions {
    IrlsOptions {
        solver: IrlsSolver::Qr,
        ..no_ridge_options()
    }
}

#[test]
fn irls_intercept_only_agrees_with_intercept_shortcut_for_equal_size_factors() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![4, 4, 8, 8]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        1,
        vec![1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();

    let irls = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options(),
    )
    .unwrap();
    let shortcut =
        fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 1.0, 1.0], &[0.1]).unwrap();

    assert!(irls.beta_converged[0]);
    assert_relative_eq!(
        irls.beta.as_slice()[0],
        shortcut.beta.as_slice()[0],
        epsilon = 1e-8
    );
    for (actual, expected) in irls.mu.as_slice().iter().zip(shortcut.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-7);
    }
}

#[test]
fn irls_two_group_design_recovers_known_log2_fold_change() {
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
        no_ridge_options(),
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[0], 10.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-8);
    assert_eq!(fit.n_terms, 2);
    assert_eq!(
        fit.model_matrix.coefficient_names().unwrap()[1],
        "condition_B_vs_A"
    );
    assert_relative_eq!(fit.mu.as_slice()[0], 10.0, epsilon = 1e-7);
    assert_relative_eq!(fit.mu.as_slice()[2], 20.0, epsilon = 1e-7);
    assert!(fit.beta_iter[0] < no_ridge_options().maxit);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_cols(), 4);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[0].sqrt(),
        fit.beta_se.row(0).unwrap()[0],
        epsilon = 1e-12
    );
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
    assert!(fit
        .hat_diagonal
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
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
    assert!(fit_fixed_dispersion_irls(
        &counts,
        &bad_design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options()
    )
    .is_err());

    let bad_weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0, -1.0, 1.0]).unwrap();
    assert!(fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        Some(&bad_weights),
        no_ridge_options()
    )
    .is_err());

    assert!(fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options().ridge_lambda_by_coefficient(vec![0.0, 0.0])
    )
    .is_err());
    assert!(fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options().ridge_lambda_by_coefficient(vec![-1.0])
    )
    .is_err());
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
