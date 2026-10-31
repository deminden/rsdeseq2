use super::*;

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
fn fit_irls_dispatches_intercept_only_shortcut() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![6, 6, 12, 12]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        1,
        vec![1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();

    let wrapped = fit_irls(
        &counts,
        &design,
        &[1.0, 1.0, 2.0, 2.0],
        &[0.2],
        IrlsOptions::default(),
    )
    .unwrap();
    let shortcut =
        fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 2.0, 2.0], &[0.2]).unwrap();

    assert_eq!(wrapped, shortcut);
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
fn fit_irls_dispatches_general_irls_for_multi_coefficient_design() {
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

    let wrapped = fit_irls(
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
fn estimate_beta_wraps_fixed_dispersion_beta_dispatch() {
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

    let beta_fit = estimate_beta(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        options.clone(),
    )
    .unwrap();
    let irls_fit = fit_irls(&counts, &design, &[1.0, 1.0, 1.0, 1.0], &[0.05], options).unwrap();

    assert_eq!(beta_fit, irls_fit);
}
