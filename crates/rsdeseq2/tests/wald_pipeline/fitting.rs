use super::*;

#[test]
fn fixed_dispersion_wald_pipeline_uses_intercept_shortcut() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        3,
        vec![2, 4, 6, 10, 10, 10],
        Some(vec!["gene_a".into(), "gene_b".into()]),
        None,
    )
    .unwrap();
    let design =
        DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], Some(vec!["Intercept".into()]))
            .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1, 0.2], 0)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap().n_coefficients(), 1);
    assert_eq!(fit.dispersion.as_deref(), Some(&[0.1, 0.2][..]));
    assert_eq!(fit.beta_converged.as_ref().unwrap(), &vec![true, true]);
    assert_eq!(fit.wald.as_ref().unwrap().stat.len(), 2);
    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[0],
        4.0_f64.log2(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        4.0_f64.log2(),
        epsilon = 1e-12
    );
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(fit.hat_diagonal.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.hat_diagonal.as_ref().unwrap().n_cols(), 3);
    assert_wald_likelihood_state(&fit, &counts);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_rows(), 2);
    assert_eq!(beta_covariance.n_cols(), 1);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[0].sqrt(),
        fit.beta_se.as_ref().unwrap().row(0).unwrap()[0],
        epsilon = 1e-12
    );
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), 2);
    assert_eq!(fit.max_cooks.as_ref().unwrap().len(), 2);
    assert!(results.rows[0].max_cooks.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_uses_normalization_factors_for_offsets() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        3,
        vec![10, 20, 40],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 2.0, 4.0]).unwrap();
    let design =
        DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], Some(vec!["Intercept".into()]))
            .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![100.0, 100.0, 100.0])
        .normalization_factors(normalization_factors.clone())
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 0)
        .unwrap();

    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_eq!(fit.size_factors, vec![100.0, 100.0, 100.0]);
    assert_relative_eq!(fit.base_mean[0], 10.0, epsilon = 1e-12);
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[0],
        10.0_f64.log2(),
        epsilon = 1e-12
    );
    for (actual, expected) in fit
        .mu
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .zip([10.0, 20.0, 40.0])
    {
        assert_relative_eq!(*actual, expected, epsilon = 1e-12);
    }
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_wald_likelihood_state(&fit, &counts);
    assert_relative_eq!(results.rows[0].base_mean, 10.0, epsilon = 1e-12);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        10.0_f64.log2(),
        epsilon = 1e-12
    );
}

#[test]
fn fixed_dispersion_wald_pipeline_uses_irls_for_two_group_design() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), 2);
    assert!(fit.beta_converged.as_ref().unwrap()[0]);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_rows(), 1);
    assert_eq!(beta_covariance.n_cols(), 4);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[3].sqrt(),
        fit.beta_se.as_ref().unwrap().row(0).unwrap()[1],
        epsilon = 1e-8
    );
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[1],
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_wald_likelihood_state(&fit, &counts);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(results.rows[0].pvalue, fit.wald.as_ref().unwrap().pvalue[0]);
    assert_eq!(results.rows[0].dispersion, Some(0.05));
    assert_eq!(fit.cooks.as_ref().unwrap().n_cols(), 4);
    assert_eq!(results.rows[0].max_cooks, None);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_use_qr_irls_solver() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            solver: IrlsSolver::Qr,
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    assert!(fit.beta_converged.as_ref().unwrap()[0]);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_eq!(results.rows[0].pvalue, fit.wald.as_ref().unwrap().pvalue[0]);
}

#[test]
fn fixed_dispersion_wald_pipeline_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();
    let rank_deficient = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald(&counts, &design, &[], 0)
            .is_err()
    );
    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
            .is_err()
    );
    assert!(
        DeseqBuilder::new()
            .fit_fixed_dispersion_wald(&counts, &rank_deficient, &[0.1], 1)
            .is_err()
    );
}
