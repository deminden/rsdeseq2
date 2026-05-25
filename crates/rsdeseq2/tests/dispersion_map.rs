use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn two_group_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "conditionB".into()]),
    )
    .unwrap()
}

fn one_gene_counts() -> CountMatrix {
    CountMatrix::from_row_major_u32(1, 4, vec![10, 30, 10, 30]).unwrap()
}

fn one_gene_mu() -> RowMajorMatrix<f64> {
    RowMajorMatrix::from_row_major(1, 4, vec![20.0, 20.0, 20.0, 20.0]).unwrap()
}

#[test]
fn map_dispersion_initial_value_matches_deseq2_rule() {
    assert_relative_eq!(
        map_dispersion_initial_value(0.5, 1.0).unwrap(),
        0.5,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        map_dispersion_initial_value(0.05, 1.0).unwrap(),
        1.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        map_dispersion_initial_value(f64::NAN, 1.0).unwrap(),
        1.0,
        epsilon = 1e-12
    );
}

#[test]
fn map_dispersion_outlier_matches_log_residual_rule() {
    assert!(map_dispersion_outlier(1.0, 0.1, 2.0, 0.01));
    assert!(!map_dispersion_outlier(0.12, 0.1, 2.0, 0.01));
    assert!(!map_dispersion_outlier(f64::NAN, 0.1, 2.0, 0.01));
}

#[test]
fn map_dispersion_outlier_rejects_overflowed_threshold_arithmetic() {
    assert!(!map_dispersion_outlier(10.0, 0.1, f64::MAX, f64::MAX));
}

#[test]
fn estimate_map_dispersions_runs_prior_aware_line_search() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            use_cox_reid: false,
            disp_tol: 1e-10,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert!(output.disp_iter[0] > 0);
    assert!(output.converged[0]);
    assert!(!output.disp_outlier[0]);
    assert_relative_eq!(output.disp_init[0], 0.5, epsilon = 1e-12);
    assert!(output.disp_map[0] < 0.5);
    assert_relative_eq!(output.dispersion[0], output.disp_map[0], epsilon = 1e-12);
}

#[test]
fn estimate_map_dispersions_accepts_observation_weights() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.25, 1.0, 0.5]).unwrap();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: Some(&weights),
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            disp_tol: 1e-10,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert!(output.disp_iter[0] > 0);
    assert!(output.disp_map[0].is_finite());
    assert!(output.dispersion[0].is_finite());
    assert_relative_eq!(output.disp_init[0], 0.5, epsilon = 1e-12);
}

#[test]
fn estimate_map_dispersions_unit_weights_match_unweighted() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0; 4]).unwrap();

    let weighted = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: Some(&weights),
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            disp_tol: 1e-10,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();
    let unweighted = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            disp_tol: 1e-10,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert_relative_eq!(
        weighted.disp_map[0],
        unweighted.disp_map[0],
        epsilon = 1e-12
    );
    assert_relative_eq!(
        weighted.dispersion[0],
        unweighted.dispersion[0],
        epsilon = 1e-12
    );
}

#[test]
fn estimate_map_dispersions_uses_grid_fallback_when_line_search_does_not_converge() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            use_cox_reid: false,
            maxit: 1,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert_eq!(output.disp_iter[0], 1);
    assert!(!output.converged[0]);
    assert!(output.disp_map[0].is_finite());
    assert!(output.disp_map[0] >= 1e-8);
}

#[test]
fn estimate_map_dispersions_replaces_outlier_final_dispersion_with_gene_estimate() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[10.0],
            disp_fit: &[0.1],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: 0.25,
            var_log_disp_estimates: 0.01,
        },
        MapDispersionOptions {
            use_cox_reid: false,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert!(output.disp_outlier[0]);
    assert_relative_eq!(output.dispersion[0], 10.0, epsilon = 1e-12);
    assert_ne!(output.dispersion[0], output.disp_map[0]);
}

#[test]
fn estimate_map_dispersions_expands_all_zero_rows_as_missing() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![0, 0, 0, 0, 10, 30, 10, 30]).unwrap();
    let design = two_group_design();
    let mu = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            f64::NAN,
            f64::NAN,
            f64::NAN,
            f64::NAN,
            20.0,
            20.0,
            20.0,
            20.0,
        ],
    )
    .unwrap();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[f64::NAN, 0.5],
            disp_fit: &[f64::NAN, 0.02],
            all_zero: &[true, false],
            observation_weights: None,
            disp_prior_var: 0.05,
            var_log_disp_estimates: 100.0,
        },
        MapDispersionOptions {
            use_cox_reid: false,
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert!(output.dispersion[0].is_nan());
    assert!(output.disp_map[0].is_nan());
    assert_eq!(output.disp_iter[0], 0);
    assert!(!output.converged[0]);
    assert!(!output.disp_outlier[0]);
    assert!(output.dispersion[1].is_finite());
}

#[test]
fn estimate_map_dispersions_validates_dimensions() {
    let counts = one_gene_counts();
    let design = two_group_design();
    let mu = one_gene_mu();

    let err = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5, 0.6],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: 0.05,
            var_log_disp_estimates: 1.0,
        },
        MapDispersionOptions::default(),
    );

    assert!(err.is_err());

    let bad_weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0; 3]).unwrap();
    let err = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &[0.5],
            disp_fit: &[0.02],
            all_zero: &[false],
            observation_weights: Some(&bad_weights),
            disp_prior_var: 0.05,
            var_log_disp_estimates: 1.0,
        },
        MapDispersionOptions::default(),
    );

    assert!(err.is_err());
}
