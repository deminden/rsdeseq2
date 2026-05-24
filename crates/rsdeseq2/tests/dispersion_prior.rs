use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

#[test]
fn trigamma_matches_known_values() {
    assert_relative_eq!(
        trigamma(0.5).unwrap(),
        std::f64::consts::PI.powi(2) / 2.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        trigamma(2.0).unwrap(),
        std::f64::consts::PI.powi(2) / 6.0 - 1.0,
        epsilon = 1e-12
    );
}

#[test]
fn mad_squared_matches_r_default_constant() {
    let residuals = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let expected = 1.4826_f64.powi(2);

    assert_relative_eq!(mad_squared(&residuals).unwrap(), expected, epsilon = 1e-12);
}

#[test]
fn log_dispersion_residuals_use_above_min_rule_and_skip_missing_rows() {
    let (residuals, above_min) = log_dispersion_residuals_above_min(
        &[1e-7, 2e-6, f64::NAN, 0.4],
        &[1e-7, 1e-6, 1.0, 0.2],
        1e-8,
    )
    .unwrap();

    assert_eq!(above_min, vec![false, true, false, true]);
    assert_eq!(residuals.len(), 2);
    assert_relative_eq!(residuals[0], 2.0_f64.ln(), epsilon = 1e-12);
    assert_relative_eq!(residuals[1], 2.0_f64.ln(), epsilon = 1e-12);
}

#[test]
fn prior_variance_subtracts_sampling_variance_and_floors_at_quarter() {
    let residuals = [-2.0_f64, -1.0, 0.0, 1.0, 2.0];
    let disp_fit = vec![1.0; residuals.len()];
    let disp_gene_est = residuals
        .iter()
        .map(|residual| residual.exp())
        .collect::<Vec<_>>();

    let output =
        estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 10, 2).unwrap();
    let expected_var_log = 1.4826_f64.powi(2);
    let expected_sampling = std::f64::consts::PI.powi(2) / 6.0 - 1.0 - 0.25 - 1.0 / 9.0;

    assert_eq!(output.residual_degrees_of_freedom, 8);
    assert_relative_eq!(
        output.var_log_disp_estimates,
        expected_var_log,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        output.expected_log_dispersion_variance,
        expected_sampling,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        output.disp_prior_var,
        expected_var_log - expected_sampling,
        epsilon = 1e-12
    );
}

#[test]
fn estimate_dispersion_prior_wraps_prior_variance_stage() {
    let residuals = [-2.0_f64, -1.0, 0.0, 1.0, 2.0];
    let disp_fit = vec![1.0; residuals.len()];
    let disp_gene_est = residuals
        .iter()
        .map(|residual| residual.exp())
        .collect::<Vec<_>>();

    let stage_output = estimate_dispersion_prior(&disp_gene_est, &disp_fit, 1e-8, 10, 2).unwrap();
    let variance_output =
        estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 10, 2).unwrap();

    assert_eq!(stage_output, variance_output);
}

#[test]
fn prior_variance_uses_quarter_floor_when_sampling_variance_is_large() {
    let residuals = [-0.2_f64, 0.0, 0.2, 0.4, -0.4];
    let disp_fit = vec![1.0; residuals.len()];
    let disp_gene_est = residuals
        .iter()
        .map(|residual| residual.exp())
        .collect::<Vec<_>>();

    let output =
        estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 10, 2).unwrap();

    assert_relative_eq!(output.disp_prior_var, 0.25, epsilon = 1e-12);
}

#[test]
fn prior_variance_with_saturated_model_does_not_subtract_or_floor() {
    let residuals = [-0.2_f64, 0.0, 0.2, 0.4, -0.4];
    let disp_fit = vec![1.0; residuals.len()];
    let disp_gene_est = residuals
        .iter()
        .map(|residual| residual.exp())
        .collect::<Vec<_>>();

    let output = estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 5, 5).unwrap();
    let expected = (0.2 * 1.4826_f64).powi(2);

    assert_eq!(output.residual_degrees_of_freedom, 0);
    assert_relative_eq!(
        output.expected_log_dispersion_variance,
        0.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(output.disp_prior_var, expected, epsilon = 1e-12);
}

#[test]
fn low_df_prior_variance_uses_deterministic_histogram_match() {
    let residuals = [-1.0_f64, 0.0, 1.0, 2.0];
    let disp_fit = vec![1.0; residuals.len()];
    let disp_gene_est = residuals
        .iter()
        .map(|residual| residual.exp())
        .collect::<Vec<_>>();

    let output_a = estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 4, 2)
        .expect("df 1..=3 should use the deterministic low-df branch");
    let output_b = estimate_dispersion_prior_variance(&disp_gene_est, &disp_fit, 1e-8, 4, 2)
        .expect("low-df branch should be deterministic");

    assert_eq!(output_a.residual_degrees_of_freedom, 2);
    assert_relative_eq!(
        output_a.expected_log_dispersion_variance,
        trigamma(1.0).unwrap(),
        epsilon = 1e-12
    );
    assert!(output_a.disp_prior_var >= 0.25);
    assert!(output_a.disp_prior_var <= 8.0);
    assert_relative_eq!(
        output_a.disp_prior_var,
        output_b.disp_prior_var,
        epsilon = 0.0
    );
}

#[test]
fn low_df_prior_variance_validates_supported_df_range() {
    let err = estimate_low_df_prior_variance(&[0.0, 1.0], 4).unwrap_err();
    assert!(matches!(err, DeseqError::InvalidDimensions { .. }));
}

#[test]
fn prior_variance_errors_when_no_rows_above_min_disp() {
    let err = estimate_dispersion_prior_variance(&[1e-8, 2e-8], &[1e-8, 2e-8], 1e-8, 10, 2);

    assert!(err.is_err());
}
