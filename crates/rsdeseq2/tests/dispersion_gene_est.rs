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

fn centered_difference<F>(x: f64, f: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = 1e-5;
    (f(x + h) - f(x - h)) / (2.0 * h)
}

#[test]
fn linear_model_mu_projects_rows_onto_design_cells() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 4, vec![10.0, 30.0, 10.0, 30.0, 5.0, 5.0, 9.0, 9.0])
            .unwrap();
    let mu = linear_model_mu(&normalized, &two_group_design()).unwrap();

    assert_eq!(mu.n_rows(), 2);
    assert_relative_eq!(mu.row(0).unwrap()[0], 20.0, epsilon = 1e-12);
    assert_relative_eq!(mu.row(0).unwrap()[1], 20.0, epsilon = 1e-12);
    assert_relative_eq!(mu.row(0).unwrap()[2], 20.0, epsilon = 1e-12);
    assert_relative_eq!(mu.row(0).unwrap()[3], 20.0, epsilon = 1e-12);
    assert_relative_eq!(mu.row(1).unwrap()[0], 5.0, epsilon = 1e-12);
    assert_relative_eq!(mu.row(1).unwrap()[2], 9.0, epsilon = 1e-12);
}

#[test]
fn rough_dispersion_matches_deseq2_formula() {
    let normalized = RowMajorMatrix::from_row_major(1, 4, vec![10.0, 30.0, 10.0, 30.0]).unwrap();
    let rough = rough_dispersion_estimates(&normalized, &two_group_design()).unwrap();

    assert_relative_eq!(rough[0], 0.4, epsilon = 1e-12);
}

#[test]
fn moments_dispersion_matches_deseq2_formula() {
    let moments =
        moments_dispersion_estimates(&[20.0], &[400.0 / 3.0], &[1.0, 1.0, 1.0, 1.0]).unwrap();

    assert_relative_eq!(moments[0], 17.0 / 60.0, epsilon = 1e-12);
}

#[test]
fn moments_dispersion_uses_normalization_factor_column_means() {
    let normalization_factors =
        RowMajorMatrix::from_row_major(2, 4, vec![1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 3.0, 1.0]).unwrap();
    let moments = moments_dispersion_estimates_with_normalization_factors(
        &[10.0, 20.0],
        &[30.0, 100.0],
        &normalization_factors,
        Some(&[false, false]),
    )
    .unwrap();

    assert_relative_eq!(moments[0], 29.0 / 120.0, epsilon = 1e-12);
    assert_relative_eq!(moments[1], 53.0 / 240.0, epsilon = 1e-12);
}

#[test]
fn initial_dispersion_uses_min_of_rough_and_moments_then_bounds() {
    let initial =
        initial_dispersion_estimates(&[0.4, 0.0, 100.0], &[0.25, 0.2, 50.0], 0.01, 10.0).unwrap();

    assert_relative_eq!(initial[0], 0.25, epsilon = 1e-12);
    assert_relative_eq!(initial[1], 0.01, epsilon = 1e-12);
    assert_relative_eq!(initial[2], 10.0, epsilon = 1e-12);
}

#[test]
fn builder_rejects_rank_deficient_linear_mu_design() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let rank_deficient = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    assert!(DeseqBuilder::new()
        .size_factors(vec![1.0; 4])
        .fit_gene_wise_dispersions_linear_mu(&counts, &rank_deficient)
        .is_err());
}

#[test]
fn dispersion_grid_prefers_minimum_for_perfect_fit() {
    let options = GeneWiseDispersionOptions::default();
    let (dispersion, evaluations) = fit_dispersion_grid_no_cr(
        &[10, 10, 20, 20],
        &[10.0, 10.0, 20.0, 20.0],
        0.1,
        options,
        4,
    )
    .unwrap();

    assert_relative_eq!(dispersion, options.min_disp, epsilon = 1e-15);
    assert_eq!(evaluations, options.grid_points * 2);
}

#[test]
fn dispersion_kernel_differs_from_full_likelihood_by_alpha_constant() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let log_alpha_a = 0.1_f64.ln();
    let log_alpha_b = 0.4_f64.ln();

    let full_a = nbinom_log_likelihood(&counts, &mu, log_alpha_a.exp()).unwrap();
    let full_b = nbinom_log_likelihood(&counts, &mu, log_alpha_b.exp()).unwrap();
    let kernel_a = dispersion_nb_log_likelihood_kernel(&counts, &mu, log_alpha_a).unwrap();
    let kernel_b = dispersion_nb_log_likelihood_kernel(&counts, &mu, log_alpha_b).unwrap();

    assert_relative_eq!(full_a - kernel_a, full_b - kernel_b, epsilon = 1e-12);
}

#[test]
fn cox_reid_adjustment_matches_intercept_hand_formula() {
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
    let mu = [10.0, 10.0, 20.0, 20.0];
    let alpha = 0.1_f64;
    let expected_weight_sum = mu
        .iter()
        .map(|value| (1.0 / value + alpha).recip())
        .sum::<f64>();

    let adjustment = cox_reid_adjustment(&design, &mu, alpha.ln()).unwrap();

    assert_relative_eq!(adjustment, -0.5 * expected_weight_sum.ln(), epsilon = 1e-12);
}

#[test]
fn cox_reid_log_posterior_adds_log_determinant_penalty() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let design = two_group_design();
    let log_alpha = 0.25_f64.ln();

    let without_cr = dispersion_log_posterior(&counts, &mu, None, log_alpha, false).unwrap();
    let with_cr = dispersion_log_posterior(&counts, &mu, Some(&design), log_alpha, true).unwrap();
    let adjustment = cox_reid_adjustment(&design, &mu, log_alpha).unwrap();

    assert_relative_eq!(with_cr - without_cr, adjustment, epsilon = 1e-12);
}

#[test]
fn likelihood_derivative_matches_finite_difference() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_nb_log_likelihood_kernel_derivative(&counts, &mu, log_alpha).unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_nb_log_likelihood_kernel(&counts, &mu, value).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn weighted_likelihood_kernel_matches_unweighted_for_unit_weights() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 1.0, 1.0, 1.0];
    let log_alpha = 0.25_f64.ln();

    let unweighted = dispersion_nb_log_likelihood_kernel(&counts, &mu, log_alpha).unwrap();
    let weighted =
        dispersion_nb_log_likelihood_kernel_weighted(&counts, &mu, log_alpha, Some(&weights))
            .unwrap();
    let unweighted_derivative =
        dispersion_nb_log_likelihood_kernel_derivative(&counts, &mu, log_alpha).unwrap();
    let weighted_derivative = dispersion_nb_log_likelihood_kernel_derivative_weighted(
        &counts,
        &mu,
        log_alpha,
        Some(&weights),
    )
    .unwrap();

    assert_relative_eq!(weighted, unweighted, epsilon = 1e-12);
    assert_relative_eq!(weighted_derivative, unweighted_derivative, epsilon = 1e-12);
}

#[test]
fn weighted_likelihood_derivative_matches_finite_difference() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_nb_log_likelihood_kernel_derivative_weighted(
        &counts,
        &mu,
        log_alpha,
        Some(&weights),
    )
    .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_nb_log_likelihood_kernel_weighted(&counts, &mu, value, Some(&weights)).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn likelihood_second_derivative_matches_finite_difference() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let log_alpha = 0.25_f64.ln();

    let analytic =
        dispersion_nb_log_likelihood_kernel_second_derivative(&counts, &mu, log_alpha).unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_nb_log_likelihood_kernel_derivative(&counts, &mu, value).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-4, max_relative = 1e-4);
}

#[test]
fn weighted_likelihood_second_derivative_matches_finite_difference() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
        &counts,
        &mu,
        log_alpha,
        Some(&weights),
    )
    .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_nb_log_likelihood_kernel_derivative_weighted(&counts, &mu, value, Some(&weights))
            .unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-4, max_relative = 1e-4);
}

#[test]
fn weighted_cox_reid_adjustment_matches_intercept_hand_formula() {
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
    let mu = [10.0, 10.0, 20.0, 20.0];
    let weights = [1.0, 0.5, 1.0, 0.25];
    let alpha = 0.1_f64;
    let expected_weight_sum = mu
        .iter()
        .copied()
        .map(|mu| (1.0 / mu + alpha).recip())
        .sum::<f64>();

    let adjustment =
        cox_reid_adjustment_weighted(&design, &mu, alpha.ln(), Some(&weights)).unwrap();

    assert_relative_eq!(adjustment, -0.5 * expected_weight_sum.ln(), epsilon = 1e-12);
}

#[test]
fn weighted_cox_reid_adjustment_uses_thresholded_design_subset() {
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
    let mu = [10.0, 10.0, 20.0, 20.0];
    let weights = [1.0, 0.005, 1.0, 0.0];
    let alpha = 0.1_f64;
    let expected_weight_sum = [mu[0], mu[2]]
        .iter()
        .copied()
        .map(|mu| (1.0 / mu + alpha).recip())
        .sum::<f64>();

    let adjustment =
        cox_reid_adjustment_weighted(&design, &mu, alpha.ln(), Some(&weights)).unwrap();

    assert_relative_eq!(adjustment, -0.5 * expected_weight_sum.ln(), epsilon = 1e-12);
}

#[test]
fn cox_reid_derivative_matches_finite_difference() {
    let design = two_group_design();
    let mu = [20.0, 20.0, 20.0, 20.0];
    let log_alpha = 0.25_f64.ln();

    let analytic = cox_reid_adjustment_derivative(&design, &mu, log_alpha).unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        cox_reid_adjustment(&design, &mu, value).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn cox_reid_second_derivative_matches_finite_difference() {
    let design = two_group_design();
    let mu = [20.0, 20.0, 20.0, 20.0];
    let log_alpha = 0.25_f64.ln();

    let analytic = cox_reid_adjustment_second_derivative(&design, &mu, log_alpha).unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        cox_reid_adjustment_derivative(&design, &mu, value).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-4, max_relative = 1e-4);
}

#[test]
fn weighted_cox_reid_second_derivative_matches_finite_difference() {
    let design = two_group_design();
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let log_alpha = 0.25_f64.ln();

    let analytic =
        cox_reid_adjustment_second_derivative_weighted(&design, &mu, log_alpha, Some(&weights))
            .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        cox_reid_adjustment_derivative_weighted(&design, &mu, value, Some(&weights)).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-4, max_relative = 1e-4);
}

#[test]
fn log_posterior_derivative_matches_finite_difference_with_cox_reid() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let design = two_group_design();
    let log_alpha = 0.25_f64.ln();

    let analytic =
        dispersion_log_posterior_derivative(&counts, &mu, Some(&design), log_alpha, true).unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_log_posterior(&counts, &mu, Some(&design), value, true).unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn weighted_posterior_derivative_matches_finite_difference() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let design = two_group_design();
    let prior = DispersionPrior::new(0.2_f64.ln(), 0.5).unwrap();
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_log_posterior_derivative_with_prior_and_weights(
        &counts,
        &mu,
        Some(&design),
        log_alpha,
        true,
        Some(prior),
        Some(&weights),
    )
    .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_log_posterior_with_prior_and_weights(
            &counts,
            &mu,
            Some(&design),
            value,
            true,
            Some(prior),
            Some(&weights),
        )
        .unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn dispersion_prior_kernel_matches_deseq2_formula() {
    let prior = DispersionPrior::new(0.5, 2.0).unwrap();
    let log_alpha = 0.25_f64.ln();
    let expected_kernel = -0.5 * (log_alpha - prior.log_mean).powi(2) / prior.variance;
    let expected_derivative = -(log_alpha - prior.log_mean) / prior.variance;
    let expected_second_derivative = -prior.variance.recip();

    assert_relative_eq!(
        dispersion_prior_log_density(log_alpha, prior).unwrap(),
        expected_kernel,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        dispersion_prior_derivative(log_alpha, prior).unwrap(),
        expected_derivative,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        dispersion_prior_second_derivative(log_alpha, prior).unwrap(),
        expected_second_derivative,
        epsilon = 1e-12
    );
}

#[test]
fn dispersion_prior_validates_finite_positive_parameters() {
    assert!(DispersionPrior::new(f64::NAN, 1.0).is_err());
    assert!(DispersionPrior::new(0.0, 0.0).is_err());
    assert!(DispersionPrior::new(0.0, -1.0).is_err());
}

#[test]
fn prior_log_posterior_adds_expected_penalty() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let prior = DispersionPrior::new(0.5, 2.0).unwrap();
    let log_alpha = 0.25_f64.ln();

    let without_prior = dispersion_log_posterior(&counts, &mu, None, log_alpha, false).unwrap();
    let with_prior =
        dispersion_log_posterior_with_prior(&counts, &mu, None, log_alpha, false, Some(prior))
            .unwrap();

    assert_relative_eq!(
        with_prior - without_prior,
        dispersion_prior_log_density(log_alpha, prior).unwrap(),
        epsilon = 1e-12
    );
}

#[test]
fn prior_derivative_matches_finite_difference_with_cox_reid() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let design = two_group_design();
    let prior = DispersionPrior::new(0.05_f64.ln(), 0.2).unwrap();
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_log_posterior_derivative_with_prior(
        &counts,
        &mu,
        Some(&design),
        log_alpha,
        true,
        Some(prior),
    )
    .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_log_posterior_with_prior(&counts, &mu, Some(&design), value, true, Some(prior))
            .unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-5, max_relative = 1e-5);
}

#[test]
fn posterior_second_derivative_matches_finite_difference_with_prior_and_weights() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let design = two_group_design();
    let prior = DispersionPrior::new(0.05_f64.ln(), 0.2).unwrap();
    let log_alpha = 0.25_f64.ln();

    let analytic = dispersion_log_posterior_second_derivative_with_prior_and_weights(
        &counts,
        &mu,
        Some(&design),
        log_alpha,
        true,
        Some(prior),
        Some(&weights),
    )
    .unwrap();
    let numeric = centered_difference(log_alpha, |value| {
        dispersion_log_posterior_derivative_with_prior_and_weights(
            &counts,
            &mu,
            Some(&design),
            value,
            true,
            Some(prior),
            Some(&weights),
        )
        .unwrap()
    });

    assert_relative_eq!(analytic, numeric, epsilon = 1e-4, max_relative = 1e-4);
}

#[test]
fn line_search_with_prior_moves_toward_prior_mean() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let options = GeneWiseDispersionOptions {
        use_cox_reid: false,
        disp_tol: 1e-10,
        ..GeneWiseDispersionOptions::default()
    };
    let prior = DispersionPrior::new(0.02_f64.ln(), 0.05).unwrap();

    let no_prior =
        fit_dispersion_line_search_no_cr(&counts, &mu, 0.5, options, counts.len()).unwrap();
    let with_prior = fit_dispersion_line_search_no_cr_with_prior(
        &counts,
        &mu,
        0.5,
        options,
        counts.len(),
        prior,
    )
    .unwrap();

    assert!(with_prior.last_lp >= with_prior.initial_lp);
    assert!(with_prior.last_d2lp.is_finite());
    assert!(
        (with_prior.log_alpha - prior.log_mean).abs() < (no_prior.log_alpha - prior.log_mean).abs()
    );
}

#[test]
fn weighted_line_search_and_grid_with_prior_return_finite_estimates() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let weights = [1.0, 0.25, 1.0, 0.5];
    let design = two_group_design();
    let prior = DispersionPrior::new(0.2_f64.ln(), 0.5).unwrap();
    let options = GeneWiseDispersionOptions::default();

    let line_search =
        fit_dispersion_line_search_with_prior_and_weights(WeightedDispersionFitInput {
            counts: &counts,
            mu: &mu,
            design: &design,
            initial_dispersion: 0.4,
            options,
            n_samples: counts.len(),
            prior,
            weights: &weights,
        })
        .unwrap();
    let (grid, _) = fit_dispersion_grid_with_prior_and_weights(WeightedDispersionFitInput {
        counts: &counts,
        mu: &mu,
        design: &design,
        initial_dispersion: line_search.dispersion,
        options,
        n_samples: counts.len(),
        prior,
        weights: &weights,
    })
    .unwrap();

    assert!(line_search.dispersion.is_finite());
    assert!(line_search.dispersion >= options.min_disp);
    assert!(grid.is_finite());
    assert!(grid >= options.min_disp);
}

#[test]
fn dispersion_grid_with_cox_reid_returns_bounded_estimate() {
    let options = GeneWiseDispersionOptions::default();
    let (dispersion, evaluations) = fit_dispersion_grid(
        &[10, 30, 10, 30],
        &[20.0, 20.0, 20.0, 20.0],
        &two_group_design(),
        0.1,
        options,
        4,
    )
    .unwrap();

    assert!(dispersion >= options.min_disp);
    assert!(dispersion <= 10.0);
    assert_eq!(evaluations, options.grid_points * 2);
}

#[test]
fn dispersion_line_search_improves_objective() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let design = two_group_design();
    let options = GeneWiseDispersionOptions::default();

    let output =
        fit_dispersion_line_search(&counts, &mu, &design, 0.1, options, counts.len()).unwrap();

    assert!(output.iter > 0);
    assert!(output.iter <= options.maxit);
    assert!(output.iter_accept > 0);
    assert!(output.last_lp >= output.initial_lp);
    assert!(output.dispersion >= options.min_disp);
    assert!(output.dispersion <= 10.0);
}

#[test]
fn dispersion_line_search_no_cr_matches_grid_scale() {
    let counts = [10, 30, 10, 30];
    let mu = [20.0, 20.0, 20.0, 20.0];
    let options = GeneWiseDispersionOptions {
        use_cox_reid: false,
        ..GeneWiseDispersionOptions::default()
    };

    let line_search =
        fit_dispersion_line_search_no_cr(&counts, &mu, 0.1, options, counts.len()).unwrap();
    let (grid, _) = fit_dispersion_grid_no_cr(&counts, &mu, 0.1, options, counts.len()).unwrap();

    assert!(line_search.dispersion >= options.min_disp);
    assert!(line_search.dispersion <= 10.0);
    assert!((line_search.dispersion.ln() - grid.ln()).abs() < 1.0);
}

#[test]
fn linear_mu_gene_wise_dispersion_estimator_expands_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        4,
        vec![
            0, 0, 0, 0, //
            10, 10, 20, 20, //
            10, 30, 10, 30,
        ],
        Some(vec!["zero".into(), "fit".into(), "variable".into()]),
        None,
    )
    .unwrap();
    let size_factors = vec![1.0, 1.0, 1.0, 1.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = counts.all_zero_flags();

    let output = estimate_gene_wise_dispersions_linear_mu(
        GeneWiseDispersionInput {
            counts: &counts,
            design: &two_group_design(),
            size_factors: &size_factors,
            normalization_factors: None,
            normalized_counts: &normalized,
            base_mean: &base_mean,
            base_var: &base_var,
            all_zero: &all_zero,
            observation_weights: None,
        },
        GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            ..GeneWiseDispersionOptions::default()
        },
    )
    .unwrap();

    assert!(output.disp_gene_est[0].is_nan());
    assert!(!output.converged[0]);
    assert!(output.mu.row(0).unwrap()[0].is_nan());
    assert_relative_eq!(output.disp_gene_est[1], 1e-8, epsilon = 1e-15);
    assert!(output.converged[1]);
    assert!(output.disp_gene_est[2] > 0.01);
    assert!(output.disp_gene_est[2] < 10.0);
    assert_eq!(output.disp_iter[2], 40);
}

#[test]
fn linear_mu_gene_wise_dispersion_estimator_uses_normalization_factor_offsets() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 20, 20, 40]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(1, 4, vec![1.0, 2.0, 1.0, 2.0]).unwrap();

    let fit = DeseqBuilder::new()
        .normalization_factors(normalization_factors.clone())
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_linear_mu(&counts, &two_group_design())
        .unwrap();

    assert_eq!(fit.size_factors, vec![1.0, 1.0, 1.0, 1.0]);
    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_relative_eq!(fit.base_mean[0], 15.0, epsilon = 1e-12);
    for (actual, expected) in fit
        .mu
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .zip([10.0, 20.0, 20.0, 40.0])
    {
        assert_relative_eq!(*actual, expected, epsilon = 1e-12);
    }
    assert_relative_eq!(
        fit.disp_gene_est.as_ref().unwrap()[0],
        GeneWiseDispersionOptions::default().min_disp,
        epsilon = 1e-15
    );
}

#[test]
fn glm_mu_gene_wise_dispersion_estimator_refits_means_and_expands_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            0, 0, 0, 0, //
            10, 10, 20, 20, //
            10, 30, 10, 30,
        ],
    )
    .unwrap();
    let size_factors = vec![1.0, 1.0, 1.0, 1.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = counts.all_zero_flags();

    let output = estimate_gene_wise_dispersions_glm_mu(
        GeneWiseDispersionInput {
            counts: &counts,
            design: &two_group_design(),
            size_factors: &size_factors,
            normalization_factors: None,
            normalized_counts: &normalized,
            base_mean: &base_mean,
            base_var: &base_var,
            all_zero: &all_zero,
            observation_weights: None,
        },
        GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        },
        IrlsOptions::default(),
    )
    .unwrap();

    assert!(output.disp_gene_est[0].is_nan());
    assert!(!output.converged[0]);
    assert!(output.mu.row(0).unwrap()[0].is_nan());
    assert_relative_eq!(output.mu.row(1).unwrap()[0], 10.0, epsilon = 1e-4);
    assert_relative_eq!(output.mu.row(1).unwrap()[2], 20.0, epsilon = 1e-4);
    assert_relative_eq!(
        output.disp_gene_est[1],
        GeneWiseDispersionOptions::default().min_disp,
        epsilon = 1e-15
    );
    assert!(output.disp_gene_est[2] >= GeneWiseDispersionOptions::default().min_disp);
    assert!(output.disp_gene_est[2] <= 10.0);
    assert!(output.disp_iter[1] > 0);
}

#[test]
fn builder_attaches_linear_mu_gene_wise_dispersion_state() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 10, 30, 10, 30]).unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_linear_mu(&counts, &two_group_design())
        .unwrap();

    let dispersions = fit.disp_gene_est.as_ref().unwrap();
    let gene_iterations = fit.disp_gene_iter.as_ref().unwrap();
    assert_eq!(dispersions.len(), 2);
    assert_eq!(gene_iterations.len(), 2);
    assert_relative_eq!(dispersions[0], 1e-8, epsilon = 1e-15);
    assert!(dispersions[1] > 0.01);
    assert!(gene_iterations.iter().all(|iterations| *iterations > 0));
    assert_eq!(fit.dispersion, None);
    assert_eq!(fit.disp_iter, None);
    assert_eq!(
        fit.dispersion_converged.as_ref().unwrap(),
        &vec![true, true]
    );
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), 4);
}

#[test]
fn builder_attaches_glm_mu_gene_wise_dispersion_state() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 10, 30, 10, 30]).unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_glm_mu(&counts, &two_group_design())
        .unwrap();

    let dispersions = fit.disp_gene_est.as_ref().unwrap();
    let gene_iterations = fit.disp_gene_iter.as_ref().unwrap();
    assert_eq!(dispersions.len(), 2);
    assert_eq!(gene_iterations.len(), 2);
    assert_relative_eq!(dispersions[0], 1e-8, epsilon = 1e-15);
    assert!(gene_iterations.iter().all(|iterations| *iterations > 0));
    assert_eq!(fit.dispersion, None);
    assert_eq!(fit.disp_iter, None);
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), 4);
    assert_relative_eq!(
        fit.mu.as_ref().unwrap().row(0).unwrap()[0],
        10.0,
        epsilon = 1e-4
    );
    assert_eq!(fit.dispersion_converged.as_ref().unwrap().len(), 2);
}

#[test]
fn weighted_glm_mu_gene_wise_dispersion_matches_deseq2_fitidx_weight_rows() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![
            10, 12, 20, 24, //
            0, 0, 5, 7, //
            100, 80, 90, 120, //
            3, 6, 9, 12,
        ],
    )
    .unwrap();
    let observation_weights = RowMajorMatrix::from_row_major(
        4,
        4,
        vec![
            1.0, 2.0, 1.0, 2.0, //
            1.0, 1.0, 0.0, 0.0, //
            2.0, 1.0, 2.0, 1.0, //
            1.0, 0.5, 1.0, 0.5,
        ],
    )
    .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![
            0.645497224367903,
            0.829777303051304,
            1.29099444873581,
            1.54919333848297,
        ])
        .observation_weights(observation_weights)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_glm_mu(&counts, &two_group_design())
        .unwrap();

    assert_eq!(
        fit.weights_fail.as_ref().unwrap(),
        &vec![false, true, false, false]
    );
    assert_eq!(fit.all_zero, vec![false, true, false, false]);
    let dispersions = fit.disp_gene_est.as_ref().unwrap();
    assert_relative_eq!(dispersions[0], 1.0e-8, epsilon = 1e-14);
    assert!(dispersions[1].is_nan());
    assert_relative_eq!(dispersions[2], 0.0279348053798306, epsilon = 1e-12);
    assert_relative_eq!(dispersions[3], 1.0e-8, epsilon = 1e-14);
}

#[test]
fn weighted_glm_mu_gene_wise_dispersion_matches_deseq2_weighted_cox_reid() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![
            10, 12, 20, 24, //
            0, 0, 5, 7, //
            100, 80, 90, 120, //
            3, 6, 9, 12,
        ],
    )
    .unwrap();
    let observation_weights = RowMajorMatrix::from_row_major(
        4,
        4,
        vec![
            1.0, 2.0, 1.0, 2.0, //
            1.0, 1.0, 0.0, 0.0, //
            2.0, 1.0, 2.0, 1.0, //
            1.0, 0.5, 1.0, 0.5,
        ],
    )
    .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![
            0.645497224367903,
            0.829777303051304,
            1.29099444873581,
            1.54919333848297,
        ])
        .observation_weights(observation_weights)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_glm_mu(&counts, &two_group_design())
        .unwrap();

    assert_eq!(
        fit.weights_fail.as_ref().unwrap(),
        &vec![false, true, false, false]
    );
    assert_eq!(fit.all_zero, vec![false, true, false, false]);
    let dispersions = fit.disp_gene_est.as_ref().unwrap();
    assert_relative_eq!(dispersions[0], 1.00000006870857e-8, epsilon = 1e-14);
    assert!(dispersions[1].is_nan());
    assert_relative_eq!(dispersions[2], 0.101575541631976, epsilon = 1e-12);
    assert_relative_eq!(dispersions[3], 1.0e-8, epsilon = 1e-14);
}

#[test]
fn dispersion_estimator_validates_residual_degrees_of_freedom() {
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap();
    let saturated = DesignMatrix::from_row_major(2, 2, vec![1.0, 0.0, 0.0, 1.0], None).unwrap();

    assert!(rough_dispersion_estimates(&normalized, &saturated).is_err());
}

#[test]
fn gene_wise_dispersion_options_reject_zero_niter() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let size_factors = vec![1.0, 1.0, 1.0, 1.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = counts.all_zero_flags();

    assert!(estimate_gene_wise_dispersions_glm_mu(
        GeneWiseDispersionInput {
            counts: &counts,
            design: &two_group_design(),
            size_factors: &size_factors,
            normalization_factors: None,
            normalized_counts: &normalized,
            base_mean: &base_mean,
            base_var: &base_var,
            all_zero: &all_zero,
            observation_weights: None,
        },
        GeneWiseDispersionOptions {
            niter: 0,
            ..GeneWiseDispersionOptions::default()
        },
        IrlsOptions::default(),
    )
    .is_err());
}
