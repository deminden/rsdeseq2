// Private numerical regression tests.

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde::Deserialize;

    #[test]
    fn log2_output_scaling_rejects_nonfinite_values() {
        assert!(matches!(
            checked_log2_scale(f64::MAX, 0, "test beta scale"),
            Err(DeseqError::NonFiniteValue { context, index, .. })
                if context == "test beta scale" && index == Some(0)
        ));
        assert!(matches!(
            checked_log2_covariance(f64::MAX, 1, "test covariance scale"),
            Err(DeseqError::NonFiniteValue { context, index, .. })
                if context == "test covariance scale" && index == Some(1)
        ));
    }

    #[test]
    fn optim_start_keeps_finite_beta_on_bounds() {
        let x = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let counts = [1_u32, 2_u32];
        let nf = [1.0, 1.0];

        let start = optim_start_beta_log2(&[30.0, -30.0], &x, &counts, &nf, 30.0).unwrap();

        assert_eq!(start, vec![30.0, -30.0]);
    }

    #[test]
    fn optim_start_preserves_backup_qr_natural_log_values() {
        let x = DMatrix::from_element(2, 1, 1.0);
        let counts = [9, 9];
        let nf = [1.0, 1.0];
        let start = optim_start_beta_log2(&[f64::NAN], &x, &counts, &nf, 30.0).unwrap();

        assert_relative_eq!(start[0], 9.1_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn optim_fallback_control_matches_r_optim_lbfgsb_defaults() {
        let options = IrlsOptions::default();
        assert_eq!(options.max_beta_abs, 30.0);
        let control = optim_fallback_control(3, 100, 0.0);

        assert_eq!(control.maxit, 100);
        assert_eq!(control.factr, 1.0e7);
        assert_eq!(control.pgtol, 0.0);
        assert_eq!(control.lmm, 5);
        assert_eq!(control.fnscale, 1.0);
        assert_eq!(control.parscale, vec![1.0; 3]);
        assert_eq!(control.ndeps, vec![1.0e-3; 3]);
    }

    #[test]
    fn r_optim_candidate_accepts_objective_gap_within_factr_band() {
        let analytic = 427.391043312651;
        let tolerance = 10.0 * 1.0e7 * f64::EPSILON * analytic;

        assert!(r_optim_candidate_within_analytic_tolerance(
            analytic + tolerance,
            analytic,
            1.0e7,
        ));
    }

    #[test]
    fn r_optim_candidate_rejects_material_objective_gap() {
        let analytic = 45.11460017889003;
        let tolerance = 10.0 * 1.0e7 * f64::EPSILON * analytic;

        assert!(!r_optim_candidate_within_analytic_tolerance(
            analytic + 1.01 * tolerance,
            analytic,
            1.0e7,
        ));
        assert!(!r_optim_candidate_within_analytic_tolerance(
            f64::INFINITY,
            analytic,
            1.0e7,
        ));
    }

    #[test]
    fn r_optim_candidate_disagreement_is_scaled_by_numeric_resolution() {
        let ndeps = [1.0e-3, 1.0e-3];
        let parscale = [1.0, 2.0];

        assert!(!r_optim_candidates_materially_disagree(
            &[1.0, 2.0],
            &[1.009, 2.019],
            &ndeps,
            &parscale,
        ));
        assert!(r_optim_candidates_materially_disagree(
            &[1.0, 2.0],
            &[1.010_001, 2.0],
            &ndeps,
            &parscale,
        ));
    }

    #[test]
    fn beta_objective_uses_log2_mu_and_nb_kernel_without_gamma_constants() {
        let x = DMatrix::from_row_slice(1, 1, &[1.0]);
        let counts = [3_u32];
        let nf = [2.0];
        let ridge_lambda = [0.0];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.5,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[1.0]).unwrap();
        let size = 2.0_f64;
        let mu = 4.0_f64;
        let expected = (3.0 + size) * (size + mu).ln() - size * size.ln() - 3.0 * mu.ln();
        let expected_gradient = std::f64::consts::LN_2 * (mu - 3.0) / (1.0 + 0.5 * mu);
        let expected_hessian =
            std::f64::consts::LN_2.powi(2) * mu * (1.0 + 0.5 * 3.0) / (1.0 + 0.5 * mu).powi(2);

        assert_relative_eq!(objective, expected, epsilon = 1e-12);
        assert_relative_eq!(gradient[0], expected_gradient, epsilon = 1e-12);
        assert_relative_eq!(hessian[(0, 0)], expected_hessian, epsilon = 1e-12);
    }

    #[test]
    fn beta_objective_converts_natural_log_ridge_to_log2_penalty() {
        let x = DMatrix::from_row_slice(1, 1, &[0.0]);
        let counts = [0_u32];
        let nf = [1.0];
        let natural_log_ridge = 2.5;
        let ridge_lambda = [natural_log_ridge];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.5,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let beta = [3.0];
        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &beta).unwrap();
        let base = beta_log2_nb_objective_term(0, 1.0, 0.5).unwrap();
        let ridge_log2 = natural_log_ridge * std::f64::consts::LN_2.powi(2);

        assert_relative_eq!(
            objective,
            base + 0.5 * ridge_log2 * beta[0] * beta[0],
            epsilon = 1e-12
        );
        assert_relative_eq!(gradient[0], ridge_log2 * beta[0], epsilon = 1e-12);
        assert_relative_eq!(hessian[(0, 0)], ridge_log2, epsilon = 1e-12);
    }

    #[test]
    fn beta_numeric_gradient_uses_r_optim_bounded_difference() {
        let x = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let counts = [1_u32, 8_u32];
        let nf = [1.0, 1.5];
        let ridge_lambda = [0.25];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.3,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };
        let beta = [30.0];
        let gradient = beta_log2_numeric_gradient(&input, &beta, -30.0, 30.0).unwrap();
        let forward_value = beta_log2_objective(&input, &[30.0]).unwrap();
        let backward_value = beta_log2_objective(&input, &[29.999]).unwrap();

        assert_relative_eq!(
            gradient[0],
            (forward_value - backward_value) / 0.001,
            epsilon = 1e-5
        );
    }

    #[test]
    fn optim_fallback_replays_portable_objective_only_r_optim_gtex_cases() {
        let fixture: HardROptimFixture = serde_json::from_str(include_str!(
            "fixtures/gtex_optim_reference.json"
        ))
        .unwrap();
        assert_eq!(fixture.contrasts.len(), 1);
        let contrast = &fixture.contrasts[0];
        assert_eq!(contrast.contrast, "GTEx_kidney_blocked");
        assert_eq!(contrast.cases.len(), 6);
        for case in &contrast.cases {
            assert_eq!(case.case_kind, "GTEx_sparse_effect_expression_extreme");
            assert_eq!(case.initial_par.len(), contrast.coefficients);
            assert_eq!(case.result.par.len(), contrast.coefficients);
            assert_eq!(case.lower, vec![-30.0; contrast.coefficients]);
            assert_eq!(case.upper, vec![30.0; contrast.coefficients]);
            assert_eq!(case.control.maxit, 100);
            assert_eq!(case.control.factr, 1.0e7);
            assert_eq!(case.control.pgtol, 0.0);
            assert_eq!(case.control.lmm, 5);
            assert_eq!(case.control.ndeps, vec![1.0e-3; contrast.coefficients]);

            let design_values = contrast
                .design
                .iter()
                .flat_map(|row| row.iter().copied())
                .collect::<Vec<_>>();
            let x =
                DMatrix::from_row_slice(contrast.samples, contrast.coefficients, &design_values);
            let counts = case
                .counts
                .iter()
                .copied()
                .map(|value| {
                    assert!(value.fract() == 0.0 && value >= 0.0 && value <= f64::from(u32::MAX));
                    value as u32
                })
                .collect::<Vec<_>>();
            let ridge_lambda = vec![0.0; contrast.coefficients];
            let input = BetaOptimInput {
                x: &x,
                counts: &counts,
                nf: &contrast.size_factors,
                dispersion: case.dispersion,
                weights: None,
                ridge_lambda: &ridge_lambda,
            };

            let output = minimize_beta_log2_lbfgsb(
                input,
                &case.initial_par,
                case.lower[0],
                case.upper[0],
                case.control.maxit,
                case.control.pgtol,
                case.r_optim_compat,
            )
            .unwrap();

            assert!(output.converged, "{}", case.fixture);
            assert!(
                output.iterations.abs_diff(case.result.counts.gradient) <= 4,
                "{} optimizer count drift: actual={} expected={}",
                case.fixture,
                output.iterations,
                case.result.counts.gradient
            );
            assert_relative_eq!(output.value, case.result.value, epsilon = 5e-6);
            for (actual, expected) in output.parameters.iter().zip(case.result.par.iter()) {
                assert_relative_eq!(*actual, *expected, epsilon = 2e-3);
            }
        }
    }

    #[test]
    fn bounded_beta_numeric_helpers_reject_nonfinite_accumulation() {
        assert_relative_eq!(
            projected_gradient_norm(&[0.0, 0.0], &[f64::MAX / 2.0, f64::MAX / 2.0], -30.0, 30.0)
                .unwrap(),
            f64::MAX / 2.0 * 2.0_f64.sqrt(),
            epsilon = 1e292
        );
        assert_eq!(
            projected_gradient_norm(&[0.0, 0.0], &[f64::MAX, f64::MAX], -30.0, 30.0),
            None
        );
        assert_eq!(checked_div2(1.0, 0.0), None);
        assert_eq!(checked_div2(f64::NAN, 1.0), None);
    }

    #[test]
    fn irls_deviance_convergence_stat_rejects_nonfinite_arithmetic() {
        assert_eq!(irls_deviance_convergence_stat(f64::MAX, -f64::MAX), None);
        assert_eq!(irls_deviance_convergence_stat(f64::INFINITY, 1.0), None);
        assert_eq!(irls_deviance_convergence_stat(0.0, 0.0), Some(0.0));
    }

    #[test]
    fn weighted_least_squares_helpers_reject_nonfinite_accumulation() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        assert!(xtwx_and_xtwz(&x, &[2.0], &[1.0], &[0.0]).is_none());
        assert!(xtwx_without_ridge(&x, &[2.0]).is_none());

        let unit_x = DMatrix::from_row_slice(1, 1, &[1.0]);
        assert!(solve_weighted_least_squares_qr(&unit_x, &[4.0], &[f64::MAX], &[0.0]).is_none());
    }

    #[test]
    fn beta_objective_returns_penalty_for_nonfinite_accumulation() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let counts = [1_u32];
        let nf = [1.0];
        let ridge_lambda = [0.0];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.2,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[20.0]).unwrap();

        assert_eq!(objective, 1.0e300);
        assert_eq!(gradient, vec![0.0]);
        assert_relative_eq!(hessian[(0, 0)], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn beta_objective_returns_penalty_for_overflowed_ridge_term() {
        let x = DMatrix::from_row_slice(1, 1, &[1.0]);
        let counts = [1_u32];
        let nf = [1.0];
        let ridge_lambda = [f64::MAX];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.2,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[20.0]).unwrap();

        assert_eq!(objective, 1.0e300);
        assert_eq!(gradient, vec![0.0]);
        assert_relative_eq!(hessian[(0, 0)], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn fitted_mu_rejects_overflowed_linear_predictor() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let beta = DVector::from_vec(vec![2.0]);
        let nf = [1.0];

        let err = fitted_mu(&x, &beta, &nf, 0.5).unwrap_err();

        match err {
            DeseqError::NonFiniteValue { context, index, .. } => {
                assert_eq!(context, "IRLS linear predictor");
                assert_eq!(index, Some(0));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn fitted_mu_saturates_overflowed_exponential_mean() {
        let x = DMatrix::from_row_slice(1, 1, &[1000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [1.0];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert_eq!(mu, vec![f64::MAX]);
    }

    #[test]
    fn fitted_mu_avoids_intermediate_exponential_overflow_with_small_factor() {
        let x = DMatrix::from_row_slice(1, 1, &[1000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [f64::MIN_POSITIVE];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert!(mu[0].is_finite());
        assert!(mu[0] > 1.0e126);
    }

    #[test]
    fn fitted_mu_saturates_final_overflowed_mean_with_small_factor() {
        let x = DMatrix::from_row_slice(1, 1, &[2000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [f64::MIN_POSITIVE];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert_eq!(mu, vec![f64::MAX]);
    }

    #[test]
    fn fitted_mu_log2_rejects_overflowed_linear_predictor() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let nf = [1.0];

        let err = fitted_mu_log2_unfloored(&x, &[2.0], &nf).unwrap_err();

        match err {
            DeseqError::NonFiniteValue { context, index, .. } => {
                assert_eq!(context, "optim fallback linear predictor");
                assert_eq!(index, Some(0));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn working_weights_reject_nonfinite_arithmetic() {
        assert_eq!(working_weights(&[f64::MAX], 2.0, None).unwrap(), vec![0.5]);

        let err = working_weights(&[f64::MAX], f64::MIN_POSITIVE, Some(&[10.0])).unwrap_err();
        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "IRLS weighted working weight" && index == Some(0)
        ));
    }

    #[test]
    fn working_response_avoids_normalized_mean_overflow() {
        let z = working_response(&[f64::MAX], &[f64::MIN_POSITIVE], &[1.0]).unwrap();

        assert!(z[0].is_finite());
        assert!(z[0] > 1400.0);
    }

    #[test]
    fn working_response_rejects_nonfinite_arithmetic() {
        let err = working_response(&[f64::MIN_POSITIVE], &[1.0], &[f64::MAX]).unwrap_err();
        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "IRLS working response scaled residual" && index == Some(0)
        ));
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimFixture {
        contrasts: Vec<HardROptimContrast>,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimContrast {
        contrast: String,
        samples: usize,
        coefficients: usize,
        design: Vec<Vec<f64>>,
        size_factors: Vec<f64>,
        cases: Vec<HardROptimCase>,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimCase {
        r_optim_compat: bool,
        fixture: String,
        case_kind: String,
        dispersion: f64,
        counts: Vec<f64>,
        initial_par: Vec<f64>,
        lower: Vec<f64>,
        upper: Vec<f64>,
        control: HardROptimControl,
        result: HardROptimResult,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimControl {
        maxit: usize,
        ndeps: Vec<f64>,
        factr: f64,
        pgtol: f64,
        lmm: usize,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimResult {
        par: Vec<f64>,
        value: f64,
        counts: HardROptimCounts,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimCounts {
        gradient: usize,
    }
}
