#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn bounded_log_alpha_proposal_keeps_unclamped_step() {
        let (proposal, effective_step) =
            bounded_log_alpha_proposal(0.0, 2.0, 0.5, -30.0, 10.0).unwrap();

        assert_relative_eq!(proposal, 1.0, epsilon = 1e-12);
        assert_relative_eq!(effective_step, 0.5, epsilon = 1e-12);
    }

    #[test]
    fn bounded_log_alpha_proposal_reports_clamped_step() {
        let (proposal, effective_step) =
            bounded_log_alpha_proposal(9.5, 2.0, 1.0, -30.0, 10.0).unwrap();

        assert_relative_eq!(proposal, 10.0, epsilon = 1e-12);
        assert_relative_eq!(effective_step, 0.25, epsilon = 1e-12);
    }

    #[test]
    fn bounded_log_alpha_proposal_rejects_no_movement_at_bound() {
        assert!(bounded_log_alpha_proposal(10.0, 2.0, 1.0, -30.0, 10.0).is_none());
    }

    #[test]
    fn bounded_log_alpha_proposal_rejects_overflowed_step() {
        assert!(bounded_log_alpha_proposal(0.0, f64::MAX, 2.0, -30.0, 10.0).is_none());
    }

    #[test]
    fn bounded_log_alpha_proposal_rejects_overflowed_effective_step() {
        assert!(
            bounded_log_alpha_proposal(-f64::MAX, f64::MIN_POSITIVE, 1.0, -30.0, 10.0).is_none()
        );
    }

    #[test]
    fn line_search_armijo_bound_rejects_nonfinite_arithmetic() {
        let err = checked_line_search_armijo_bound(0.0, 1.0, 1.0, f64::MAX).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "dispersion line-search Armijo scale" && index == Some(0)
        ));
    }

    #[test]
    fn line_search_armijo_bound_matches_finite_formula() {
        let observed = checked_line_search_armijo_bound(-10.0, 0.5, 1.0e-4, 2.0).unwrap();
        let expected = 10.0 - 0.5 * 1.0e-4 * 4.0;

        assert_relative_eq!(observed, expected, epsilon = 1e-15);
    }

    #[test]
    fn dispersion_grid_linspace_rejects_nonfinite_arithmetic() {
        let endpoint_err = linspace(f64::INFINITY, 1.0, 3).unwrap_err();
        assert!(matches!(
            endpoint_err,
            DeseqError::NonFiniteValue { context, .. } if context == "dispersion grid endpoint"
        ));

        let span_err = linspace(-f64::MAX, f64::MAX, 3).unwrap_err();
        assert!(matches!(
            span_err,
            DeseqError::NonFiniteValue { context, .. } if context == "dispersion grid span"
        ));

        assert_eq!(linspace(2.0, 2.0, 3).unwrap(), vec![2.0, 2.0, 2.0]);
    }

    #[test]
    fn checked_div_rejects_nonfinite_dispersion_summaries() {
        let zero_err = checked_div(1.0, 0.0, 7, "test dispersion division").unwrap_err();
        assert!(matches!(
            zero_err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test dispersion division" && index == Some(7)
        ));

        let overflow_err =
            checked_div(f64::MAX, f64::MIN_POSITIVE, 3, "test dispersion division").unwrap_err();
        assert!(matches!(
            overflow_err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test dispersion division" && index == Some(3)
        ));
    }

    #[test]
    fn checked_log_alpha_first_derivative_uses_reduced_product() {
        let value =
            checked_log_alpha_first_derivative(0.5, f64::MAX / 2.0, "test derivative").unwrap();

        assert_eq!(value, f64::MAX / 4.0);
    }

    #[test]
    fn checked_log_alpha_second_derivative_keeps_cancelling_large_terms() {
        let value = checked_log_alpha_second_derivative(
            f64::MAX / 2.0,
            0.5,
            f64::MAX / 2.0,
            1.25,
            "test second derivative",
        )
        .unwrap();

        assert!((value - 1.25).abs() < 1e-12);
    }

    #[test]
    fn checked_cox_reid_log_alpha_second_derivative_rejects_alpha_square_overflow() {
        let err = checked_cox_reid_log_alpha_second_derivative(1.0, f64::MAX, 0.0).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "Cox-Reid log-alpha second derivative alpha square"
                    && index == Some(0)
        ));
    }

    #[test]
    fn mu_alpha_terms_keep_overflowed_product_finite() {
        let terms = mu_alpha_terms(f64::MAX / 2.0, 4.0, 0, "test mean-dispersion terms").unwrap();

        assert!(terms.log1p.is_finite());
        assert_eq!(terms.ratio, 1.0);
        assert_eq!(terms.inv_one_plus_squared, 0.0);
        assert!(terms.alpha_over_one_plus.is_finite());
        assert_eq!(terms.mu_squared_alpha_over_one_plus_squared, 0.25);
    }

    #[test]
    fn dispersion_derivatives_keep_overflowed_mu_alpha_products_finite() {
        let counts = [0, 10];
        let mu = [f64::MAX / 2.0, f64::MAX / 3.0];
        let log_alpha = 4.0_f64.ln();

        let derivative =
            dispersion_nb_log_likelihood_kernel_derivative(&counts, &mu, log_alpha).unwrap();
        let second_derivative =
            dispersion_nb_log_likelihood_kernel_second_derivative(&counts, &mu, log_alpha).unwrap();

        assert!(derivative.is_finite());
        assert!(second_derivative.is_finite());
    }

    #[test]
    fn dispersion_digamma_matches_r_positive_path() {
        let cases = [
            (5.845_272_899_427_574, 1.677_662_134_240_236_7),
            (6.845_272_899_427_574, 1.848_740_545_326_557_3),
            (15.845_272_899_427_574, 2.730_984_286_464_481_8),
            (100_005.845_272_899_43, 11.512_978_916_274_841),
            (0.125, -8.388_492_663_295_855),
            (0.5, -1.963_510_026_021_423),
            (1.0, -0.577_215_664_901_532_3),
            (8.75, 2.110_823_820_758_861_7),
            (9.0, 2.140_641_477_955_61),
            (100_000_000.0, 18.420_680_738_952_367),
            (1.0e18, 41.446_531_673_892_82),
        ];

        for (x, expected) in cases {
            assert_relative_eq!(
                dispersion_digamma(x),
                expected,
                epsilon = 2e-15,
                max_relative = 2e-15
            );
        }
    }

    #[test]
    fn dispersion_lgamma_matches_r_stirling_thresholds() {
        let cases = [
            (10.5, 13.940_625_219_403_763),
            (11.0, 15.104_412_573_075_52),
            (100.0, 359.134_205_369_575_4),
            (4_934_721.0, 71_118_238.355_323_93),
            (100_000_000.0, 1_742_068_066.103_834_9),
            (1.0e18, 4.044_653_167_389_282e19),
        ];

        for (x, expected) in cases {
            assert_relative_eq!(
                dispersion_lgamma(x),
                expected,
                epsilon = 1e-12,
                max_relative = 1e-15
            );
        }
    }

    #[test]
    fn cox_reid_weight_terms_keep_overflowed_mu_alpha_product_finite() {
        let terms = cox_reid_weight_terms(f64::MAX / 2.0, 4.0, 0).unwrap();

        assert_eq!(terms.weight, 0.25);
        assert_eq!(terms.d_weight, -0.0625);
        assert_eq!(terms.d2_weight, 0.03125);
    }

    #[test]
    fn cox_reid_weight_terms_reject_overflowed_square() {
        let err = cox_reid_weight_terms(f64::MAX / 2.0, f64::MIN_POSITIVE, 0).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "Cox-Reid working weight square" && index == Some(0)
        ));
    }

    #[test]
    fn dispersion_prior_log_density_rejects_overflowed_residual_square() {
        let prior = DispersionPrior::new(0.0, 1.0).unwrap();
        let err = dispersion_prior_log_density(f64::MAX, prior).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "dispersion prior log residual square" && index == Some(0)
        ));
    }

    #[test]
    fn dispersion_kernel_keeps_small_mu_alpha_terms_stable() {
        let counts = [0, 1, 2];
        let mu = [1.0e-6, 2.0e-6, 3.0e-6];
        let log_alpha = 1.0e-8_f64.ln();

        let kernel = dispersion_nb_log_likelihood_kernel(&counts, &mu, log_alpha).unwrap();
        let derivative =
            dispersion_nb_log_likelihood_kernel_derivative(&counts, &mu, log_alpha).unwrap();
        let second_derivative =
            dispersion_nb_log_likelihood_kernel_second_derivative(&counts, &mu, log_alpha).unwrap();

        assert!(kernel.is_finite());
        assert!(derivative.is_finite());
        assert!(second_derivative.is_finite());
    }
}
