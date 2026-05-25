use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn trend_data() -> (Vec<f64>, Vec<f64>) {
    let means = vec![10.0, 20.0, 40.0, 80.0, 160.0, 320.0];
    let disps = means
        .iter()
        .map(|mean| 0.05 + 2.0 / mean)
        .collect::<Vec<_>>();
    (means, disps)
}

#[test]
fn parametric_trend_evaluates_deseq2_formula() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: 2.0,
    };

    assert_relative_eq!(trend.evaluate(20.0).unwrap(), 0.15, epsilon = 1e-12);
    assert_relative_eq!(trend.evaluate(80.0).unwrap(), 0.075, epsilon = 1e-12);
}

#[test]
fn parametric_trend_rejects_overflowed_tiny_mean_fit() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: f64::MAX,
    };

    let err = trend.evaluate(f64::MIN_POSITIVE).unwrap_err();
    assert!(err
        .to_string()
        .contains("parametric dispersion trend produced"));
    assert!(trend
        .evaluate_many_allow_missing(&[0.0, f64::MIN_POSITIVE])
        .is_err());
}

#[test]
fn dispersion_trends_report_generic_mean_validation_context() {
    let parametric = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: 2.0,
    };
    let mean = MeanDispersionTrend { mean_disp: 0.2 };
    let local = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 1,
        log_means: Vec::new(),
        log_disps: Vec::new(),
        weights: Vec::new(),
    };

    for err in [
        parametric.evaluate(0.0).unwrap_err(),
        mean.evaluate(0.0).unwrap_err(),
        local.evaluate(0.0).unwrap_err(),
    ] {
        assert!(err.to_string().contains("dispersion trend mean"));
    }
}

#[test]
fn dispersion_trends_validate_fit_state_before_allow_missing_evaluation() {
    let invalid_parametric = ParametricDispersionTrend {
        asympt_disp: f64::NAN,
        extra_pois: 2.0,
    };
    let invalid_mean = MeanDispersionTrend { mean_disp: 0.0 };
    let invalid_local = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 1,
        log_means: vec![10.0_f64.ln()],
        log_disps: vec![f64::NAN],
        weights: vec![1.0],
    };
    let missing_means = [0.0, f64::NAN];

    assert!(invalid_parametric
        .evaluate_many_allow_missing(&missing_means)
        .is_err());
    assert!(invalid_mean
        .evaluate_many_allow_missing(&missing_means)
        .is_err());
    assert!(invalid_local
        .evaluate_many_allow_missing(&missing_means)
        .is_err());
}

#[test]
fn parametric_trend_use_for_fit_matches_min_disp_rule() {
    let use_for_fit = parametric_trend_use_for_fit(
        &[0.0, 10.0, 20.0, f64::NAN, 40.0],
        &[1.0, 1e-7, 2e-6, 1.0, f64::NAN],
        1e-8,
    )
    .unwrap();

    assert_eq!(use_for_fit, vec![false, false, true, false, false]);
}

#[test]
fn trend_selection_helpers_reject_invalid_min_disp() {
    for min_disp in [0.0, f64::NAN, f64::INFINITY] {
        assert!(parametric_trend_use_for_fit(&[10.0], &[0.1], min_disp).is_err());
        assert!(local_trend_use_for_fit(&[10.0], &[0.1], min_disp).is_err());
        assert!(mean_trend_use_for_fit(&[10.0], &[0.1], min_disp).is_err());
        assert!(mean_trend_use_for_mean(&[10.0], &[0.1], min_disp).is_err());
    }
}

#[test]
fn trend_selection_helpers_report_helper_specific_dimension_contexts() {
    let parametric = parametric_trend_use_for_fit(&[10.0, 20.0], &[0.1], 1e-8).unwrap_err();
    let local = local_trend_use_for_fit(&[10.0, 20.0], &[0.1], 1e-8).unwrap_err();
    let mean_fit = mean_trend_use_for_fit(&[10.0, 20.0], &[0.1], 1e-8).unwrap_err();
    let mean_value = mean_trend_use_for_mean(&[10.0, 20.0], &[0.1], 1e-8).unwrap_err();

    assert!(parametric
        .to_string()
        .contains("parametric dispersion trend rows"));
    assert!(local.to_string().contains("local dispersion trend rows"));
    assert!(mean_fit
        .to_string()
        .contains("mean dispersion trend fit rows"));
    assert!(mean_value
        .to_string()
        .contains("mean dispersion trend rows"));
}

#[test]
fn fit_dispersion_trends_reject_invalid_min_disp_options() {
    let (means, disps) = trend_data();

    for min_disp in [0.0, f64::NAN, f64::INFINITY] {
        assert!(fit_parametric_dispersion_trend(
            &means,
            &disps,
            ParametricDispersionTrendOptions {
                min_disp,
                ..ParametricDispersionTrendOptions::default()
            },
        )
        .is_err());
        assert!(fit_local_dispersion_trend(
            &means,
            &disps,
            LocalDispersionTrendOptions {
                min_disp,
                ..LocalDispersionTrendOptions::default()
            },
        )
        .is_err());
        assert!(fit_mean_dispersion_trend(
            &means,
            &disps,
            MeanDispersionTrendOptions {
                min_disp,
                ..MeanDispersionTrendOptions::default()
            },
        )
        .is_err());
    }
}

#[test]
fn fit_dispersion_trends_reject_empty_inputs() {
    assert!(
        fit_parametric_dispersion_trend(&[], &[], ParametricDispersionTrendOptions::default(),)
            .is_err()
    );
    assert!(fit_local_dispersion_trend(&[], &[], LocalDispersionTrendOptions::default()).is_err());
    assert!(fit_mean_dispersion_trend(&[], &[], MeanDispersionTrendOptions::default()).is_err());
}

#[test]
fn fit_parametric_dispersion_trend_recovers_exact_curve() {
    let (means, disps) = trend_data();

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();

    assert_eq!(fit.genes_used, means.len());
    assert!(fit.converged);
    assert_relative_eq!(fit.trend.asympt_disp, 0.05, epsilon = 1e-9);
    assert_relative_eq!(fit.trend.extra_pois, 2.0, epsilon = 1e-8);
    for (mean, fitted) in means.iter().copied().zip(fit.disp_fit.iter().copied()) {
        assert_relative_eq!(fitted, 0.05 + 2.0 / mean, epsilon = 1e-9);
    }
}

#[test]
fn fit_parametric_dispersion_trend_robustly_ignores_large_residual_outlier() {
    let (mut means, mut disps) = trend_data();
    means.push(640.0);
    disps.push(1000.0);

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();

    assert_eq!(fit.genes_used, means.len());
    assert!(fit.converged);
    assert_relative_eq!(fit.trend.asympt_disp, 0.05, epsilon = 1e-8);
    assert_relative_eq!(fit.trend.extra_pois, 2.0, epsilon = 1e-7);
}

#[test]
fn fit_parametric_dispersion_trend_returns_nan_for_missing_fitted_rows() {
    let (mut means, mut disps) = trend_data();
    means.push(0.0);
    disps.push(f64::NAN);

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();

    assert!(fit.disp_fit.last().unwrap().is_nan());
}

#[test]
fn fit_parametric_dispersion_trend_errors_when_all_estimates_at_minimum() {
    let err = fit_parametric_dispersion_trend(
        &[10.0, 20.0, 30.0],
        &[1e-8, 2e-8, 3e-8],
        ParametricDispersionTrendOptions::default(),
    );

    assert!(err.is_err());
}

#[test]
fn fit_parametric_dispersion_trend_validates_dimensions() {
    let err = fit_parametric_dispersion_trend(
        &[10.0, 20.0],
        &[0.1],
        ParametricDispersionTrendOptions::default(),
    );

    assert!(err.is_err());
}

#[test]
fn fit_parametric_dispersion_trend_rejects_zero_iteration_limits() {
    let (means, disps) = trend_data();

    let outer_err = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions {
            max_outer_iter: 0,
            ..ParametricDispersionTrendOptions::default()
        },
    );
    let irls_err = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions {
            max_irls_iter: 0,
            ..ParametricDispersionTrendOptions::default()
        },
    );

    assert!(outer_err.is_err());
    assert!(irls_err.is_err());
}

#[test]
fn fit_parametric_dispersion_trend_rejects_invalid_numeric_controls() {
    let (means, disps) = trend_data();

    for options in [
        ParametricDispersionTrendOptions {
            min_residual: 0.0,
            ..ParametricDispersionTrendOptions::default()
        },
        ParametricDispersionTrendOptions {
            min_residual: 2.0,
            max_residual: 1.0,
            ..ParametricDispersionTrendOptions::default()
        },
        ParametricDispersionTrendOptions {
            coefficient_tol: f64::NAN,
            ..ParametricDispersionTrendOptions::default()
        },
        ParametricDispersionTrendOptions {
            glm_tol: 0.0,
            ..ParametricDispersionTrendOptions::default()
        },
    ] {
        assert!(fit_parametric_dispersion_trend(&means, &disps, options).is_err());
    }
}

#[test]
fn fit_dispersion_trend_dispatches_default_fit_types() {
    let (means, disps) = trend_data();

    let parametric = fit_dispersion_trend(&means, &disps, FitType::Parametric).unwrap();
    assert!(matches!(parametric, DispersionTrendFit::Parametric(_)));
    assert_eq!(parametric.disp_fit().len(), means.len());
    assert_eq!(parametric.use_for_fit(), vec![true; means.len()].as_slice());
    assert!(parametric.use_for_mean().is_none());
    assert_eq!(parametric.genes_used_for_mean(), None);
    assert_eq!(parametric.used_min_disp_floor(), None);
    assert_eq!(parametric.genes_used(), means.len());

    let local = fit_dispersion_trend(&means, &disps, FitType::Local).unwrap();
    assert!(matches!(local, DispersionTrendFit::Local(_)));
    assert_eq!(local.disp_fit().len(), means.len());
    assert!(local.use_for_mean().is_none());
    assert_eq!(local.genes_used_for_mean(), None);
    assert_eq!(local.used_min_disp_floor(), Some(false));

    let mean = fit_dispersion_trend(&means, &disps, FitType::Mean).unwrap();
    assert!(matches!(mean, DispersionTrendFit::Mean(_)));
    assert_eq!(mean.disp_fit().len(), means.len());
    assert_eq!(mean.use_for_fit(), vec![true; means.len()].as_slice());
    assert_eq!(
        mean.use_for_mean().unwrap(),
        vec![true; means.len()].as_slice()
    );
    assert_eq!(mean.genes_used_for_mean(), Some(means.len()));
    assert_eq!(mean.used_min_disp_floor(), None);
}

#[test]
fn fit_dispersion_trend_rejects_glm_gam_poi_until_supported() {
    let (means, disps) = trend_data();
    let err = fit_dispersion_trend(&means, &disps, FitType::GlmGamPoi).unwrap_err();

    assert!(matches!(err, DeseqError::UnsupportedFeature { .. }));
}

#[test]
fn mean_trend_evaluates_constant_for_positive_means() {
    let trend = MeanDispersionTrend { mean_disp: 0.2 };

    assert_relative_eq!(trend.evaluate(10.0).unwrap(), 0.2, epsilon = 1e-12);
    let fitted = trend
        .evaluate_many_allow_missing(&[10.0, 0.0, f64::NAN, 100.0])
        .unwrap();
    assert_relative_eq!(fitted[0], 0.2, epsilon = 1e-12);
    assert!(fitted[1].is_nan());
    assert!(fitted[2].is_nan());
    assert_relative_eq!(fitted[3], 0.2, epsilon = 1e-12);
}

#[test]
fn mean_trend_reports_trend_specific_invalid_fit_context() {
    let err = MeanDispersionTrend { mean_disp: 0.0 }
        .evaluate(10.0)
        .unwrap_err();

    assert!(err.to_string().contains("mean dispersion trend value"));
}

#[test]
fn local_trend_use_for_fit_matches_deseq2_threshold() {
    let use_for_fit = local_trend_use_for_fit(
        &[0.0, 10.0, 20.0, f64::NAN, 40.0],
        &[1.0, 5e-8, 1e-7, 1.0, f64::NAN],
        1e-8,
    )
    .unwrap();

    assert_eq!(use_for_fit, vec![false, false, true, false, false]);
}

#[test]
fn fit_local_dispersion_trend_recovers_log_linear_curve() {
    let means = vec![10.0, 20.0, 40.0, 80.0, 160.0, 320.0];
    let disps = means
        .iter()
        .map(|mean| 0.8 * f64::powf(*mean, -0.35))
        .collect::<Vec<_>>();

    let fit = fit_local_dispersion_trend(
        &means,
        &disps,
        LocalDispersionTrendOptions {
            degree: 1,
            ..LocalDispersionTrendOptions::default()
        },
    )
    .unwrap();

    assert_eq!(fit.genes_used, means.len());
    assert!(!fit.used_min_disp_floor);
    for (expected, fitted) in disps.iter().copied().zip(fit.disp_fit.iter().copied()) {
        assert_relative_eq!(fitted, expected, epsilon = 1e-12);
    }
}

#[test]
fn fit_local_dispersion_trend_sorts_fit_rows_once() {
    let means = vec![160.0, 10.0, 80.0, 20.0, 40.0];
    let disps = means
        .iter()
        .map(|mean| 0.8 * f64::powf(*mean, -0.35))
        .collect::<Vec<_>>();

    let fit = fit_local_dispersion_trend(
        &means,
        &disps,
        LocalDispersionTrendOptions {
            degree: 1,
            ..LocalDispersionTrendOptions::default()
        },
    )
    .unwrap();

    assert!(fit
        .trend
        .log_means
        .windows(2)
        .all(|window| window[0] <= window[1]));
    for (mean, expected) in means.iter().copied().zip(disps.iter().copied()) {
        assert_relative_eq!(fit.trend.evaluate(mean).unwrap(), expected, epsilon = 1e-12);
    }
}

#[test]
fn fit_local_dispersion_trend_handles_duplicated_means_with_lower_degree_fallback() {
    let means = vec![10.0, 10.0, 10.0, 40.0, 80.0];
    let disps = vec![0.08, 0.12, 0.16, 0.2, 0.28];

    let fit = fit_local_dispersion_trend(
        &means,
        &disps,
        LocalDispersionTrendOptions {
            span: 0.6,
            degree: 2,
            ..LocalDispersionTrendOptions::default()
        },
    )
    .unwrap();

    assert_eq!(fit.genes_used, means.len());
    assert!(!fit.used_min_disp_floor);
    assert_eq!(fit.trend.degree, 2);
    assert!(fit.disp_fit.iter().all(|value| value.is_finite()));
    let expected_duplicate_fit = ((0.08_f64.ln() + 0.12_f64.ln() + 0.16_f64.ln()) / 3.0).exp();
    assert_relative_eq!(
        fit.trend.evaluate(10.0).unwrap(),
        expected_duplicate_fit,
        epsilon = 1e-12
    );
}

#[test]
fn local_dispersion_trend_rejects_manually_unsorted_fit_rows() {
    let trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 0.7,
        degree: 1,
        log_means: vec![2.0, 1.0],
        log_disps: vec![-1.0, -2.0],
        weights: vec![10.0, 20.0],
    };

    assert!(trend.evaluate(10.0).is_err());
}

#[test]
fn local_dispersion_trend_uses_weighted_constant_fallback_for_rank_degenerate_window() {
    let log_mean = 10.0_f64.ln();
    let log_disps = vec![0.08_f64.ln(), 0.12_f64.ln(), 0.2_f64.ln()];
    let weights = vec![1.0, 2.0, 7.0];
    let trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 2,
        log_means: vec![log_mean, log_mean, log_mean],
        log_disps,
        weights,
    };

    let expected = ((1.0 * 0.08_f64.ln() + 2.0 * 0.12_f64.ln() + 7.0 * 0.2_f64.ln()) / 10.0).exp();

    assert_relative_eq!(trend.evaluate(10.0).unwrap(), expected, epsilon = 1e-12);
}

#[test]
fn local_dispersion_trend_keeps_large_weighted_constant_fit_finite() {
    let log_mean = 10.0_f64.ln();
    let trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 2,
        log_means: vec![log_mean, log_mean, log_mean],
        log_disps: vec![600.0, 600.0, 600.0],
        weights: vec![10.0, 20.0, 30.0],
    };

    assert_relative_eq!(
        trend.evaluate(10.0).unwrap(),
        600.0_f64.exp(),
        epsilon = 1e-12
    );
}

#[test]
fn local_dispersion_trend_uses_nearest_fallback_when_manual_weights_overflow() {
    let log_mean = 10.0_f64.ln();
    let trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 2,
        log_means: vec![log_mean, log_mean, log_mean],
        log_disps: vec![0.1_f64.ln(), 0.2_f64.ln(), 0.3_f64.ln()],
        weights: vec![f64::MAX, f64::MAX, f64::MAX],
    };

    assert_relative_eq!(trend.evaluate(10.0).unwrap(), 0.1, epsilon = 1e-12);
}

#[test]
fn local_dispersion_trend_rejects_overflowed_fitted_value() {
    let trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 0,
        log_means: vec![10.0_f64.ln(), 20.0_f64.ln()],
        log_disps: vec![1000.0, 1000.0],
        weights: vec![10.0, 20.0],
    };

    let err = trend.evaluate(10.0).unwrap_err();
    assert!(err.to_string().contains("local dispersion trend produced"));
    assert!(trend.evaluate_many_allow_missing(&[0.0, 10.0]).is_err());
}

#[test]
fn local_dispersion_trend_rejects_invalid_manual_fit_values_and_weights() {
    let mut trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 1,
        log_means: vec![10.0_f64.ln(), 20.0_f64.ln()],
        log_disps: vec![0.1_f64.ln(), 0.2_f64.ln()],
        weights: vec![10.0, 20.0],
    };

    trend.log_disps[1] = f64::NAN;
    assert!(trend.evaluate(10.0).is_err());

    trend.log_disps[1] = 0.2_f64.ln();
    trend.weights[0] = 0.0;
    assert!(trend.evaluate(10.0).is_err());

    trend.weights[0] = f64::INFINITY;
    assert!(trend.evaluate(10.0).is_err());
}

#[test]
fn local_dispersion_trend_rejects_invalid_manual_shape_options() {
    let mut trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 1,
        log_means: vec![10.0_f64.ln(), 20.0_f64.ln()],
        log_disps: vec![0.1_f64.ln(), 0.2_f64.ln()],
        weights: vec![10.0, 20.0],
    };

    trend.span = 0.0;
    assert!(trend.evaluate(10.0).is_err());

    trend.span = f64::NAN;
    assert!(trend.evaluate(10.0).is_err());

    trend.span = 1.0;
    trend.degree = 3;
    assert!(trend.evaluate(10.0).is_err());
}

#[test]
fn local_dispersion_trend_rejects_inconsistent_empty_manual_fit_state() {
    let mut trend = LocalDispersionTrend {
        min_disp: 1e-8,
        span: 1.0,
        degree: 1,
        log_means: Vec::new(),
        log_disps: Vec::new(),
        weights: Vec::new(),
    };

    assert_relative_eq!(trend.evaluate(10.0).unwrap(), 1e-8, epsilon = 1e-15);

    trend.log_disps.push(0.1_f64.ln());
    assert!(trend.evaluate(10.0).is_err());

    trend.log_disps.clear();
    trend.weights.push(1.0);
    assert!(trend.evaluate(10.0).is_err());
}

#[test]
fn fit_local_dispersion_trend_uses_min_disp_floor_when_all_estimates_are_near_minimum() {
    let fit = fit_local_dispersion_trend(
        &[10.0, 20.0, 0.0],
        &[1e-8, 2e-8, 3e-8],
        LocalDispersionTrendOptions::default(),
    )
    .unwrap();

    assert_eq!(fit.genes_used, 0);
    assert!(fit.used_min_disp_floor);
    assert_relative_eq!(fit.disp_fit[0], 1e-8, epsilon = 1e-15);
    assert_relative_eq!(fit.disp_fit[1], 1e-8, epsilon = 1e-15);
    assert!(fit.disp_fit[2].is_nan());
}

#[test]
fn local_dispersion_trend_enum_reports_min_disp_floor_usage() {
    let fit =
        fit_dispersion_trend(&[10.0, 20.0, 0.0], &[1e-8, 2e-8, 3e-8], FitType::Local).unwrap();

    assert!(matches!(fit, DispersionTrendFit::Local(_)));
    assert_eq!(fit.used_min_disp_floor(), Some(true));
    assert_eq!(fit.genes_used(), 0);
}

#[test]
fn fit_local_dispersion_trend_returns_nan_for_missing_fitted_rows() {
    let fit = fit_local_dispersion_trend(
        &[10.0, 20.0, f64::NAN, 40.0],
        &[0.1, 0.15, 0.2, 0.3],
        LocalDispersionTrendOptions::default(),
    )
    .unwrap();

    assert!(fit.disp_fit[0].is_finite());
    assert!(fit.disp_fit[1].is_finite());
    assert!(fit.disp_fit[2].is_nan());
    assert!(fit.disp_fit[3].is_finite());
}

#[test]
fn fit_local_dispersion_trend_uses_constant_fit_for_single_usable_row() {
    let fit = fit_local_dispersion_trend(
        &[10.0, 20.0, 40.0],
        &[1e-8, 0.2, 1e-8],
        LocalDispersionTrendOptions::default(),
    )
    .unwrap();

    assert_eq!(fit.genes_used, 1);
    assert_eq!(fit.use_for_fit, vec![false, true, false]);
    for value in fit.disp_fit {
        assert_relative_eq!(value, 0.2, epsilon = 1e-12);
    }
}

#[test]
fn fit_local_dispersion_trend_validates_inputs() {
    assert!(fit_local_dispersion_trend(
        &[10.0, 20.0],
        &[0.1],
        LocalDispersionTrendOptions::default(),
    )
    .is_err());
    for span in [0.0, 1.1, f64::NAN, f64::INFINITY] {
        assert!(fit_local_dispersion_trend(
            &[10.0, 20.0],
            &[0.1, 0.2],
            LocalDispersionTrendOptions {
                span,
                ..LocalDispersionTrendOptions::default()
            },
        )
        .is_err());
    }
    assert!(fit_local_dispersion_trend(
        &[10.0, 20.0],
        &[0.1, 0.2],
        LocalDispersionTrendOptions {
            degree: 3,
            ..LocalDispersionTrendOptions::default()
        },
    )
    .is_err());
}

#[test]
fn mean_trend_use_for_mean_matches_deseq2_threshold() {
    let use_for_mean = mean_trend_use_for_mean(
        &[0.0, 10.0, 20.0, f64::NAN, 40.0],
        &[1.0, 5e-8, 2e-7, 1.0, f64::NAN],
        1e-8,
    )
    .unwrap();

    assert_eq!(use_for_mean, vec![false, false, true, false, false]);
}

#[test]
fn mean_trend_use_for_fit_matches_deseq2_viability_threshold() {
    let use_for_fit = mean_trend_use_for_fit(
        &[0.0, 10.0, 20.0, 40.0, f64::NAN],
        &[1.0, 1e-6, 1.0000000001e-6, f64::NAN, 2e-6],
        1e-8,
    )
    .unwrap();

    assert_eq!(use_for_fit, vec![false, false, true, false, false]);

    let err = mean_trend_use_for_fit(&[10.0, 20.0], &[0.1], 1e-8).unwrap_err();
    assert!(err.to_string().contains("mean dispersion trend fit rows"));
}

#[test]
fn fit_mean_dispersion_trend_uses_deseq2_trimmed_mean() {
    let fit = fit_mean_dispersion_trend(
        &[10.0, 20.0, 30.0, 40.0],
        &[0.1, 0.2, 0.3, 1000.0],
        MeanDispersionTrendOptions {
            trim: 0.25,
            ..MeanDispersionTrendOptions::default()
        },
    )
    .unwrap();

    assert_eq!(fit.genes_used, 4);
    assert_eq!(fit.genes_used_for_mean, 4);
    assert_eq!(fit.use_for_fit, vec![true, true, true, true]);
    assert_relative_eq!(fit.trend.mean_disp, 0.25, epsilon = 1e-12);
    for value in fit.disp_fit {
        assert_relative_eq!(value, 0.25, epsilon = 1e-12);
    }
}

#[test]
fn fit_mean_dispersion_trend_rejects_invalid_trim_controls() {
    let (means, disps) = trend_data();

    for trim in [-0.1, 0.5001, f64::NAN, f64::INFINITY] {
        assert!(fit_mean_dispersion_trend(
            &means,
            &disps,
            MeanDispersionTrendOptions {
                trim,
                ..MeanDispersionTrendOptions::default()
            },
        )
        .is_err());
    }
}

#[test]
fn fit_mean_dispersion_trend_keeps_large_finite_averages_finite() {
    let large = f64::MAX / 2.0;
    let fit = fit_mean_dispersion_trend(
        &[10.0, 20.0, 30.0, 40.0],
        &[large, large, large, large],
        MeanDispersionTrendOptions::default(),
    )
    .unwrap();

    assert!(fit.trend.mean_disp.is_finite());
    assert_relative_eq!(fit.trend.mean_disp, large, max_relative = 1e-12);
    for value in fit.disp_fit {
        assert_relative_eq!(value, large, max_relative = 1e-12);
    }
}

#[test]
fn fit_mean_dispersion_trend_keeps_large_even_median_finite() {
    let large = f64::MAX / 2.0;
    let fit = fit_mean_dispersion_trend(
        &[10.0, 20.0, 30.0, 40.0],
        &[large, large, large, large],
        MeanDispersionTrendOptions {
            trim: 0.5,
            ..MeanDispersionTrendOptions::default()
        },
    )
    .unwrap();

    assert!(fit.trend.mean_disp.is_finite());
    assert_relative_eq!(fit.trend.mean_disp, large, max_relative = 1e-12);
}

#[test]
fn fit_mean_dispersion_trend_returns_nan_for_missing_fitted_rows() {
    let fit = fit_mean_dispersion_trend(
        &[10.0, 0.0, f64::NAN, 40.0],
        &[0.1, 0.2, 0.3, 0.4],
        MeanDispersionTrendOptions::default(),
    )
    .unwrap();

    assert_relative_eq!(fit.disp_fit[0], fit.trend.mean_disp, epsilon = 1e-12);
    assert!(fit.disp_fit[1].is_nan());
    assert!(fit.disp_fit[2].is_nan());
    assert_relative_eq!(fit.disp_fit[3], fit.trend.mean_disp, epsilon = 1e-12);
}

#[test]
fn fit_mean_dispersion_trend_keeps_deseq2_viability_gate() {
    let err = fit_mean_dispersion_trend(
        &[10.0, 20.0, 30.0],
        &[1e-8, 2e-8, 3e-8],
        MeanDispersionTrendOptions::default(),
    );

    assert!(err.is_err());
}

#[test]
fn builder_mean_fit_type_dispatch_populates_constant_disp_fit() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 30, 10, 30, //
            20, 60, 25, 55, //
            0, 0, 0, 0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let fit = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 4])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_dispersion_trend_linear_mu(&counts, &design)
        .unwrap();

    let disp_fit = fit.disp_fit.unwrap();
    assert_eq!(disp_fit.len(), counts.n_genes());
    assert!(disp_fit[0].is_finite());
    assert_relative_eq!(disp_fit[0], disp_fit[1], epsilon = 1e-12);
    assert!(disp_fit[2].is_nan());
}

#[test]
fn builder_local_fit_type_dispatch_populates_local_disp_fit() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![
            10, 30, 10, 30, //
            20, 60, 25, 55, //
            40, 95, 45, 105, //
            0, 0, 0, 0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let fit = DeseqBuilder::new()
        .fit_type(FitType::Local)
        .size_factors(vec![1.0; 4])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_dispersion_trend_linear_mu(&counts, &design)
        .unwrap();

    let disp_fit = fit.disp_fit.unwrap();
    assert_eq!(disp_fit.len(), counts.n_genes());
    assert!(disp_fit[0].is_finite());
    assert!(disp_fit[1].is_finite());
    assert!(disp_fit[2].is_finite());
    assert!(disp_fit[3].is_nan());
}
