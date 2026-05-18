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
