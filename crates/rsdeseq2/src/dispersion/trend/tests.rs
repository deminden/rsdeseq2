//! Numerical edge cases for trend fitting.

use super::numerics::{
    checked_add_assign, checked_div2, checked_product_difference, stable_midpoint,
};
use super::parametric::{
    ParametricDispersionTrend, ParametricDispersionTrendOptions, fitted_values_are_positive,
    gamma_deviance, interpolate_parametric_coefficient, parametric_coefficient_change,
    robust_parametric_trend_rows, weighted_gamma_identity_least_squares,
};
use approx::assert_relative_eq;

#[test]
fn gamma_deviance_keeps_near_perfect_fit_precision() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: 2.0,
    };
    let means = [20.0, 40.0, 80.0, 160.0];
    let relative_delta = 1.0e-8;
    let disps = means
        .iter()
        .map(|mean| trend.evaluate(*mean).unwrap() * (1.0 + relative_delta))
        .collect::<Vec<_>>();
    let good = [true; 4];

    let deviance = gamma_deviance(&means, &disps, &good, trend).unwrap();
    let expected = means.len() as f64 * 2.0 * (relative_delta - relative_delta.ln_1p());

    assert!(deviance.is_finite());
    assert!(deviance > 0.0);
    assert_relative_eq!(deviance, expected, max_relative = 1e-12);
}

#[test]
fn gamma_deviance_rejects_nonfinite_accumulation() {
    let trend = ParametricDispersionTrend {
        asympt_disp: f64::MIN_POSITIVE,
        extra_pois: f64::MIN_POSITIVE,
    };
    let means = [1.0, 2.0];
    let disps = [f64::MAX, f64::MAX];
    let good = [true, true];

    let err = gamma_deviance(&means, &disps, &good, trend).unwrap_err();

    assert!(
        err.to_string()
            .contains("parametric dispersion gamma deviance")
    );
}

#[test]
fn checked_div2_rejects_zero_and_nonfinite_quotients() {
    let zero_err = checked_div2(1.0, 0.0, "test division").unwrap_err();
    let overflow_err = checked_div2(f64::MAX, f64::MIN_POSITIVE, "test division").unwrap_err();

    assert!(zero_err.to_string().contains("test division"));
    assert!(overflow_err.to_string().contains("test division"));
}

#[test]
fn stable_midpoint_handles_opposite_extreme_endpoints() {
    assert_eq!(stable_midpoint(-f64::MAX, f64::MAX), 0.0);
    assert_eq!(stable_midpoint(f64::MAX, f64::MAX), f64::MAX);
}

#[test]
fn parametric_work_helpers_validate_aligned_rows() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: 2.0,
    };
    let means = [20.0, 40.0];
    let disps = [0.15];
    let good = [true, true];

    assert!(
        robust_parametric_trend_rows(
            &means,
            &disps,
            trend,
            ParametricDispersionTrendOptions::default(),
        )
        .is_err()
    );
    assert!(gamma_deviance(&means, &disps, &good, trend).is_err());
    assert!(weighted_gamma_identity_least_squares(&means, &disps, &good, trend).is_err());
    assert!(fitted_values_are_positive(&means, &[true], trend).is_err());
}

#[test]
fn weighted_gamma_identity_least_squares_rejects_nonfinite_coefficients() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.05,
        extra_pois: 2.0,
    };
    let means = [20.0, 40.0];
    let disps = [f64::MAX, f64::MAX];
    let good = [true, true];

    assert!(weighted_gamma_identity_least_squares(&means, &disps, &good, trend).is_err());
}

#[test]
fn weighted_gamma_identity_least_squares_rejects_nonfinite_weights() {
    let trend = ParametricDispersionTrend {
        asympt_disp: f64::MIN_POSITIVE,
        extra_pois: f64::MIN_POSITIVE,
    };
    let means = [1.0, 2.0];
    let disps = [0.1, 0.2];
    let good = [true, true];

    let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

    assert!(
        err.to_string()
            .contains("parametric dispersion weighted least-squares")
    );
}

#[test]
fn weighted_gamma_identity_least_squares_rejects_nonfinite_determinant_product() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 1.0,
        extra_pois: f64::MIN_POSITIVE,
    };
    let means = [1.0e-154, 2.0e-154];
    let disps = [0.1, 0.2];
    let good = [true, true];

    let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

    assert!(
        err.to_string()
            .contains("parametric dispersion weighted least-squares determinant")
    );
}

#[test]
fn weighted_gamma_identity_least_squares_rejects_nonfinite_numerator_product() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 1.0,
        extra_pois: f64::MIN_POSITIVE,
    };
    let means = [1.0, 0.5];
    let disps = [f64::MAX / 4.0, f64::MAX / 4.0];
    let good = [true, true];

    let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

    assert!(
        err.to_string()
            .contains("parametric dispersion weighted least-squares coefficient numerator")
    );
}

#[test]
fn checked_add_assign_rejects_nonfinite_sum() {
    let mut total = f64::MAX;

    let err = checked_add_assign(&mut total, f64::MAX, "test accumulation").unwrap_err();

    assert!(err.to_string().contains("test accumulation"));
    assert_eq!(total, f64::MAX);
}

#[test]
fn checked_product_difference_rejects_nonfinite_product() {
    let err = checked_product_difference(f64::MAX, 2.0, 1.0, 1.0, "test product").unwrap_err();

    assert!(err.to_string().contains("test product"));
}

#[test]
fn checked_product_difference_matches_finite_difference() {
    let observed = checked_product_difference(6.0, 7.0, 5.0, 3.0, "test product").unwrap();

    assert_eq!(observed, 27.0);
}

#[test]
fn interpolate_parametric_coefficient_avoids_difference_overflow() {
    let interpolated = interpolate_parametric_coefficient(f64::MAX, -f64::MAX, 0.5);

    assert_eq!(interpolated, 0.0);
}

#[test]
fn parametric_coefficient_change_rejects_nonfinite_ratio() {
    let old = ParametricDispersionTrend {
        asympt_disp: f64::MIN_POSITIVE,
        extra_pois: 1.0,
    };
    let new = ParametricDispersionTrend {
        asympt_disp: f64::MAX,
        extra_pois: 1.0,
    };

    let err = parametric_coefficient_change(old, new).unwrap_err();

    assert!(
        err.to_string()
            .contains("parametric dispersion coefficient change")
    );
}

#[test]
fn parametric_coefficient_change_matches_log_squared_update() {
    let old = ParametricDispersionTrend {
        asympt_disp: 0.1,
        extra_pois: 2.0,
    };
    let new = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.0,
    };

    let observed = parametric_coefficient_change(old, new).unwrap();
    let expected = 2.0 * 2.0_f64.ln().powi(2);

    assert_relative_eq!(observed, expected, epsilon = 1e-15);
}
