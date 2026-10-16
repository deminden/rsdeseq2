// Private numerical regression tests.

#[cfg(test)]
mod tests {
    use super::{
        checked_div2, checked_mean_dispersion_working_weight, checked_scaled_sum,
        intercept_working_weights, stable_interpolate, weighted_mean,
    };
    use crate::errors::DeseqError;

    #[test]
    fn stable_interpolate_preserves_equal_extreme_endpoints() {
        let interpolated = stable_interpolate(f64::MAX, f64::MAX, 0.5);

        assert_eq!(interpolated, f64::MAX);
    }

    #[test]
    fn stable_interpolate_rejects_overflowed_delta() {
        let interpolated = stable_interpolate(-f64::MAX, f64::MAX, 0.5);

        assert!(interpolated.is_nan());
    }

    #[test]
    fn quantile_type7_rejects_overflowed_interpolation() {
        let err = super::quantile_type7(vec![-f64::MAX, f64::MAX], 0.5).unwrap_err();

        assert!(matches!(
            err,
            crate::errors::DeseqError::InvalidOptions { reason }
                if reason == "beta prior variance quantile produced non-finite quantile"
        ));
    }

    #[test]
    fn checked_scaled_sum_rejects_overflowed_rescale() {
        assert_eq!(checked_scaled_sum([f64::MAX, f64::MAX]), None);
    }

    #[test]
    fn checked_div2_rejects_zero_and_nonfinite_inputs() {
        assert_eq!(checked_div2(1.0, 0.0), None);
        assert_eq!(checked_div2(f64::NAN, 1.0), None);
        assert_eq!(checked_div2(4.0, 2.0), Some(2.0));
    }

    #[test]
    fn mean_dispersion_weight_rejects_nonfinite_arithmetic() {
        let err = checked_mean_dispersion_working_weight(
            f64::MIN_POSITIVE,
            f64::MAX,
            2,
            "test mean-dispersion weight",
        )
        .unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test mean-dispersion weight" && index == Some(2)
        ));
    }

    #[test]
    fn intercept_working_weights_reject_nonfinite_weight_scaling() {
        let err =
            intercept_working_weights(&[f64::MAX], f64::MIN_POSITIVE, Some(&[10.0])).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "intercept working weight" && index == Some(0)
        ));
    }

    #[test]
    fn weighted_mean_rejects_nonfinite_products() {
        let err = weighted_mean(&[f64::MAX], &[2.0], 0).unwrap_err();

        assert!(
            err.to_string()
                .contains("gene 0 has non-finite weighted normalized mean")
        );
    }
}
