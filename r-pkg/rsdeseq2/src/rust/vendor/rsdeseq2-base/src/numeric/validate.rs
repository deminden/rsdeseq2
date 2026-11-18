//! Scalar input checks shared by the GLM stages.

use crate::errors::DeseqError;

/// Reject a value that is not finite and strictly positive.
///
/// The error is `NonFiniteValue` with the given context, index and the value.
#[inline]
pub fn positive_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

/// Reject a value that is not finite and non-negative.
///
/// The error is `NonFiniteValue` with the given context, index and the value.
#[inline]
pub fn nonnegative_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundaries_and_errors() {
        assert!(positive_finite(f64::MIN_POSITIVE, "x", 0).is_ok());
        assert!(positive_finite(0.0, "x", 0).is_err());
        assert!(nonnegative_finite(0.0, "x", 0).is_ok());
        assert!(nonnegative_finite(-0.0, "x", 0).is_ok());
        assert!(nonnegative_finite(-f64::MIN_POSITIVE, "x", 0).is_err());
        for value in [f64::NAN, f64::INFINITY] {
            assert!(positive_finite(value, "x", 0).is_err());
            assert!(nonnegative_finite(value, "x", 0).is_err());
        }
        let error = positive_finite(-1.0, "size factor", 3).unwrap_err();
        assert!(matches!(
            error,
            DeseqError::NonFiniteValue { ref context, index: Some(3), value }
                if context == "size factor" && value == -1.0
        ));
    }
}
