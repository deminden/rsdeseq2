//! Checked floating-point arithmetic shared by every statistical stage.
//!
//! Each operation evaluates exactly one IEEE-754 expression in a fixed order and
//! accepts the result only when it is finite. Two forms exist:
//!
//! - `add`, `sub`, `mul`, `div`, `sum`, `scaled_sum`, `mean` return
//!   `Option<f64>` for callers that only need success or failure.
//! - `try_*` return `Result<f64, NonFinite>`, keeping the offending value and,
//!   for reductions, the element index, so error messages can report them.
//!
//! [`CheckedExt`] turns either form into a [`DeseqError`] with a context string
//! and an optional index.
//!
//! # Equivalence of input and result checks
//!
//! For `+`, `-` and `*`, a non-finite operand always yields a non-finite result
//! (NaN propagates, `inf - inf` and `inf * 0` are NaN, and any other infinite
//! operand gives an infinite result). Checking only the result is therefore the
//! same as also checking both operands, so one definition serves both of the
//! historical styles.
//!
//! Division differs: `x / inf` is a finite zero. [`try_div`] rejects non-finite
//! operands and a zero divisor; [`try_div_finite_result`] checks the quotient
//! only and is kept for the callers that historically accepted `x / inf`.
//! A zero divisor with finite operands already yields an infinite or NaN
//! quotient, so the explicit zero test in [`try_div`] never changes a result.

use crate::errors::DeseqError;

/// A rejected value, with the reduction index that produced it when known.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonFinite {
    pub value: f64,
    pub index: Option<usize>,
}

/// Result of a checked operation that keeps the rejected value.
pub type Checked = Result<f64, NonFinite>;

/// Rejection of `value` without a reduction index.
#[inline(always)]
pub fn reject(value: f64) -> Checked {
    Err(NonFinite { value, index: None })
}

#[inline(always)]
fn finite(value: f64) -> Checked {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NonFinite { value, index: None })
    }
}

/// `left + right`, finite.
#[inline(always)]
pub fn try_add(left: f64, right: f64) -> Checked {
    finite(left + right)
}

/// `left - right`, finite.
#[inline(always)]
pub fn try_sub(left: f64, right: f64) -> Checked {
    finite(left - right)
}

/// `left * right`, finite.
#[inline(always)]
pub fn try_mul(left: f64, right: f64) -> Checked {
    finite(left * right)
}

/// `left / right` with finite operands, a nonzero divisor and a finite quotient.
#[inline(always)]
pub fn try_div(left: f64, right: f64) -> Checked {
    let quotient = left / right;
    if left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite() {
        Ok(quotient)
    } else {
        Err(NonFinite {
            value: quotient,
            index: None,
        })
    }
}

/// `left / right` when the quotient is finite; accepts `x / inf == 0`.
#[inline(always)]
pub fn try_div_finite_result(left: f64, right: f64) -> Checked {
    finite(left / right)
}

/// Sequential sum from `0.0` in input order.
///
/// On failure the error holds the first non-finite partial sum and the index of
/// the element that produced it.
#[inline]
pub fn try_sum(values: impl IntoIterator<Item = f64>) -> Checked {
    let mut sum = 0.0;
    for (index, value) in values.into_iter().enumerate() {
        sum = try_add(sum, value).map_err(|error| NonFinite {
            index: Some(index),
            ..error
        })?;
    }
    Ok(sum)
}

/// Largest absolute value, rejecting the first non-finite element.
#[inline]
fn max_abs(values: &[f64]) -> Checked {
    let mut scale = 0.0_f64;
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(NonFinite {
                value,
                index: Some(index),
            });
        }
        scale = scale.max(value.abs());
    }
    Ok(scale)
}

/// Scale-normalized sum: `(Σ values[i] / scale) * scale`, `scale = max |values|`.
///
/// Dividing by the largest magnitude first keeps the running sum bounded by the
/// element count, so only the final rescaling can overflow. All-zero and empty
/// inputs return `0.0`.
pub fn try_scaled_sum(values: &[f64]) -> Checked {
    let scale = max_abs(values)?;
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized = try_sum(values.iter().map(|value| value / scale))?;
    finite(normalized * scale)
}

/// Scale-normalized mean: `((Σ values[i] / scale) / n) * scale`.
///
/// All-zero and empty inputs return `0.0`.
pub fn try_mean(values: &[f64]) -> Checked {
    let scale = max_abs(values)?;
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized = try_sum(values.iter().map(|value| value / scale))?;
    finite(normalized / values.len() as f64 * scale)
}

/// `left + right`, finite.
#[inline(always)]
pub fn add(left: f64, right: f64) -> Option<f64> {
    try_add(left, right).ok()
}

/// `left - right`, finite.
#[inline(always)]
pub fn sub(left: f64, right: f64) -> Option<f64> {
    try_sub(left, right).ok()
}

/// `left * right`, finite.
#[inline(always)]
pub fn mul(left: f64, right: f64) -> Option<f64> {
    try_mul(left, right).ok()
}

/// `left / right` with finite operands, a nonzero divisor and a finite quotient.
#[inline(always)]
pub fn div(left: f64, right: f64) -> Option<f64> {
    try_div(left, right).ok()
}

/// Sequential sum from `0.0` in input order.
#[inline]
pub fn sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    try_sum(values).ok()
}

/// See [`try_scaled_sum`].
pub fn scaled_sum(values: &[f64]) -> Option<f64> {
    try_scaled_sum(values).ok()
}

/// See [`try_mean`].
pub fn mean(values: &[f64]) -> Option<f64> {
    try_mean(values).ok()
}

/// Conversion of checked results into [`DeseqError`].
pub trait CheckedExt {
    /// `NonFiniteValue` with the index recorded by the operation, if any.
    fn nonfinite(self, context: &str) -> Result<f64, DeseqError>;

    /// `NonFiniteValue` with an explicit index, replacing any recorded one.
    fn nonfinite_at(self, context: &str, index: Option<usize>) -> Result<f64, DeseqError>;

    /// `InvalidDispersion` with the reason "`{context}` must remain finite".
    fn invalid_dispersion(self, context: &str) -> Result<f64, DeseqError>;
}

#[cold]
#[inline(never)]
fn nonfinite_error(context: &str, index: Option<usize>, value: f64) -> DeseqError {
    DeseqError::NonFiniteValue {
        context: context.to_string(),
        index,
        value,
    }
}

#[cold]
#[inline(never)]
fn dispersion_error(context: &str) -> DeseqError {
    DeseqError::InvalidDispersion {
        reason: format!("{context} must remain finite"),
    }
}

impl CheckedExt for Checked {
    #[inline(always)]
    fn nonfinite(self, context: &str) -> Result<f64, DeseqError> {
        self.map_err(|error| nonfinite_error(context, error.index, error.value))
    }

    #[inline(always)]
    fn nonfinite_at(self, context: &str, index: Option<usize>) -> Result<f64, DeseqError> {
        self.map_err(|error| nonfinite_error(context, index, error.value))
    }

    #[inline(always)]
    fn invalid_dispersion(self, context: &str) -> Result<f64, DeseqError> {
        self.map_err(|_| dispersion_error(context))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPECIALS: [f64; 9] = [
        0.0,
        -0.0,
        1.0,
        -2.5,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];

    /// The historical strict form, which also checked both operands.
    fn strict(op: fn(f64, f64) -> f64, left: f64, right: f64) -> Option<f64> {
        let value = op(left, right);
        (left.is_finite() && right.is_finite() && value.is_finite()).then_some(value)
    }

    fn same(left: Option<f64>, right: Option<f64>) -> bool {
        left.map(f64::to_bits) == right.map(f64::to_bits)
    }

    #[test]
    fn result_only_checks_equal_operand_checks_for_add_sub_mul() {
        for left in SPECIALS {
            for right in SPECIALS {
                assert!(same(add(left, right), strict(|a, b| a + b, left, right)));
                assert!(same(sub(left, right), strict(|a, b| a - b, left, right)));
                assert!(same(mul(left, right), strict(|a, b| a * b, left, right)));
                // The explicit zero-divisor test never changes the outcome.
                assert!(same(div(left, right), strict(|a, b| a / b, left, right)));
            }
        }
    }

    #[test]
    fn division_variants_differ_only_for_infinite_divisor() {
        assert_eq!(div(1.0, f64::INFINITY), None);
        assert_eq!(try_div_finite_result(1.0, f64::INFINITY), Ok(0.0));
        assert_eq!(div(1.0, 0.0), None);
        assert!(try_div_finite_result(1.0, 0.0).is_err());
    }

    #[test]
    fn sum_reports_partial_value_and_index() {
        assert_eq!(try_sum([1.0, 2.0, 3.0]), Ok(6.0));
        assert_eq!(try_sum([]), Ok(0.0));
        let error = try_sum([1.0, f64::MAX, f64::MAX]).unwrap_err();
        assert_eq!(error.index, Some(2));
        assert_eq!(error.value, f64::INFINITY);
    }

    #[test]
    fn scaled_reductions_match_reference_order() {
        let values = [3.0, -1.5, 0.25, 7.0];
        let scale: f64 = 7.0;
        let normalized = (((0.0 + 3.0 / scale) + -1.5 / scale) + 0.25 / scale) + 7.0 / scale;
        assert_eq!(
            scaled_sum(&values).map(f64::to_bits),
            Some((normalized * scale).to_bits())
        );
        assert_eq!(
            mean(&values).map(f64::to_bits),
            Some((normalized / 4.0 * scale).to_bits())
        );
        assert_eq!(mean(&[]), Some(0.0));
        assert_eq!(scaled_sum(&[0.0, -0.0]), Some(0.0));
        assert_eq!(mean(&[f64::MAX, f64::MAX]), Some(f64::MAX));
        let overflow = try_scaled_sum(&[f64::MAX, f64::MAX]).unwrap_err();
        assert_eq!((overflow.value, overflow.index), (f64::INFINITY, None));
        let input = try_mean(&[1.0, f64::NAN]).unwrap_err();
        assert!(input.value.is_nan());
        assert_eq!(input.index, Some(1));
    }

    #[test]
    fn conversions_keep_context_index_and_value() {
        let error = try_mul(f64::MAX, 2.0)
            .nonfinite_at("product", Some(4))
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            DeseqError::NonFiniteValue {
                context: "product".to_string(),
                index: Some(4),
                value: f64::INFINITY,
            }
            .to_string()
        );
        let error = try_sum([f64::INFINITY]).nonfinite("total").unwrap_err();
        assert!(matches!(
            error,
            DeseqError::NonFiniteValue { index: Some(0), value, .. } if value == f64::INFINITY
        ));
        let error = try_div(1.0, 0.0).invalid_dispersion("trend").unwrap_err();
        assert!(matches!(
            error,
            DeseqError::InvalidDispersion { ref reason } if reason == "trend must remain finite"
        ));
        let error = try_mean(&[1.0, f64::NEG_INFINITY])
            .nonfinite_at("mean", None)
            .unwrap_err();
        assert!(matches!(
            error,
            DeseqError::NonFiniteValue { index: None, value, .. } if value == f64::NEG_INFINITY
        ));
    }
}
