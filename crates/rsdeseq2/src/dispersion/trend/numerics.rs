//! Checked arithmetic and shared trend input validation.

use crate::errors::DeseqError;

pub(super) fn checked_add_assign(
    total: &mut f64,
    value: f64,
    context: &str,
) -> Result<(), DeseqError> {
    let next = *total + value;
    if !value.is_finite() || !next.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    *total = next;
    Ok(())
}

pub(super) fn checked_product2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let product = left * right;
    if !left.is_finite() || !right.is_finite() || !product.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(product)
}

pub(super) fn checked_product3(
    left: f64,
    middle: f64,
    right: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    checked_product2(checked_product2(left, middle, context)?, right, context)
}

pub(super) fn checked_sum2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let sum = left + right;
    if !left.is_finite() || !right.is_finite() || !sum.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(sum)
}

pub(super) fn checked_sub2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let difference = left - right;
    if !left.is_finite() || !right.is_finite() || !difference.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(difference)
}

pub(super) fn checked_div2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let quotient = left / right;
    if !left.is_finite() || !right.is_finite() || right == 0.0 || !quotient.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(quotient)
}

pub(super) fn checked_product_difference(
    left_a: f64,
    left_b: f64,
    right_a: f64,
    right_b: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    let left = left_a * left_b;
    let right = right_a * right_b;
    let difference = left - right;
    if !left.is_finite() || !right.is_finite() || !difference.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(difference)
}

pub(super) fn stable_midpoint(left: f64, right: f64) -> f64 {
    left / 2.0 + right / 2.0
}

pub(super) fn validate_positive_mean(mean: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion trend mean".to_string(),
            index,
            value: mean,
        });
    }
    Ok(())
}

pub(super) fn validate_positive_dispersion(
    dispersion: f64,
    index: Option<usize>,
) -> Result<(), DeseqError> {
    if !dispersion.is_finite() || dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: format!(
                "gene-wise dispersion at index {} must be finite and positive",
                index
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            ),
        });
    }
    Ok(())
}
