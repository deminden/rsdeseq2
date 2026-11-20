//! Checked arithmetic and shared trend input validation.

use crate::errors::{CalculationLocation, DeseqError, ErrorKind, invalid_calculation};
use crate::numeric::checked::{self, CheckedExt};

pub(super) fn checked_add_assign(
    total: &mut f64,
    value: f64,
    context: &str,
    location: CalculationLocation,
) -> Result<(), DeseqError> {
    *total = checked::try_add(*total, value).invalid_dispersion_at(context, location)?;
    Ok(())
}

pub(super) fn checked_product3(left: f64, middle: f64, right: f64) -> checked::Checked {
    checked::try_mul(checked::try_mul(left, middle)?, right)
}

pub(super) fn checked_product_difference(
    left_a: f64,
    left_b: f64,
    right_a: f64,
    right_b: f64,
) -> checked::Checked {
    let left = checked::try_mul(left_a, left_b)?;
    let right = checked::try_mul(right_a, right_b)?;
    checked::try_sub(left, right)
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

pub(super) fn validate_fitted_dispersion(
    trend: &'static str,
    gene: Option<usize>,
    value: f64,
) -> Result<f64, DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(invalid_calculation(
            ErrorKind::Dispersion,
            trend,
            gene.map(CalculationLocation::Gene),
            value,
            format_args!("{trend} produced a non-finite or non-positive fitted value"),
        ));
    }
    Ok(value)
}

pub(super) fn validate_positive_dispersion(
    dispersion: f64,
    index: Option<usize>,
) -> Result<(), DeseqError> {
    if !dispersion.is_finite() || dispersion <= 0.0 {
        return Err(invalid_calculation(
            ErrorKind::Dispersion,
            "gene-wise dispersion",
            index.map(CalculationLocation::Gene),
            dispersion,
            format_args!(
                "gene-wise dispersion at index {} must be finite and positive",
                index
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            ),
        ));
    }
    Ok(())
}
