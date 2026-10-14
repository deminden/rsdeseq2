use crate::errors::DeseqError;

/// Deterministic median over finite input values.
pub fn median(values: &[f64]) -> Result<f64, DeseqError> {
    if values.is_empty() {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "cannot compute median of an empty slice".to_string(),
        });
    }
    let mut sorted = Vec::with_capacity(values.len());
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "median input".to_string(),
                index: Some(idx),
                value,
            });
        }
        sorted.push(value);
    }
    sorted.sort_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        Ok(sorted[mid])
    } else {
        midpoint(sorted[mid - 1], sorted[mid]).ok_or_else(|| DeseqError::NonFiniteValue {
            context: "median midpoint".to_string(),
            index: None,
            value: sorted[mid - 1],
        })
    }
}

/// Median after dropping non-finite values. Returns `None` if no values remain.
pub fn median_finite(values: &[f64]) -> Option<f64> {
    let mut sorted = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    if sorted.is_empty() {
        return None;
    }
    sorted.sort_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        Some(sorted[mid])
    } else {
        midpoint(sorted[mid - 1], sorted[mid])
    }
}

fn midpoint(left: f64, right: f64) -> Option<f64> {
    let value = left / 2.0 + right / 2.0;
    value.is_finite().then_some(value)
}

#[cfg(test)]
mod tests {
    use super::{median, median_finite};

    #[test]
    fn median_odd_length() {
        assert_eq!(median(&[3.0, 1.0, 2.0]).unwrap(), 2.0);
    }

    #[test]
    fn median_even_length() {
        assert_eq!(median(&[4.0, 1.0, 2.0, 3.0]).unwrap(), 2.5);
    }

    #[test]
    fn median_even_large_values_avoids_midpoint_overflow() {
        assert_eq!(median(&[f64::MAX, f64::MAX]).unwrap(), f64::MAX);
        assert_eq!(median_finite(&[f64::MAX, f64::MAX]), Some(f64::MAX));
        assert_eq!(median(&[-f64::MAX, f64::MAX]).unwrap(), 0.0);
        assert_eq!(median_finite(&[-f64::MAX, f64::MAX]), Some(0.0));
    }

    #[test]
    fn median_finite_filters_nan() {
        assert_eq!(median_finite(&[f64::NAN, 1.0, 3.0]), Some(2.0));
    }
}
