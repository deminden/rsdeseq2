use crate::errors::Error;

/// Deterministic median over finite input values.
pub fn median(values: &[f64]) -> Result<f64, Error> {
    median_in_place(&mut values.to_vec())
}

/// Deterministic median over finite input values, partitioning the input in place.
pub fn median_in_place(values: &mut [f64]) -> Result<f64, Error> {
    if values.is_empty() {
        return Err(Error::InvalidDimensions {
            context: "median input".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(Error::NonFiniteValue {
                context: "median input".to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    let (lower, upper) = middle_values(values);
    lower
        .map_or(Some(upper), |left| midpoint(left, upper))
        .ok_or_else(|| Error::NonFiniteValue {
            context: "median midpoint".to_string(),
            index: None,
            value: lower.unwrap_or(upper),
        })
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
    let (lower, upper) = middle_values(&mut sorted);
    lower.map_or(Some(upper), |left| midpoint(left, upper))
}

// Select the upper midpoint, then find the lower one in the left partition.
// total_cmp preserves ordering of signed zero just like the former sort.
fn middle_values(values: &mut [f64]) -> (Option<f64>, f64) {
    let mid = values.len() / 2;
    let even = values.len().is_multiple_of(2);
    let (left, upper, _) = values.select_nth_unstable_by(mid, f64::total_cmp);
    let lower = if even {
        left.iter().copied().max_by(f64::total_cmp)
    } else {
        None
    };
    (lower, *upper)
}

fn midpoint(left: f64, right: f64) -> Option<f64> {
    let value = left / 2.0 + right / 2.0;
    value.is_finite().then_some(value)
}

#[cfg(test)]
mod tests {
    use super::{median, median_finite};
    use crate::errors::Error;

    #[test]
    fn median_rejects_empty_input_as_dimensions() {
        assert!(matches!(
            median(&[]),
            Err(Error::InvalidDimensions { ref context, expected: 1, actual: 0 })
                if context == "median input"
        ));
    }

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
