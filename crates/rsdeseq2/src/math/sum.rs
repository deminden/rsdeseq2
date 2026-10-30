//! Shared checked scaled sum with fixed input-order reduction.

pub(crate) fn checked_scaled_sum(values: &[f64]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut normalized_sum = 0.0;
    for value in values.iter().copied() {
        let term = value / scale;
        let next = normalized_sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        normalized_sum = next;
    }
    let sum = normalized_sum * scale;
    sum.is_finite().then_some(sum)
}
