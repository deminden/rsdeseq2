use crate::errors::DeseqError;

/// Trigamma function for positive finite inputs.
///
/// DESeq2 uses R's `trigamma((m - p) / 2)` when estimating dispersion prior
/// variance. This implementation uses recurrence to shift to a stable
/// asymptotic region and then evaluates a Bernoulli expansion.
pub fn trigamma(mut x: f64) -> Result<f64, DeseqError> {
    if !x.is_finite() || x <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "trigamma input".to_string(),
            index: None,
            value: x,
        });
    }

    let mut result = 0.0;
    while x < 8.0 {
        result += x.recip().powi(2);
        x += 1.0;
    }

    let inv = x.recip();
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    let inv5 = inv3 * inv2;
    let inv7 = inv5 * inv2;
    let inv9 = inv7 * inv2;
    let inv11 = inv9 * inv2;
    let inv13 = inv11 * inv2;

    result += inv + 0.5 * inv2 + inv3 / 6.0 - inv5 / 30.0 + inv7 / 42.0 - inv9 / 30.0
        + 5.0 * inv11 / 66.0
        - 691.0 * inv13 / 2730.0;
    Ok(result)
}
