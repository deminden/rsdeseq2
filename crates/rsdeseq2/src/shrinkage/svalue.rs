use crate::DeseqError;

/// Average ordered local error rates and restore input order (s-values).
///
/// Ties keep input order, including signed zero. Missing entries remain missing
/// and do not enter the averages. Finite rates must lie within `[0, 1]`.
/// This is the ordering rule used by apeglm and the wrapper's thresholded ashr
/// route; it is distinct from a Benjamini-Hochberg adjustment.
pub fn svalues(rates: &[Option<f64>]) -> Result<Vec<Option<f64>>, DeseqError> {
    let mut order = Vec::with_capacity(rates.len());
    for (index, rate) in rates.iter().copied().enumerate() {
        if let Some(rate) = rate {
            if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "local error rate at index {index} must be finite and within [0, 1]"
                    ),
                });
            }
            order.push((index, rate));
        }
    }
    // partial_cmp deliberately treats -0 and +0 as equal; sorting is stable.
    order.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mut out = vec![None; rates.len()];
    let mut cumulative = 0.0;
    for (rank, (index, rate)) in order.into_iter().enumerate() {
        cumulative += rate;
        out[index] = Some(cumulative / (rank + 1) as f64);
    }
    Ok(out)
}
