use crate::errors::DeseqError;

/// Benjamini-Hochberg adjusted p-values with missing-value support.
///
/// `None`, non-finite values, and values outside `[0, 1]` are omitted from
/// ranking and returned as `None`.
pub fn bh_adjust(pvalues: &[Option<f64>]) -> Vec<Option<f64>> {
    let mut indexed = pvalues
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(idx, value)| {
            value.and_then(|pvalue| {
                if pvalue.is_finite() && (0.0..=1.0).contains(&pvalue) {
                    Some((idx, pvalue))
                } else {
                    None
                }
            })
        })
        .collect::<Vec<_>>();

    indexed.sort_by(|left, right| left.1.total_cmp(&right.1));
    let m = indexed.len();
    let mut adjusted = vec![None; pvalues.len()];
    let mut running_min = 1.0_f64;
    for (rank_from_zero, (idx, pvalue)) in indexed.into_iter().enumerate().rev() {
        let rank = rank_from_zero + 1;
        let raw_adjusted = pvalue * m as f64 / rank as f64;
        running_min = running_min.min(raw_adjusted);
        adjusted[idx] = Some(running_min.min(1.0));
    }
    adjusted
}

/// Benjamini-Hochberg adjusted p-values for fully observed `f64` p-values.
pub fn bh_adjust_f64(pvalues: &[f64]) -> Result<Vec<f64>, DeseqError> {
    for (idx, value) in pvalues.iter().copied().enumerate() {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("p-value at index {idx} must be finite and within [0, 1]"),
            });
        }
    }
    let adjusted = bh_adjust(&pvalues.iter().copied().map(Some).collect::<Vec<_>>());
    let mut out = Vec::with_capacity(adjusted.len());
    for value in adjusted {
        match value {
            Some(value) => out.push(value),
            None => {
                return Err(DeseqError::InvalidSizeFactors {
                    reason: "validated p-value unexpectedly produced a missing adjustment"
                        .to_string(),
                });
            }
        }
    }
    Ok(out)
}
