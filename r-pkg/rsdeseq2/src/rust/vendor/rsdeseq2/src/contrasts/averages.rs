//! Mean-difference weights for expanded-model shrinkage list contrasts.
use crate::errors::DeseqError;
use std::collections::HashSet;

/// Build a difference between means of distinct selected coefficient groups.
/// Indices are zero-based; either side may be empty, but their union must be
/// nonempty. Selection/name alignment is separate from this shared arithmetic.
pub fn averaged_list_contrast(
    coefficients: usize,
    positive: &[usize],
    negative: &[usize],
) -> Result<Vec<f64>, DeseqError> {
    let mut selected = HashSet::new();
    if coefficients == 0
        || positive.is_empty() && negative.is_empty()
        || positive
            .iter()
            .chain(negative)
            .any(|&i| i >= coefficients || !selected.insert(i))
    {
        return Err(DeseqError::InvalidOptions {
            reason: "average contrast requires distinct in-range indices and a nonempty selection"
                .into(),
        });
    }
    let mut values = vec![0.0; coefficients];
    for (indices, numerator) in [(positive, 1.0), (negative, -1.0)] {
        if !indices.is_empty() {
            let weight = numerator / indices.len() as f64;
            for &index in indices {
                values[index] = weight;
            }
        }
    }
    Ok(values)
}
