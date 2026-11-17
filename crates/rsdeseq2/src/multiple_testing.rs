use crate::errors::DeseqError;

/// Benjamini-Hochberg adjusted p-values with missing-value support.
///
/// `None`, non-finite values, and values outside `[0, 1]` are omitted from
/// ranking and returned as `None`.
pub fn bh_adjust(pvalues: &[Option<f64>]) -> Vec<Option<f64>> {
    PreparedBh::new(pvalues).adjusted()
}

/// A stable global p-value order can be reused for any filtered population.
/// Filtering preserves the relative order and ranks of each retained subset.
pub(crate) struct PreparedBh {
    indexed: Vec<(usize, f64)>,
    length: usize,
}

impl PreparedBh {
    pub(crate) fn new(pvalues: &[Option<f64>]) -> Self {
        let mut indexed = pvalues
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(index, value)| {
                value
                    .filter(|p| p.is_finite() && (0.0..=1.0).contains(p))
                    .map(|p| (index, p))
            })
            .collect::<Vec<_>>();
        indexed.sort_by(|left, right| left.1.total_cmp(&right.1));
        Self {
            indexed,
            length: pvalues.len(),
        }
    }

    fn visit_population(&self, keep: impl Fn(usize) -> bool, mut visit: impl FnMut(usize, f64)) {
        let m = self
            .indexed
            .iter()
            .filter(|(index, _)| keep(*index))
            .count();
        let mut rank = m;
        let mut running_min = 1.0_f64;
        for &(index, pvalue) in self.indexed.iter().rev().filter(|(index, _)| keep(*index)) {
            // One canonical order for every frontend, matching DESeq2's p.adjust.
            let raw_adjusted = (m as f64 / rank as f64) * pvalue;
            running_min = running_min.min(raw_adjusted);
            visit(index, running_min.min(1.0));
            rank -= 1;
        }
    }

    fn adjusted(&self) -> Vec<Option<f64>> {
        let mut adjusted = vec![None; self.length];
        self.visit_population(|_| true, |index, pvalue| adjusted[index] = Some(pvalue));
        adjusted
    }

    pub(crate) fn adjusted_for_filter(&self, filter: &[f64], cutoff: f64) -> Vec<Option<f64>> {
        let mut adjusted = vec![None; self.length];
        self.visit_population(
            |index| filter[index] >= cutoff,
            |index, pvalue| adjusted[index] = Some(pvalue),
        );
        adjusted
    }

    pub(crate) fn rejections_for_filter(&self, filter: &[f64], cutoff: f64, alpha: f64) -> usize {
        let mut count = 0;
        self.visit_population(
            |index| filter[index] >= cutoff,
            |_, pvalue| count += usize::from(pvalue < alpha),
        );
        count
    }
}

/// Benjamini-Hochberg adjusted p-values for fully observed `f64` p-values.
pub fn bh_adjust_f64(pvalues: &[f64]) -> Result<Vec<f64>, DeseqError> {
    for (idx, value) in pvalues.iter().copied().enumerate() {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(DeseqError::InvalidOptions {
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
                return Err(DeseqError::InvalidOptions {
                    reason: "validated p-value unexpectedly produced a missing adjustment"
                        .to_string(),
                });
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::{PreparedBh, bh_adjust_f64};

    /// Found by cargo-mutants: replacing the count with a constant survived
    /// the whole suite. The count must equal the number of adjusted p-values
    /// below alpha among rows passing the filter.
    #[test]
    fn rejections_for_filter_counts_filtered_adjusted_pvalues_below_alpha() {
        let pvalues = [
            Some(0.001),
            Some(0.04),
            None,
            Some(0.02),
            Some(0.5),
            Some(0.011),
        ];
        let filter = [10.0, 1.0, 5.0, 8.0, 0.5, 3.0];
        let prepared = PreparedBh::new(&pvalues);
        for cutoff in [0.0, 1.0, 3.0, 8.0, 11.0] {
            for alpha in [0.01, 0.05, 0.1] {
                let expected = prepared
                    .adjusted_for_filter(&filter, cutoff)
                    .iter()
                    .filter(|value| value.is_some_and(|p| p < alpha))
                    .count();
                assert_eq!(
                    prepared.rejections_for_filter(&filter, cutoff, alpha),
                    expected,
                    "cutoff {cutoff}, alpha {alpha}"
                );
            }
        }
        // By hand: adjusted 0.005, 0.0275, 0.0333, 0.05, 0.5; three are below 0.05.
        assert_eq!(prepared.rejections_for_filter(&filter, 0.0, 0.05), 3);
        assert_eq!(prepared.rejections_for_filter(&filter, 11.0, 0.05), 0);
        // DESeq2 counts `padj < alpha`: an adjusted p-value equal to alpha is
        // not a rejection (0.025 * 2 is exactly 0.05).
        let boundary = PreparedBh::new(&[Some(0.025), Some(0.05)]);
        assert_eq!(boundary.rejections_for_filter(&[1.0, 1.0], 0.0, 0.05), 0);
    }

    /// Found by cargo-mutants: each invalid p-value must be reported at its
    /// index, including finite values outside [0, 1].
    #[test]
    fn bh_adjust_f64_names_each_invalid_pvalue() {
        for bad in [1.5, -0.1, f64::NAN, f64::INFINITY] {
            let error = bh_adjust_f64(&[0.2, bad]).unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("p-value at index 1 must be finite and within [0, 1]"),
                "{bad}: {error}"
            );
        }
        assert_eq!(bh_adjust_f64(&[0.01, 0.04]).unwrap(), vec![0.02, 0.04]);
    }
}
