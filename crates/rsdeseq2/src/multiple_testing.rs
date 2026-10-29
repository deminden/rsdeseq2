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
