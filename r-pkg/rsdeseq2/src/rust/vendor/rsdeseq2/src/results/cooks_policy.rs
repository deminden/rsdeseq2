//! Borrowed Cook's masking without table construction or p-value adjustment.
use super::ReplacementRows;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::{CountView, MatrixView};

/// Borrowed two-group low-count heuristic inputs in either storage order.
#[derive(Clone, Copy, Debug)]
pub struct LowCountHeuristic<'a> {
    counts: CountView<'a>,
    cooks: MatrixView<'a, f64>,
}
impl<'a> LowCountHeuristic<'a> {
    /// Validate matching complete count/Cook's matrix dimensions.
    pub fn new(counts: CountView<'a>, cooks: MatrixView<'a, f64>) -> Result<Self, DeseqError> {
        if cooks.n_rows() != counts.n_genes() {
            return Err(invalid_dimensions(
                "Cook's heuristic Cook's rows",
                counts.n_genes(),
                cooks.n_rows(),
            ));
        }
        if cooks.n_cols() != counts.n_samples() {
            return Err(invalid_dimensions(
                "Cook's heuristic Cook's columns",
                counts.n_samples(),
                cooks.n_cols(),
            ));
        }
        Ok(Self { counts, cooks })
    }
    fn spares(self, gene: usize) -> bool {
        // Keep the established finite-only total ordering and last-tie choice.
        let Some((sample, _)) = (0..self.counts.n_samples())
            .map(|i| (i, *self.cooks.at(gene, i)))
            .filter(|(_, v)| v.is_finite())
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
        else {
            return false;
        };
        let outlier_count = self.counts.at(gene, sample);
        (0..self.counts.n_samples())
            .filter(|&i| self.counts.at(gene, i) > outlier_count)
            .take(3)
            .count()
            == 3
    }
}

/// Masked p-values and diagnostic flags; no BH/filtering is performed here.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksMask {
    /// Original or replacement-aware p-values, with outliers made missing.
    pub pvalue: Vec<Option<f64>>,
    /// Missing when filtering is disabled or the maximum Cook's value is absent.
    pub outlier: Vec<Option<bool>>,
}

/// Cook's threshold and optional borrowed heuristic/replacement context.
#[derive(Clone, Copy, Debug)]
pub struct CooksPolicy<'a> {
    cutoff: Option<f64>,
    heuristic: Option<LowCountHeuristic<'a>>,
    replacement: Option<ReplacementRows<'a>>,
}
impl<'a> CooksPolicy<'a> {
    /// Validate a resolved threshold. The low-level legacy API accepts every
    /// finite threshold; high-level request validation requires nonnegative values.
    pub fn new(
        cutoff: Option<f64>,
        heuristic: Option<LowCountHeuristic<'a>>,
    ) -> Result<Self, DeseqError> {
        if let Some(value) = cutoff
            && !value.is_finite()
        {
            return Err(DeseqError::NonFiniteValue {
                context: "Cook's cutoff".into(),
                index: None,
                value,
            });
        }
        Ok(Self {
            cutoff,
            heuristic,
            replacement: None,
        })
    }
    /// Suppress original outlier maxima on refitted/new-zero rows; new-zero
    /// p-values become one before masking. This never changes stored maxima.
    pub fn with_replacement(mut self, rows: ReplacementRows<'a>) -> Self {
        self.replacement = Some(rows);
        self
    }

    /// Evaluate aligned row iterators atomically, validating before a caller
    /// mutates its table. Missingness and strict `maxCook > cutoff` are retained.
    pub fn evaluate<P, M>(&self, pvalues: P, max_cooks: M) -> Result<CooksMask, DeseqError>
    where
        P: IntoIterator<Item = Option<f64>>,
        P::IntoIter: ExactSizeIterator,
        M: IntoIterator<Item = Option<f64>>,
        M::IntoIter: ExactSizeIterator,
    {
        let pvalues = pvalues.into_iter();
        let max_cooks = max_cooks.into_iter();
        let genes = pvalues.len();
        if max_cooks.len() != genes {
            return Err(invalid_dimensions(
                "Cook's maximum rows",
                genes,
                max_cooks.len(),
            ));
        }
        if let Some(heuristic) = self.heuristic
            && heuristic.counts.n_genes() != genes
        {
            return Err(invalid_dimensions(
                "Cook's heuristic result rows",
                genes,
                heuristic.counts.n_genes(),
            ));
        }
        if let Some(rows) = self.replacement
            && rows.genes() != genes
        {
            return Err(invalid_dimensions(
                "Cook's replacement rows",
                genes,
                rows.genes(),
            ));
        }
        let mut result = CooksMask {
            pvalue: Vec::with_capacity(genes),
            outlier: Vec::with_capacity(genes),
        };
        for (gene, (mut pvalue, mut maximum)) in pvalues.zip(max_cooks).enumerate() {
            if let Some(rows) = self.replacement {
                maximum = rows.max_cook(gene, maximum);
                if rows.is_zero(gene) {
                    pvalue = Some(1.0);
                }
            }
            let outlier = if let Some(cutoff) = self.cutoff {
                if let Some(value) = maximum {
                    if !value.is_finite() {
                        return Err(DeseqError::NonFiniteValue {
                            context: "result maxCooks".into(),
                            index: Some(gene),
                            value,
                        });
                    }
                    Some(value > cutoff && !self.heuristic.is_some_and(|h| h.spares(gene)))
                } else {
                    None
                }
            } else {
                None
            };
            if outlier == Some(true) {
                pvalue = None;
            }
            if let Some(p) = pvalue
                && (!p.is_finite() || !(0.0..=1.0).contains(&p))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "result p-value at index {gene} must be finite and within [0, 1]"
                    ),
                });
            }
            result.pvalue.push(pvalue);
            result.outlier.push(outlier);
        }
        Ok(result)
    }
}

/// Test whether a two-level label partition exactly matches a two-column
/// design's second column. Linear time and constant scratch replace an equality outer
/// product; label spelling/order and coefficient scaling do not matter.
pub fn two_level_design_matches<T: Eq>(
    design: MatrixView<'_, f64>,
    labels: &[T],
    declared_levels: usize,
) -> Result<bool, DeseqError> {
    if design.n_cols() != 2 || declared_levels != 2 {
        return Ok(false);
    }
    if labels.len() != design.n_rows() {
        return Err(invalid_dimensions(
            "Cook's factor samples",
            design.n_rows(),
            labels.len(),
        ));
    }
    let Some(first) = labels.first() else {
        return Ok(true);
    };
    let first_value = *design.at(0, 1);
    let mut second = None;
    for (sample, label) in labels.iter().enumerate() {
        let value = *design.at(sample, 1);
        if !value.is_finite() {
            return Ok(false);
        }
        if label == first {
            if value != first_value {
                return Ok(false);
            }
        } else if let Some((other, other_value)) = second {
            if label != other || value != other_value {
                return Ok(false);
            }
        } else {
            if value == first_value {
                return Ok(false);
            }
            second = Some((label, value));
        }
    }
    Ok(true)
}
