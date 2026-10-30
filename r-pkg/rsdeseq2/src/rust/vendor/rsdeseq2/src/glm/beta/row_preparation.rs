//! Shared raw-count row selection before normal-prior refitting.
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::CountRows;
use crate::results::RowFlags;

/// Nonzero rows to marshal into a normal-prior fit, in original row order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalShrinkageRows {
    indices: Vec<usize>,
}
impl NormalShrinkageRows {
    /// Validate final/fitted dispersions on raw nonzero rows and select them.
    /// Invalid values on raw-zero or explicitly excluded fitted rows are unused.
    /// Keep raw zeros distinct from fitted exclusions such as weight failures;
    /// a nonzero excluded row must also stay out of prior estimation.
    pub fn new(
        rows: &CountRows<'_>,
        dispersions: &[f64],
        fitted: &[f64],
        excluded: Option<RowFlags<'_>>,
    ) -> Result<Self, DeseqError> {
        let genes = rows.counts().n_genes();
        for (name, values) in [
            ("normal shrinkage dispersions", dispersions),
            ("normal shrinkage fitted trends", fitted),
        ] {
            if values.len() != genes {
                return Err(invalid_dimensions(name, genes, values.len()));
            }
        }
        if let Some(flags) = excluded
            && flags.len() != genes
        {
            return Err(invalid_dimensions(
                "normal shrinkage excluded rows",
                genes,
                flags.len(),
            ));
        }
        let mut indices = Vec::new();
        for (gene, &zero) in rows.all_zero().iter().enumerate() {
            if zero || excluded.is_some_and(|flags| flags.at(gene)) {
                continue;
            }
            if !dispersions[gene].is_finite()
                || dispersions[gene] <= 0.0
                || !fitted[gene].is_finite()
                || fitted[gene] <= 0.0
            {
                return Err(DeseqError::InvalidOptions {reason:"normal lfcShrink requires finite positive dispersions and fitted trends for non-zero genes".into()});
            }
            indices.push(gene);
        }
        if indices.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "normal lfcShrink cannot estimate a beta prior when every gene is all zero or excluded"
                    .into(),
            });
        }
        Ok(Self { indices })
    }
    /// Zero-based source rows, in ascending input order.
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }
}
