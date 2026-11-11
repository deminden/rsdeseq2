//! Replacement-row presentation shared by standalone and wrapped results.
use crate::errors::{DeseqError, invalid_dimensions};

#[derive(Clone, Copy, Debug)]
enum Flags<'a> {
    Boolean(&'a [bool]),
    Integer(&'a [i32]),
    Indices { genes: usize, rows: &'a [usize] },
}

/// Borrowed, validated row indicators without a foreign-to-native flag copy.
#[derive(Clone, Copy, Debug)]
pub struct RowFlags<'a>(Flags<'a>);
impl<'a> RowFlags<'a> {
    /// Borrow native boolean flags.
    pub fn from_bools(values: &'a [bool]) -> Self {
        Self(Flags::Boolean(values))
    }
    /// Borrow integer flags, rejecting any value other than zero or one.
    pub fn from_integers(values: &'a [i32]) -> Result<Self, DeseqError> {
        if let Some(index) = values.iter().position(|&v| v != 0 && v != 1) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("row flag at index {index} must be zero or one"),
            });
        }
        Ok(Self(Flags::Integer(values)))
    }
    fn from_indices(genes: usize, rows: &'a [usize]) -> Result<Self, DeseqError> {
        if rows.iter().any(|&row| row >= genes) || rows.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(DeseqError::InvalidOptions {
                reason: "replacement row indices must be in range and strictly increasing".into(),
            });
        }
        Ok(Self(Flags::Indices { genes, rows }))
    }
    pub(crate) fn len(self) -> usize {
        match self.0 {
            Flags::Boolean(v) => v.len(),
            Flags::Integer(v) => v.len(),
            Flags::Indices { genes, .. } => genes,
        }
    }
    pub(crate) fn at(self, gene: usize) -> bool {
        match self.0 {
            Flags::Boolean(v) => v[gene],
            Flags::Integer(v) => v[gene] == 1,
            Flags::Indices { rows, .. } => rows.binary_search(&gene).is_ok(),
        }
    }
}

/// Row context for successful refits and rows made all-zero by replacement.
#[derive(Clone, Copy, Debug)]
pub struct ReplacementRows<'a> {
    genes: usize,
    refitted: Option<RowFlags<'a>>,
    new_zero: Option<RowFlags<'a>>,
}
impl<'a> ReplacementRows<'a> {
    /// Borrow the sorted row indices retained by a standalone replacement plan.
    /// No gene-length flag vectors or count matrices are copied.
    pub fn from_plan(plan: &'a crate::cooks::CooksRefitPlan) -> Result<Self, DeseqError> {
        let genes = plan.replacement.replaced_counts.n_genes();
        if !plan.has_replaceable_samples() {
            return Self::new(genes, None, None);
        }
        Self::new(
            genes,
            Some(RowFlags::from_indices(genes, &plan.refit_rows)?),
            Some(RowFlags::from_indices(genes, &plan.new_all_zero_rows)?),
        )
    }
    /// Validate every optional indicator length before applying row policy.
    pub fn new(
        genes: usize,
        refitted: Option<RowFlags<'a>>,
        new_zero: Option<RowFlags<'a>>,
    ) -> Result<Self, DeseqError> {
        for (name, flags) in [
            ("refitted result rows", refitted),
            ("replacement zero rows", new_zero),
        ] {
            if let Some(flags) = flags
                && flags.len() != genes
            {
                return Err(invalid_dimensions(name, genes, flags.len()));
            }
        }
        Ok(Self {
            genes,
            refitted,
            new_zero,
        })
    }
    pub(super) fn genes(self) -> usize {
        self.genes
    }
    /// Whether effect/SE/stat columns need writable replacement storage.
    pub fn has_new_zero(self) -> bool {
        (0..self.genes).any(|gene| self.is_zero(gene))
    }
    pub(super) fn is_zero(self, gene: usize) -> bool {
        self.new_zero.is_some_and(|v| v.at(gene))
    }
    pub(super) fn max_cook(self, gene: usize, original: Option<f64>) -> Option<f64> {
        if self.is_zero(gene) || self.refitted.is_some_and(|v| v.at(gene)) {
            Some(0.0)
        } else {
            original
        }
    }

    /// Set the established new-zero presentation values: effect/SE/stat/p = 0/0/0/1.
    /// Both native optional columns and foreign double columns use this helper.
    pub fn set_zero_values<T: From<f64>>(effect: &mut T, se: &mut T, stat: &mut T, pvalue: &mut T) {
        *effect = 0.0.into();
        *se = 0.0.into();
        *stat = 0.0.into();
        *pvalue = 1.0.into();
    }

    /// Apply overrides after validating all four column lengths. Other values,
    /// including foreign missing-value payloads, remain untouched.
    pub fn apply_columns<T: From<f64>>(
        &self,
        effect: &mut [T],
        se: &mut [T],
        stat: &mut [T],
        pvalue: &mut [T],
    ) -> Result<(), DeseqError> {
        for (name, length) in [
            ("replacement effects", effect.len()),
            ("replacement SEs", se.len()),
            ("replacement statistics", stat.len()),
            ("replacement p-values", pvalue.len()),
        ] {
            if length != self.genes {
                return Err(invalid_dimensions(name, self.genes, length));
            }
        }
        for gene in 0..self.genes {
            if self.is_zero(gene) {
                Self::set_zero_values(
                    &mut effect[gene],
                    &mut se[gene],
                    &mut stat[gene],
                    &mut pvalue[gene],
                );
            }
        }
        Ok(())
    }
}
