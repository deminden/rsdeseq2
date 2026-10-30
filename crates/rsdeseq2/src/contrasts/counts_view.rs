//! Borrowed counts and reusable count-zero invariants for contrast requests.
use super::all_zero::{factor_level_samples, numeric_samples};
use crate::core::CountMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::{MatrixOrder, MatrixView};

#[derive(Clone, Copy, Debug)]
enum Storage<'a> {
    Unsigned(MatrixView<'a, u32>),
    Double(MatrixView<'a, f64>),
}

/// Validated borrowed counts with one cached global-zero flag per gene.
///
/// Matrix storage is never copied. Validation and global-zero detection happen
/// once; the same view can evaluate multiple numeric or factor contrasts.
#[derive(Debug)]
pub struct ContrastCounts<'a> {
    storage: Storage<'a>,
    all_zero: Vec<bool>,
}

impl<'a> ContrastCounts<'a> {
    /// Borrow unsigned counts and cache global-zero rows.
    pub fn from_unsigned(values: MatrixView<'a, u32>) -> Self {
        let mut all_zero = vec![true; values.n_rows()];
        for (index, value) in values.as_slice().iter().enumerate() {
            let gene = match values.order() {
                MatrixOrder::RowMajor => index / values.n_cols(),
                MatrixOrder::ColumnMajor => index % values.n_rows(),
            };
            all_zero[gene] &= *value == 0;
        }
        Self {
            storage: Storage::Unsigned(values),
            all_zero,
        }
    }

    /// Borrow doubles representing finite, nonnegative integer u32 counts.
    /// Validate the complete buffer and cache global-zero rows in one pass.
    pub fn from_doubles(values: MatrixView<'a, f64>) -> Result<Self, DeseqError> {
        let mut all_zero = vec![true; values.n_rows()];
        for (index, &value) in values.as_slice().iter().enumerate() {
            if !value.is_finite()
                || value < 0.0
                || value.fract() != 0.0
                || value > f64::from(u32::MAX)
            {
                return Err(DeseqError::InvalidCounts {
                    reason: format!(
                        "contrast count at storage index {index} is not a valid u32 count"
                    ),
                });
            }
            let gene = match values.order() {
                MatrixOrder::RowMajor => index / values.n_cols(),
                MatrixOrder::ColumnMajor => index % values.n_rows(),
            };
            all_zero[gene] &= value == 0.0;
        }
        Ok(Self {
            storage: Storage::Double(values),
            all_zero,
        })
    }

    /// Gene count.
    pub fn n_genes(&self) -> usize {
        self.all_zero.len()
    }

    /// Sample count.
    pub fn n_samples(&self) -> usize {
        match self.storage {
            Storage::Unsigned(v) => v.n_cols(),
            Storage::Double(v) => v.n_cols(),
        }
    }

    /// Cached flags for genes with zero counts in every sample.
    pub fn all_zero(&self) -> &[bool] {
        &self.all_zero
    }

    /// Numeric-contrast zero flags, including globally zero rows when the
    /// contrast has both signs. One-sided requests return all false and do not
    /// require a model matrix. Mixed-sign requests require aligned design data.
    pub fn numeric_flags(
        &self,
        design: Option<MatrixView<'_, f64>>,
        contrast: &[f64],
    ) -> Result<Vec<bool>, DeseqError> {
        match numeric_samples(self.n_samples(), design, contrast)? {
            Some(selected) => self.selected_flags(&selected),
            None => Ok(vec![false; self.n_genes()]),
        }
    }

    /// Factor-contrast zero flags using the existing native alias policy.
    pub fn factor_flags<S: AsRef<str>>(
        &self,
        sample_levels: &[S],
        numerator: &str,
        denominator: &str,
    ) -> Result<Vec<bool>, DeseqError> {
        let selected =
            factor_level_samples(self.n_samples(), sample_levels, numerator, denominator)?;
        self.selected_flags(&selected)
    }

    pub(super) fn selected_flags(&self, selected: &[bool]) -> Result<Vec<bool>, DeseqError> {
        if selected.len() != self.n_samples() {
            return Err(invalid_dimensions(
                "contrastAllZero selected samples",
                self.n_samples(),
                selected.len(),
            ));
        }
        Ok((0..self.n_genes())
            .map(|gene| {
                selected.iter().enumerate().all(|(sample, &selected)| {
                    !selected
                        || match self.storage {
                            Storage::Unsigned(values) => *values.at(gene, sample) == 0,
                            Storage::Double(values) => *values.at(gene, sample) == 0.0,
                        }
                })
            })
            .collect())
    }
}

impl<'a> From<&'a CountMatrix> for ContrastCounts<'a> {
    fn from(counts: &'a CountMatrix) -> Self {
        // CountMatrix's private constructors guarantee this complete shape.
        let matrix = MatrixView::new(
            counts.n_genes(),
            counts.n_samples(),
            counts.as_slice(),
            MatrixOrder::RowMajor,
        )
        .expect("validated CountMatrix shape");
        Self::from_unsigned(matrix)
    }
}
