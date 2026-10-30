//! Borrowed unsigned counts in native or foreign double storage.
use super::{MatrixOrder, MatrixView};
use crate::core::CountMatrix;
use crate::errors::DeseqError;

#[derive(Clone, Copy, Debug)]
enum Storage<'a> {
    Unsigned(MatrixView<'a, u32>),
    Double(MatrixView<'a, f64>),
}

/// Shape-checked count storage with finite, nonnegative integer u32 values.
/// No matrix is copied or transposed, and double values convert exactly to u32.
#[derive(Clone, Copy, Debug)]
pub struct CountView<'a>(Storage<'a>);

impl<'a> CountView<'a> {
    /// Borrow a shape-checked unsigned matrix without scanning its values.
    pub fn from_unsigned(values: MatrixView<'a, u32>) -> Self {
        Self(Storage::Unsigned(values))
    }

    /// Validate the complete double buffer before borrowing it as counts.
    pub fn from_doubles(values: MatrixView<'a, f64>) -> Result<Self, DeseqError> {
        Self::from_doubles_with(values, |_, _| {})
    }

    // Fuse count validation with a caller's invariant preparation. On failure,
    // callers discard their partial local scratch; no view is returned.
    pub(crate) fn from_doubles_with(
        values: MatrixView<'a, f64>,
        mut visit: impl FnMut(usize, u32),
    ) -> Result<Self, DeseqError> {
        for (index, &value) in values.as_slice().iter().enumerate() {
            if !value.is_finite()
                || value < 0.0
                || value.fract() != 0.0
                || value > f64::from(u32::MAX)
            {
                return Err(DeseqError::InvalidCounts {
                    reason: format!("count at storage index {index} is not a valid u32 count"),
                });
            }
            visit(index, value as u32);
        }
        Ok(Self(Storage::Double(values)))
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        match self.0 {
            Storage::Unsigned(v) => v.n_rows(),
            Storage::Double(v) => v.n_rows(),
        }
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        match self.0 {
            Storage::Unsigned(v) => v.n_cols(),
            Storage::Double(v) => v.n_cols(),
        }
    }

    /// Read a count after checking both logical indices.
    pub fn get(&self, gene: usize, sample: usize) -> Result<u32, DeseqError> {
        match self.0 {
            Storage::Unsigned(v) => v.get(gene, sample).copied(),
            Storage::Double(v) => v.get(gene, sample).map(|&v| v as u32),
        }
    }

    #[inline]
    pub(crate) fn at(&self, gene: usize, sample: usize) -> u32 {
        match self.0 {
            Storage::Unsigned(v) => *v.at(gene, sample),
            Storage::Double(v) => *v.at(gene, sample) as u32,
        }
    }
}

impl<'a> From<&'a CountMatrix> for CountView<'a> {
    fn from(counts: &'a CountMatrix) -> Self {
        Self::from_unsigned(
            MatrixView::new(
                counts.n_genes(),
                counts.n_samples(),
                counts.as_slice(),
                MatrixOrder::RowMajor,
            )
            .expect("CountMatrix has validated shape"),
        )
    }
}
