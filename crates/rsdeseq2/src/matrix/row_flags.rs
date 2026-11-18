//! Borrowed per-row boolean flags.

use crate::errors::DeseqError;

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
    pub(crate) fn from_indices(genes: usize, rows: &'a [usize]) -> Result<Self, DeseqError> {
        if rows.iter().any(|&row| row >= genes) || rows.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(DeseqError::InvalidOptions {
                reason: "replacement row indices must be in range and strictly increasing".into(),
            });
        }
        Ok(Self(Flags::Indices { genes, rows }))
    }
    #[doc(hidden)]
    #[allow(clippy::len_without_is_empty)] // the gene count, not a collection length
    pub fn len(self) -> usize {
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
