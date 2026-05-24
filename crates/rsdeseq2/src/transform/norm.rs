use crate::errors::DeseqError;
use crate::matrix::RowMajorMatrix;

/// Apply DESeq2's `normTransform` to normalized counts.
///
/// This transformation is `log2(q + 1)` for each normalized count `q`. It is a
/// lightweight visualization transform and is not used by the differential
/// expression GLM fit.
pub fn norm_transform(
    normalized_counts: &RowMajorMatrix<f64>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let values = normalized_counts
        .as_slice()
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, count)| norm_transform_value(count, idx))
        .collect::<Result<Vec<_>, _>>()?;
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// Apply DESeq2's `normTransform` to one normalized count.
pub fn norm_transform_value(normalized_count: f64, index: usize) -> Result<f64, DeseqError> {
    if !normalized_count.is_finite() || normalized_count < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "normTransform normalized count".to_string(),
            index: Some(index),
            value: normalized_count,
        });
    }
    Ok((normalized_count + 1.0).log2())
}
