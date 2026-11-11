//! Canonical normal-mixture capture layout and component recycling, in Rust.
use super::BuiltinLocationLikelihood;
use crate::{
    DeseqError,
    matrix::{MatrixOrder, MatrixView},
};
fn invalid(s: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("normal-mixture captures: {s}"),
    }
}
fn columns(m: MatrixView<'_, f64>) -> Vec<f64> {
    (0..m.n_cols())
        .flat_map(|j| (0..m.n_rows()).map(move |i| *m.at(i, j)))
        .collect()
}
pub(super) fn import(
    weights: MatrixView<'_, f64>,
    scales: MatrixView<'_, f64>,
    scales_are_matrix: bool,
) -> Result<BuiltinLocationLikelihood, DeseqError> {
    // R iterates atomic captures in column order, independent of their dimensions.
    let weights = columns(weights);
    if weights.len() == scales.as_slice().len() {
        if scales.n_rows() == 1 {
            return BuiltinLocationLikelihood::normal_mixture_by_observation(weights, scales);
        }
        // Equal total lengths flatten the SD capture into one global mixture.
        return BuiltinLocationLikelihood::normal_mixture(weights, columns(scales));
    }
    if !scales_are_matrix {
        return Err(invalid("unequal capture lengths require an SD matrix"));
    }
    let k = scales.n_cols();
    if k == 0 || weights.is_empty() {
        return Err(invalid("nonempty weights and SD columns are required"));
    }
    if weights.len() == k {
        return BuiltinLocationLikelihood::normal_mixture_by_observation(weights, scales);
    }
    // mapply recycles component lists, retaining operation order; do not merge
    // equal SD components or combine weights, which changes floating-point sums.
    let components = weights.len().max(k);
    let len = scales
        .n_rows()
        .checked_mul(components)
        .ok_or_else(|| invalid("expanded capture size overflows"))?;
    let mut values = Vec::with_capacity(len);
    for i in 0..scales.n_rows() {
        for j in 0..components {
            values.push(*scales.at(i, j % k));
        }
    }
    let expanded = (0..components)
        .map(|j| weights[j % weights.len()])
        .collect();
    BuiltinLocationLikelihood::normal_mixture_by_observation(
        expanded,
        MatrixView::new(scales.n_rows(), components, &values, MatrixOrder::RowMajor)?,
    )
}
