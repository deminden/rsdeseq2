use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::math::qr::matrix_rank;
use crate::matrix::RowMajorMatrix;

/// Options for DESeq2-style observation-weight preprocessing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservationWeightOptions {
    /// Threshold used by DESeq2 when checking the Cox-Reid sub-design rank.
    pub weight_threshold: f64,
    /// Absolute tolerance for deterministic matrix-rank checks.
    pub rank_tolerance: f64,
}

impl Default for ObservationWeightOptions {
    fn default() -> Self {
        Self {
            weight_threshold: 1e-2,
            rank_tolerance: 1e-10,
        }
    }
}

/// Output of DESeq2-style observation-weight preprocessing.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservationWeights {
    /// Row-normalized weights, genes x samples.
    pub weights: RowMajorMatrix<f64>,
    /// Rows whose weights fail to allow parameter estimation.
    pub weights_fail: Vec<bool>,
    /// Rank of the unweighted design matrix under the configured tolerance.
    pub design_rank: usize,
}

impl ObservationWeights {
    /// Per-row flags indicating usable weights.
    pub fn weights_ok(&self) -> Vec<bool> {
        self.weights_fail.iter().map(|fail| !fail).collect()
    }
}

/// Normalize and check observation weights following DESeq2's `getAndCheckWeights` shape.
///
/// The input matrix is genes x samples. Each gene row is divided by its row
/// maximum before rank checks. Rows that cannot support parameter estimation
/// are marked in `weights_fail`; callers can then treat those rows like DESeq2
/// treats `mcols(dds)$allZero = TRUE` for failed weights.
pub fn preprocess_observation_weights(
    weights: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
) -> Result<ObservationWeights, DeseqError> {
    preprocess_observation_weights_with_options(
        weights,
        design,
        ObservationWeightOptions::default(),
    )
}

/// Normalize and check observation weights with explicit options.
pub fn preprocess_observation_weights_with_options(
    weights: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
    options: ObservationWeightOptions,
) -> Result<ObservationWeights, DeseqError> {
    validate_weight_options(options)?;
    if weights.n_cols() != design.n_samples() {
        return Err(invalid_dimensions(
            "observation weight columns",
            design.n_samples(),
            weights.n_cols(),
        ));
    }

    let normalized_weights = normalize_weight_rows(weights)?;
    let design_values = design.matrix().as_slice();
    let p = design.n_coefficients();
    let design_rank = matrix_rank(design_values, design.n_samples(), p, options.rank_tolerance)?;
    let design_full_rank = design_rank == p;
    let mut weights_fail = Vec::with_capacity(weights.n_rows());

    for gene in 0..weights.n_rows() {
        let row = normalized_weights.row(gene)?;
        let ok = if row.iter().all(|weight| *weight == 0.0) {
            false
        } else if design_full_rank {
            weighted_design_full_rank(row, design, options.rank_tolerance)?
                && cox_reid_subset_full_rank(
                    row,
                    design,
                    options.weight_threshold,
                    options.rank_tolerance,
                )?
        } else {
            weighted_design_has_no_zero_columns(row, design)?
        };
        weights_fail.push(!ok);
    }

    Ok(ObservationWeights {
        weights: normalized_weights,
        weights_fail,
        design_rank,
    })
}

fn normalize_weight_rows(weights: &RowMajorMatrix<f64>) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut normalized = Vec::with_capacity(weights.len());
    for gene in 0..weights.n_rows() {
        let row = weights.row(gene)?;
        let mut max_weight = 0.0_f64;
        for (sample, weight) in row.iter().copied().enumerate() {
            validate_nonnegative_finite(weight, "observation weight", sample)?;
            max_weight = max_weight.max(weight);
        }
        if max_weight == 0.0 {
            normalized.extend(std::iter::repeat(0.0).take(weights.n_cols()));
        } else {
            normalized.extend(row.iter().map(|weight| *weight / max_weight));
        }
    }
    RowMajorMatrix::from_row_major(weights.n_rows(), weights.n_cols(), normalized)
}

fn weighted_design_full_rank(
    weights: &[f64],
    design: &DesignMatrix,
    tolerance: f64,
) -> Result<bool, DeseqError> {
    let values = weighted_design_values(weights, design)?;
    Ok(matrix_rank(
        &values,
        design.n_samples(),
        design.n_coefficients(),
        tolerance,
    )? == design.n_coefficients())
}

fn cox_reid_subset_full_rank(
    weights: &[f64],
    design: &DesignMatrix,
    weight_threshold: f64,
    tolerance: f64,
) -> Result<bool, DeseqError> {
    let kept_rows = weights
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(sample, weight)| (weight > weight_threshold).then_some(sample))
        .collect::<Vec<_>>();
    if kept_rows.is_empty() {
        return Ok(false);
    }

    let mut keep_cols = Vec::new();
    for col in 0..design.n_coefficients() {
        let mut col_sum = 0.0;
        for sample in &kept_rows {
            let value = design
                .matrix()
                .get(*sample, col)
                .copied()
                .unwrap_or(0.0)
                .abs();
            let Some(next_sum) = checked_sum2(col_sum, value) else {
                return Err(DeseqError::NonFiniteValue {
                    context: "observation weight Cox-Reid column sum".to_string(),
                    index: Some(col),
                    value,
                });
            };
            col_sum = next_sum;
        }
        if col_sum > 0.0 {
            keep_cols.push(col);
        }
    }
    if keep_cols.is_empty() {
        return Ok(false);
    }

    let mut values = Vec::with_capacity(kept_rows.len() * keep_cols.len());
    for sample in &kept_rows {
        for col in &keep_cols {
            values.push(design.matrix().get(*sample, *col).copied().unwrap_or(0.0));
        }
    }
    Ok(matrix_rank(&values, kept_rows.len(), keep_cols.len(), tolerance)? == keep_cols.len())
}

fn weighted_design_has_no_zero_columns(
    weights: &[f64],
    design: &DesignMatrix,
) -> Result<bool, DeseqError> {
    for col in 0..design.n_coefficients() {
        let mut has_nonzero = false;
        for (sample, weight) in weights.iter().copied().enumerate() {
            let value = design.matrix().get(sample, col).copied().unwrap_or(0.0);
            if weight != 0.0 && value != 0.0 {
                has_nonzero = true;
                break;
            }
        }
        if !has_nonzero {
            return Ok(false);
        }
    }
    Ok(true)
}

fn weighted_design_values(weights: &[f64], design: &DesignMatrix) -> Result<Vec<f64>, DeseqError> {
    if weights.len() != design.n_samples() {
        return Err(invalid_dimensions(
            "observation weight row",
            design.n_samples(),
            weights.len(),
        ));
    }
    let mut values = Vec::with_capacity(design.n_samples() * design.n_coefficients());
    for (sample, weight) in weights.iter().copied().enumerate() {
        for col in 0..design.n_coefficients() {
            let value = weight * design.matrix().get(sample, col).copied().unwrap_or(0.0);
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "weighted design value".to_string(),
                    index: Some(sample * design.n_coefficients() + col),
                    value,
                });
            }
            values.push(value);
        }
    }
    Ok(values)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn validate_weight_options(options: ObservationWeightOptions) -> Result<(), DeseqError> {
    if !options.weight_threshold.is_finite() || options.weight_threshold < 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "weight threshold must be finite and non-negative".to_string(),
        });
    }
    if !options.rank_tolerance.is_finite() || options.rank_tolerance < 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rank tolerance must be finite and non-negative".to_string(),
        });
    }
    Ok(())
}

fn validate_nonnegative_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}
