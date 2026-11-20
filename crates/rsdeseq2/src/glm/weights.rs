use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions, nonfinite_gene_cell};
use crate::errors::{ErrorKind, Requirement, invalid_parameter};
use crate::math::rank::matrix_rank_in_place;
use crate::matrix::{RowMajorMatrix, try_alloc_matrix};
use crate::numeric::checked;

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
    let samples = design.n_samples();
    let coefficients = design.n_coefficients();
    let mut values = try_alloc_matrix(samples, coefficients, "weighted design")?;
    values.extend_from_slice(design.matrix().as_slice());
    let design_rank =
        matrix_rank_in_place(&mut values, samples, coefficients, options.rank_tolerance)?;
    let mut scratch = if design_rank == coefficients {
        Some(WeightRankScratch {
            values,
            kept_samples: try_alloc_matrix(samples, 1, "Cox-Reid retained samples")?,
            kept_columns: try_alloc_matrix(coefficients, 1, "Cox-Reid retained columns")?,
        })
    } else {
        None
    };
    let mut weights_fail = try_alloc_matrix(weights.n_rows(), 1, "failed weight rows")?;

    for gene in 0..weights.n_rows() {
        let row = normalized_weights.row(gene)?;
        let ok = if row.iter().all(|weight| *weight == 0.0) {
            false
        } else if let Some(scratch) = &mut scratch {
            scratch.weighted_design_full_rank(row, design, options.rank_tolerance)?
                && scratch.cox_reid_subset_full_rank(row, design, gene, options)?
        } else {
            weighted_design_has_no_zero_columns(row, design)
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
    let mut normalized = try_alloc_matrix(
        weights.n_rows(),
        weights.n_cols(),
        "normalized observation weights",
    )?;
    for gene in 0..weights.n_rows() {
        let row = weights.row(gene)?;
        let mut max_weight = 0.0_f64;
        for (sample, weight) in row.iter().copied().enumerate() {
            if !weight.is_finite() || weight < 0.0 {
                return Err(nonfinite_gene_cell(
                    "observation weight",
                    gene,
                    sample,
                    sample,
                    weight,
                ));
            }
            max_weight = max_weight.max(weight);
        }
        if max_weight == 0.0 {
            normalized.extend(std::iter::repeat_n(0.0, weights.n_cols()));
        } else {
            normalized.extend(row.iter().map(|weight| *weight / max_weight));
        }
    }
    RowMajorMatrix::from_row_major(weights.n_rows(), weights.n_cols(), normalized)
}

struct WeightRankScratch {
    values: Vec<f64>,
    kept_samples: Vec<usize>,
    kept_columns: Vec<usize>,
}

impl WeightRankScratch {
    fn weighted_design_full_rank(
        &mut self,
        weights: &[f64],
        design: &DesignMatrix,
        tolerance: f64,
    ) -> Result<bool, DeseqError> {
        let samples = design.n_samples();
        let coefficients = design.n_coefficients();
        self.values.clear();
        for (&weight, row) in weights
            .iter()
            .zip(design.matrix().as_slice().chunks_exact(coefficients))
        {
            // The normalized weight is at most one and design entries are finite.
            self.values.extend(row.iter().map(|value| weight * value));
        }
        Ok(
            matrix_rank_in_place(&mut self.values, samples, coefficients, tolerance)?
                == coefficients,
        )
    }

    fn cox_reid_subset_full_rank(
        &mut self,
        weights: &[f64],
        design: &DesignMatrix,
        gene: usize,
        options: ObservationWeightOptions,
    ) -> Result<bool, DeseqError> {
        self.kept_samples.clear();
        self.kept_samples.extend(
            weights.iter().enumerate().filter_map(|(sample, &weight)| {
                (weight > options.weight_threshold).then_some(sample)
            }),
        );
        if self.kept_samples.is_empty() {
            return Ok(false);
        }

        let coefficients = design.n_coefficients();
        let design_values = design.matrix().as_slice();
        self.kept_columns.clear();
        for column in 0..coefficients {
            let mut col_sum = 0.0;
            for &sample in &self.kept_samples {
                let term = design_values[sample * coefficients + column].abs();
                col_sum = checked::try_add(col_sum, term).map_err(|failure| {
                    DeseqError::CoxReidColumnSumOverflow {
                        gene,
                        sample,
                        column,
                        term,
                        value: failure.value,
                    }
                })?;
            }
            if col_sum > 0.0 {
                self.kept_columns.push(column);
            }
        }
        if self.kept_columns.is_empty() {
            return Ok(false);
        }

        self.values.clear();
        for &sample in &self.kept_samples {
            let start = sample * coefficients;
            self.values.extend(
                self.kept_columns
                    .iter()
                    .map(|&column| design_values[start + column]),
            );
        }
        Ok(matrix_rank_in_place(
            &mut self.values,
            self.kept_samples.len(),
            self.kept_columns.len(),
            options.rank_tolerance,
        )? == self.kept_columns.len())
    }
}

fn weighted_design_has_no_zero_columns(weights: &[f64], design: &DesignMatrix) -> bool {
    let values = design.matrix().as_slice();
    let coefficients = design.n_coefficients();
    (0..coefficients).all(|column| {
        weights
            .iter()
            .enumerate()
            .any(|(sample, &weight)| weight != 0.0 && values[sample * coefficients + column] != 0.0)
    })
}

fn validate_weight_options(options: ObservationWeightOptions) -> Result<(), DeseqError> {
    if !options.weight_threshold.is_finite() || options.weight_threshold < 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Options,
            "weight threshold",
            Requirement::FiniteAndNonNegative,
        ));
    }
    if !options.rank_tolerance.is_finite() || options.rank_tolerance < 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Options,
            "rank tolerance",
            Requirement::FiniteAndNonNegative,
        ));
    }
    Ok(())
}
