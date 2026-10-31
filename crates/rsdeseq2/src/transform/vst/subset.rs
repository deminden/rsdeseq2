use super::*;

/// Row-aligned subset for fast-VST trend fitting.
///
/// The count matrix, normalized counts, normalization factors, and observation
/// weights are all returned in the same deterministic row order selected by
/// [`fast_vst_subset_indices`].
#[derive(Clone, Debug, PartialEq)]
pub struct FastVstSubset {
    /// Zero-based row indices into the original full dataset.
    pub row_indices: Vec<usize>,
    /// Raw count subset in fast-VST order.
    pub counts: CountMatrix,
    /// Normalized count subset in fast-VST order.
    pub normalized_counts: RowMajorMatrix<f64>,
    /// Optional normalization-factor subset in fast-VST order.
    pub normalization_factors: Option<RowMajorMatrix<f64>>,
    /// Optional observation-weight subset in fast-VST order.
    pub observation_weights: Option<RowMajorMatrix<f64>>,
}

/// Metadata summary for a row-aligned fast-VST subset bundle.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastVstSubsetMetadata {
    /// Number of rows selected for the deterministic fast subset.
    pub rows: usize,
    /// Number of samples in the selected matrices.
    pub cols: usize,
    /// Original zero-based row indices selected for the fast subset.
    pub row_indices: Vec<usize>,
    /// Whether normalization factors were included in the subset.
    pub has_normalization_factors: bool,
    /// Whether observation weights were included in the subset.
    pub has_observation_weights: bool,
}

impl FastVstSubset {
    /// Metadata for the deterministic fast-VST subset.
    pub fn metadata(&self) -> FastVstSubsetMetadata {
        FastVstSubsetMetadata {
            rows: self.counts.n_genes(),
            cols: self.counts.n_samples(),
            row_indices: self.row_indices.clone(),
            has_normalization_factors: self.normalization_factors.is_some(),
            has_observation_weights: self.observation_weights.is_some(),
        }
    }
}

/// Select rows using the deterministic fast-VST rule.
///
/// With ordinary normalized-count row means, this matches DESeq2 1.52.0: keep
/// values above 5, order by those means, and take
/// `round(seq(1, n, length.out=nsub))` positions on the ordered one-based
/// index. Callers may supply another finite selection statistic, such as the
/// weighted `baseMean` used by Rust's weighted builder path. The returned
/// indices are zero-based row indices into the original matrix.
pub fn fast_vst_subset_indices(base_mean: &[f64], nsub: usize) -> Result<Vec<usize>, DeseqError> {
    if nsub == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "fast VST subset size must be positive".to_string(),
        });
    }
    let mut eligible = fast_vst_eligible_rows(base_mean)?;
    if eligible.len() < nsub {
        return Err(DeseqError::InvalidDimensions {
            context: "fast VST rows with mean normalized count above 5".to_string(),
            expected: nsub,
            actual: eligible.len(),
        });
    }
    eligible.sort_by(|left, right| {
        left.1
            .total_cmp(&right.1)
            .then_with(|| left.0.cmp(&right.0))
    });

    let last = eligible.len();
    let positions = if nsub == 1 {
        vec![1usize]
    } else {
        (0..nsub)
            .map(|idx| {
                let value = 1.0 + (last as f64 - 1.0) * idx as f64 / (nsub as f64 - 1.0);
                round_half_to_even(value) as usize
            })
            .collect()
    };
    Ok(positions
        .into_iter()
        .map(|position| eligible[position - 1].0)
        .collect())
}

/// Count rows eligible under the deterministic fast-VST rule.
///
/// Rows are eligible when their base mean is finite and greater than 5. This
/// helper validates finite input with the same checks used by
/// [`fast_vst_subset_indices`], which lets callers detect whether the default
/// `nsub` can be used before requesting the subset.
pub fn fast_vst_eligible_count(base_mean: &[f64]) -> Result<usize, DeseqError> {
    Ok(fast_vst_eligible_rows(base_mean)?.len())
}

fn fast_vst_eligible_rows(base_mean: &[f64]) -> Result<Vec<(usize, f64)>, DeseqError> {
    let mut eligible = Vec::new();
    for (idx, mean) in base_mean.iter().copied().enumerate() {
        if !mean.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "fast VST base mean".to_string(),
                index: Some(idx),
                value: mean,
            });
        }
        if mean > 5.0 {
            eligible.push((idx, mean));
        }
    }
    Ok(eligible)
}

/// Select normalized-count rows for a fast-VST trend fit.
///
/// This applies [`fast_vst_subset_indices`] to row base means and returns the
/// selected rows in deterministic subset order.
pub fn fast_vst_subset_normalized_counts(
    normalized_counts: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    nsub: usize,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    fast_vst_subset_matrix_rows(normalized_counts, base_mean, nsub, "fast VST base means")
}

/// Select rows from any gene/sample matrix aligned to the fast `vst()` subset.
///
/// This is useful for normalization factors or observation weights that must
/// stay aligned with the subset count matrix used for fast-VST trend fitting.
pub fn fast_vst_subset_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    nsub: usize,
    context: &str,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if base_mean.len() != matrix.n_rows() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: matrix.n_rows(),
            actual: base_mean.len(),
        });
    }
    let row_indices = fast_vst_subset_indices(base_mean, nsub)?;
    select_matrix_rows(matrix, &row_indices)
}

/// Build a complete row-aligned fast-VST input bundle.
///
/// This helper centralizes the subset rule so raw counts and optional
/// gene/sample matrices cannot drift out of alignment with the normalized
/// counts used to fit the dispersion trend.
pub fn fast_vst_subset(
    counts: &CountMatrix,
    normalized_counts: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    nsub: usize,
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    observation_weights: Option<&RowMajorMatrix<f64>>,
) -> Result<FastVstSubset, DeseqError> {
    validate_fast_vst_matrix_shape(normalized_counts, counts, "fast VST normalized counts")?;
    let row_indices = fast_vst_subset_indices(base_mean, nsub)?;
    let subset_counts = counts.select_rows(&row_indices)?;
    let subset_normalized = select_matrix_rows(normalized_counts, &row_indices)?;
    let subset_factors = match normalization_factors {
        Some(factors) => {
            validate_fast_vst_matrix_shape(factors, counts, "fast VST normalization factors")?;
            Some(select_matrix_rows(factors, &row_indices)?)
        }
        None => None,
    };
    let subset_weights = match observation_weights {
        Some(weights) => {
            validate_fast_vst_matrix_shape(weights, counts, "fast VST observation weights")?;
            Some(select_matrix_rows(weights, &row_indices)?)
        }
        None => None,
    };
    Ok(FastVstSubset {
        row_indices,
        counts: subset_counts,
        normalized_counts: subset_normalized,
        normalization_factors: subset_factors,
        observation_weights: subset_weights,
    })
}

fn validate_fast_vst_matrix_shape(
    matrix: &RowMajorMatrix<f64>,
    counts: &CountMatrix,
    context: &str,
) -> Result<(), DeseqError> {
    if matrix.n_rows() != counts.n_genes() || matrix.n_cols() != counts.n_samples() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: matrix.len(),
        });
    }
    Ok(())
}

fn select_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    row_indices: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if row_indices.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "selected matrix rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let mut values = Vec::with_capacity(row_indices.len() * matrix.n_cols());
    for row in row_indices {
        values.extend_from_slice(matrix.row(*row)?);
    }
    RowMajorMatrix::from_row_major(row_indices.len(), matrix.n_cols(), values)
}
