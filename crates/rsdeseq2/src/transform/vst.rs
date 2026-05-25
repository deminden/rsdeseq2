use crate::core::CountMatrix;
use crate::dispersion::{DispersionTrendFit, LocalDispersionTrend, ParametricDispersionTrend};
use crate::errors::DeseqError;
use crate::matrix::RowMajorMatrix;

/// Row-aligned subset used by DESeq2's fast `vst()` trend fit.
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
    /// DESeq2-shaped metadata view for the deterministic fast-VST subset.
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

/// Select rows for DESeq2's fast `vst()` dispersion-trend subset.
///
/// DESeq2 first keeps rows with mean normalized count above 5, orders those
/// rows by base mean, then takes `round(seq(1, n, length.out=nsub))` positions
/// on the ordered one-based index. The returned indices are zero-based row
/// indices into the original matrix.
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

/// Count rows eligible for DESeq2's fast `vst()` trend-fit subset.
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

/// Select normalized-count rows for DESeq2's fast `vst()` trend fit.
///
/// This applies [`fast_vst_subset_indices`] to row base means and returns the
/// selected rows in the same deterministic order DESeq2 uses for the subset
/// dataset passed to dispersion fitting.
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

/// Build the complete row-aligned input bundle for DESeq2's fast `vst()`.
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

/// Apply DESeq2's mean-fit variance-stabilizing transformation.
///
/// This is the closed-form fixed-dispersion branch used by DESeq2 when the
/// dispersion fit type is `mean`:
///
/// `log2(exp(2 * asinh(sqrt(alpha * q))) / (4 * alpha))`
///
/// where `q` is a normalized count and `alpha` is the mean dispersion.
pub fn vst_mean(
    normalized_counts: &RowMajorMatrix<f64>,
    mean_dispersion: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_mean_dispersion(mean_dispersion)?;
    let values = normalized_counts
        .as_slice()
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, count)| vst_mean_value(count, mean_dispersion, idx))
        .collect::<Result<Vec<_>, _>>()?;
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// Apply DESeq2's parametric-trend variance-stabilizing transformation.
///
/// This is the closed-form branch used by DESeq2 when the dispersion function
/// has `fitType="parametric"` and coefficients `asymptDisp` and `extraPois`.
pub fn vst_parametric(
    normalized_counts: &RowMajorMatrix<f64>,
    trend: ParametricDispersionTrend,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_parametric_trend(trend)?;
    let values = normalized_counts
        .as_slice()
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, count)| vst_parametric_value(count, trend, idx))
        .collect::<Result<Vec<_>, _>>()?;
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// Apply the local-trend VST numerical-integration branch.
///
/// DESeq2 uses this branch for `fitType="local"` by integrating the reciprocal
/// square root of the normalized-count variance curve, then rescaling the
/// transform so high counts follow `log2(q)`. `inverse_size_factor_mean`
/// corresponds to `mean(1 / sizeFactors)` or, for normalization factors, the
/// analogous approximation from their column geometric means.
pub fn vst_local(
    normalized_counts: &RowMajorMatrix<f64>,
    trend: &LocalDispersionTrend,
    inverse_size_factor_mean: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_inverse_size_factor_mean(inverse_size_factor_mean)?;
    let max_count = validate_normalized_counts_and_max(normalized_counts)?;
    if max_count <= 0.0 {
        return RowMajorMatrix::from_elem(
            normalized_counts.n_rows(),
            normalized_counts.n_cols(),
            f64::NAN,
        );
    }
    let row_means = normalized_count_row_means(normalized_counts)?;
    let h1 = quantile_type7(&row_means, 0.95)?;
    let h2 = quantile_type7(&row_means, 0.999)?;
    if h1 <= 0.0 || h2 <= h1 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local VST scaling quantiles must be positive and increasing".to_string(),
        });
    }
    let integral = LocalVstIntegral::fit(max_count, trend, inverse_size_factor_mean)?;
    let int_h1 = integral.evaluate(h1.asinh())?;
    let int_h2 = integral.evaluate(h2.asinh())?;
    if int_h2 <= int_h1 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local VST integral must increase across scaling quantiles".to_string(),
        });
    }
    let eta = (h2.log2() - h1.log2()) / (int_h2 - int_h1);
    let xi = h1.log2() - eta * int_h1;
    let values = normalized_counts
        .as_slice()
        .iter()
        .copied()
        .map(|count| {
            integral
                .evaluate(count.asinh())
                .map(|value| eta * value + xi)
        })
        .collect::<Result<Vec<_>, _>>()?;
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// Compute the local-VST size-factor summary from sample size factors.
///
/// This is DESeq2's `mean(1 / sizeFactors)` term used in the local VST
/// variance curve.
pub fn local_vst_inverse_size_factor_mean(size_factors: &[f64]) -> Result<f64, DeseqError> {
    if size_factors.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "local VST size factors".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let mut sum = 0.0;
    for (idx, size_factor) in size_factors.iter().copied().enumerate() {
        if !size_factor.is_finite() || size_factor <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("local VST size factor at index {idx} must be finite and positive"),
            });
        }
        sum = checked_add(
            sum,
            size_factor.recip(),
            idx,
            "local VST inverse size-factor sum",
        )?;
    }
    Ok(sum / size_factors.len() as f64)
}

/// Compute the local-VST size-factor summary from normalization factors.
///
/// This mirrors DESeq2's approximation for the local VST branch when
/// normalization factors are present:
/// `sf = exp(colMeans(log(normalizationFactors)))`, then `mean(1 / sf)`.
pub fn local_vst_inverse_size_factor_mean_from_normalization_factors(
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<f64, DeseqError> {
    let mut log_col_sums = vec![0.0; normalization_factors.n_cols()];
    for row in 0..normalization_factors.n_rows() {
        for (sample, factor) in normalization_factors.row(row)?.iter().copied().enumerate() {
            if !factor.is_finite() || factor <= 0.0 {
                return Err(DeseqError::InvalidSizeFactors {
                    reason: format!(
                        "local VST normalization factor at index {} must be finite and positive",
                        row * normalization_factors.n_cols() + sample
                    ),
                });
            }
            checked_add_assign(
                &mut log_col_sums[sample],
                factor.ln(),
                row * normalization_factors.n_cols() + sample,
                "local VST normalization-factor log column sum",
            )?;
        }
    }
    let mut inverse_sum = 0.0;
    for (sample, sum) in log_col_sums.into_iter().enumerate() {
        let column_geometric_mean = (sum / normalization_factors.n_rows() as f64).exp();
        if !column_geometric_mean.is_finite() || column_geometric_mean <= 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "local VST normalization-factor geometric mean".to_string(),
                index: Some(sample),
                value: column_geometric_mean,
            });
        }
        checked_add_assign(
            &mut inverse_sum,
            column_geometric_mean.recip(),
            sample,
            "local VST inverse normalization-factor mean sum",
        )?;
    }
    Ok(inverse_sum / normalization_factors.n_cols() as f64)
}

/// Apply VST using an already-fitted dispersion trend.
///
/// This mirrors DESeq2's `getVarianceStabilizedData` dispatch once the
/// dispersion function is known. The local branch requires
/// `inverse_size_factor_mean`; parametric and mean branches ignore it.
pub fn vst_with_dispersion_trend(
    normalized_counts: &RowMajorMatrix<f64>,
    trend_fit: &DispersionTrendFit,
    inverse_size_factor_mean: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    match trend_fit {
        DispersionTrendFit::Parametric(fit) => vst_parametric(normalized_counts, fit.trend),
        DispersionTrendFit::Local(fit) => {
            vst_local(normalized_counts, &fit.trend, inverse_size_factor_mean)
        }
        DispersionTrendFit::Mean(fit) => vst_mean(normalized_counts, fit.trend.mean_disp),
    }
}

/// Apply VST using an already-fitted dispersion trend and sample size factors.
///
/// This is the DESeq2-shaped dispatch for callers that have normalized counts
/// and ordinary sample size factors. The local branch computes its
/// `mean(1 / sizeFactors)` variance term internally; parametric and mean
/// branches ignore the size factors after validating them.
pub fn vst_with_dispersion_trend_and_size_factors(
    normalized_counts: &RowMajorMatrix<f64>,
    trend_fit: &DispersionTrendFit,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if size_factors.len() != normalized_counts.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "VST size factors".to_string(),
            expected: normalized_counts.n_cols(),
            actual: size_factors.len(),
        });
    }
    let inverse_size_factor_mean = local_vst_inverse_size_factor_mean(size_factors)?;
    vst_with_dispersion_trend(normalized_counts, trend_fit, inverse_size_factor_mean)
}

/// Apply VST using an already-fitted dispersion trend and normalization factors.
///
/// This mirrors DESeq2's local-VST normalization-factor approximation by
/// deriving column geometric mean size factors and using `mean(1 / sf)`.
/// Parametric and mean branches ignore the derived value.
pub fn vst_with_dispersion_trend_and_normalization_factors(
    normalized_counts: &RowMajorMatrix<f64>,
    trend_fit: &DispersionTrendFit,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if normalization_factors.n_rows() != normalized_counts.n_rows()
        || normalization_factors.n_cols() != normalized_counts.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "VST normalization factors".to_string(),
            expected: normalized_counts.len(),
            actual: normalization_factors.len(),
        });
    }
    let inverse_size_factor_mean =
        local_vst_inverse_size_factor_mean_from_normalization_factors(normalization_factors)?;
    vst_with_dispersion_trend(normalized_counts, trend_fit, inverse_size_factor_mean)
}

/// Apply the currently implemented variance-stabilizing transformation.
///
/// At this stage the convenience alias uses the DESeq2 `fitType="mean"`
/// closed-form branch. Parametric and local trend transforms are exposed as
/// explicit lower-level functions.
pub fn vst(
    normalized_counts: &RowMajorMatrix<f64>,
    mean_dispersion: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    vst_mean(normalized_counts, mean_dispersion)
}

/// Apply DESeq2's mean-fit VST to one normalized count.
pub fn vst_mean_value(
    normalized_count: f64,
    mean_dispersion: f64,
    index: usize,
) -> Result<f64, DeseqError> {
    validate_mean_dispersion(mean_dispersion)?;
    if !normalized_count.is_finite() || normalized_count < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "VST normalized count".to_string(),
            index: Some(index),
            value: normalized_count,
        });
    }
    Ok((2.0 * (mean_dispersion * normalized_count).sqrt().asinh()
        - mean_dispersion.ln()
        - 4.0_f64.ln())
        / std::f64::consts::LN_2)
}

/// Apply DESeq2's parametric-trend VST to one normalized count.
pub fn vst_parametric_value(
    normalized_count: f64,
    trend: ParametricDispersionTrend,
    index: usize,
) -> Result<f64, DeseqError> {
    validate_parametric_trend(trend)?;
    validate_normalized_count(normalized_count, index)?;
    let alpha = trend.asympt_disp;
    let extra = trend.extra_pois;
    let alpha_count = alpha * normalized_count;
    let numerator_root = alpha_count.sqrt() + (1.0 + extra + alpha_count).sqrt();
    Ok((2.0 * numerator_root.ln() - (4.0 * alpha).ln()) / std::f64::consts::LN_2)
}

#[derive(Clone, Debug)]
struct LocalVstIntegral {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl LocalVstIntegral {
    fn fit(
        max_count: f64,
        trend: &LocalDispersionTrend,
        inverse_size_factor_mean: f64,
    ) -> Result<Self, DeseqError> {
        let grid_len = 1000usize;
        let max_asinh = max_count.asinh();
        let mut grid = Vec::with_capacity(grid_len - 1);
        for step in 1..grid_len {
            let value = (max_asinh * step as f64 / (grid_len as f64 - 1.0)).sinh();
            grid.push(value);
        }
        let integrand = grid
            .iter()
            .copied()
            .map(|q| {
                let dispersion = trend.evaluate(q)?;
                let dispersion_q = checked_mul(dispersion, q, 0, "local VST dispersion q")?;
                let dispersion_q2 = checked_mul(dispersion_q, q, 0, "local VST dispersion q2")?;
                let poisson_q = checked_mul(inverse_size_factor_mean, q, 0, "local VST Poisson q")?;
                let variance =
                    checked_add(dispersion_q2, poisson_q, 0, "local VST variance curve")?;
                if !variance.is_finite() || variance <= 0.0 {
                    return Err(DeseqError::InvalidDispersion {
                        reason: "local VST variance curve must be finite and positive".to_string(),
                    });
                }
                Ok(variance.sqrt().recip())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut x = Vec::with_capacity(grid.len() - 1);
        let mut y = Vec::with_capacity(grid.len() - 1);
        let mut cumulative = 0.0;
        for idx in 1..grid.len() {
            let width = grid[idx] - grid[idx - 1];
            let height_sum = checked_add(
                integrand[idx],
                integrand[idx - 1],
                idx,
                "local VST integration height sum",
            )?;
            let area = checked_mul(width, height_sum, idx, "local VST integration area")? / 2.0;
            cumulative = checked_add(cumulative, area, idx, "local VST integration cumulative")?;
            x.push(((grid[idx] + grid[idx - 1]) / 2.0).asinh());
            y.push(cumulative);
        }
        if x.is_empty() {
            return Err(DeseqError::InvalidDimensions {
                context: "local VST integration grid".to_string(),
                expected: 2,
                actual: x.len(),
            });
        }
        Ok(Self { x, y })
    }

    fn evaluate(&self, target: f64) -> Result<f64, DeseqError> {
        if !target.is_finite() || target < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "local VST interpolation target".to_string(),
                index: None,
                value: target,
            });
        }
        if target <= self.x[0] {
            let value = self.y[0] * target / self.x[0];
            return finite_value(value, None, "local VST interpolation lower extrapolation");
        }
        let last = self.x.len() - 1;
        if target >= self.x[last] {
            let slope = (self.y[last] - self.y[last - 1]) / (self.x[last] - self.x[last - 1]);
            let value = self.y[last] + slope * (target - self.x[last]);
            return finite_value(value, None, "local VST interpolation upper extrapolation");
        }
        let upper = self.x.partition_point(|value| *value < target);
        let lower = upper - 1;
        let fraction = (target - self.x[lower]) / (self.x[upper] - self.x[lower]);
        finite_value(
            self.y[lower] + fraction * (self.y[upper] - self.y[lower]),
            None,
            "local VST interpolation",
        )
    }
}

fn validate_mean_dispersion(mean_dispersion: f64) -> Result<(), DeseqError> {
    if !mean_dispersion.is_finite() || mean_dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "mean dispersion for VST must be finite and positive".to_string(),
        });
    }
    Ok(())
}

fn validate_parametric_trend(trend: ParametricDispersionTrend) -> Result<(), DeseqError> {
    if !trend.asympt_disp.is_finite() || trend.asympt_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric VST asymptotic dispersion must be finite and positive".to_string(),
        });
    }
    if !trend.extra_pois.is_finite() || trend.extra_pois < 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric VST extra-Poisson coefficient must be finite and non-negative"
                .to_string(),
        });
    }
    Ok(())
}

fn round_half_to_even(value: f64) -> f64 {
    let floor = value.floor();
    let fraction = value - floor;
    if (fraction - 0.5).abs() < f64::EPSILON {
        if (floor as u64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        value.round()
    }
}

fn validate_normalized_count(normalized_count: f64, index: usize) -> Result<(), DeseqError> {
    if !normalized_count.is_finite() || normalized_count < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "VST normalized count".to_string(),
            index: Some(index),
            value: normalized_count,
        });
    }
    Ok(())
}

fn validate_inverse_size_factor_mean(value: f64) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "local VST inverse size-factor mean must be finite and positive".to_string(),
        });
    }
    Ok(())
}

fn validate_normalized_counts_and_max(
    normalized_counts: &RowMajorMatrix<f64>,
) -> Result<f64, DeseqError> {
    let mut max_count = 0.0;
    for (idx, count) in normalized_counts.as_slice().iter().copied().enumerate() {
        validate_normalized_count(count, idx)?;
        if count > max_count {
            max_count = count;
        }
    }
    Ok(max_count)
}

fn normalized_count_row_means(
    normalized_counts: &RowMajorMatrix<f64>,
) -> Result<Vec<f64>, DeseqError> {
    let mut means = Vec::with_capacity(normalized_counts.n_rows());
    for row in 0..normalized_counts.n_rows() {
        let values = normalized_counts.row(row)?;
        means
            .push(checked_sum(values.iter().copied(), "local VST row mean")? / values.len() as f64);
    }
    Ok(means)
}

fn quantile_type7(values: &[f64], probability: f64) -> Result<f64, DeseqError> {
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "local VST quantile values".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    if sorted.iter().any(|value| !value.is_finite()) {
        return Err(DeseqError::NonFiniteValue {
            context: "local VST row means".to_string(),
            index: None,
            value: f64::NAN,
        });
    }
    if sorted.len() == 1 {
        return Ok(sorted[0]);
    }
    let h = (sorted.len() as f64 - 1.0) * probability + 1.0;
    let lower = h.floor() as usize;
    let fraction = h - lower as f64;
    let lower_idx = lower.saturating_sub(1);
    if lower_idx + 1 >= sorted.len() {
        return Ok(sorted[sorted.len() - 1]);
    }
    Ok(sorted[lower_idx] + fraction * (sorted[lower_idx + 1] - sorted[lower_idx]))
}

fn checked_add(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left + right;
    finite_value(value, Some(index), context)
}

fn checked_add_assign(
    sum: &mut f64,
    term: f64,
    index: usize,
    context: &str,
) -> Result<(), DeseqError> {
    *sum = checked_add(*sum, term, index, context)?;
    Ok(())
}

fn checked_mul(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left * right;
    finite_value(value, Some(index), context)
}

fn checked_sum(values: impl IntoIterator<Item = f64>, context: &str) -> Result<f64, DeseqError> {
    let mut sum = 0.0;
    for (idx, value) in values.into_iter().enumerate() {
        sum = checked_add(sum, value, idx, context)?;
    }
    Ok(sum)
}

fn finite_value(value: f64, index: Option<usize>, context: &str) -> Result<f64, DeseqError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index,
            value,
        })
    }
}
