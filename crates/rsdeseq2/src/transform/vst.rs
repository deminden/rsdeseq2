use crate::core::CountMatrix;
use crate::dispersion::{DispersionTrendFit, LocalDispersionTrend, ParametricDispersionTrend};
use crate::errors::DeseqError;
use crate::matrix::RowMajorMatrix;
use crate::numeric::checked::{self, CheckedExt};

// Fast-VST row selection and aligned subsets.
mod frozen;
mod subset;
pub use frozen::{
    FrozenLocalDispersion, SavedLocalDispersionFit, frozen_local_vst_grid_from_counts,
    vst_frozen_local,
};
pub use subset::{
    FastVstSubset, FastVstSubsetMetadata, fast_vst_eligible_count, fast_vst_subset,
    fast_vst_subset_indices, fast_vst_subset_matrix_rows, fast_vst_subset_normalized_counts,
};

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
    vst_local_with_integral(normalized_counts, &integral, h1, h2)
}

/// Apply local VST using dispersion-function values sampled on the DESeq2 grid.
///
/// This supports frozen local dispersion functions whose fitted representation
/// is owned by the caller (for example, an R closure). The input grid must be
/// the increasing 999-point grid used by DESeq2, without the initial zero.
pub fn vst_local_from_dispersion_grid(
    normalized_counts: &RowMajorMatrix<f64>,
    grid: &[f64],
    dispersions: &[f64],
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
    let integral = LocalVstIntegral::fit_from_grid(grid, dispersions, inverse_size_factor_mean)?;
    vst_local_with_integral(normalized_counts, &integral, h1, h2)
}

fn vst_local_with_integral(
    normalized_counts: &RowMajorMatrix<f64>,
    integral: &LocalVstIntegral,
    h1: f64,
    h2: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let int_h1 = integral.evaluate(h1.asinh())?;
    let int_h2 = integral.evaluate(h2.asinh())?;
    if int_h2 <= int_h1 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local VST integral must increase across scaling quantiles".to_string(),
        });
    }
    let log_delta = checked::try_sub(h2.log2(), h1.log2())
        .nonfinite_at("local VST scaling log delta", Some(0))?;
    let integral_delta = checked::try_sub(int_h2, int_h1)
        .nonfinite_at("local VST scaling integral delta", Some(0))?;
    let eta = checked::try_div_finite_result(log_delta, integral_delta)
        .nonfinite_at("local VST scaling eta", Some(0))?;
    let eta_integral =
        checked::try_mul(eta, int_h1).nonfinite_at("local VST scaling eta integral", Some(0))?;
    let xi =
        checked::try_sub(h1.log2(), eta_integral).nonfinite_at("local VST scaling xi", Some(0))?;
    let values = normalized_counts
        .as_slice()
        .iter()
        .copied()
        .map(|count| {
            integral.evaluate(count.asinh()).and_then(|value| {
                checked::try_add(
                    checked::try_mul(eta, value).nonfinite_at("local VST scaled value", Some(0))?,
                    xi,
                )
                .nonfinite_at("local VST scaled value", Some(0))
            })
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
    let mut inverse_size_factors = Vec::with_capacity(size_factors.len());
    for (idx, size_factor) in size_factors.iter().copied().enumerate() {
        if !size_factor.is_finite() || size_factor <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("local VST size factor at index {idx} must be finite and positive"),
            });
        }
        inverse_size_factors.push(size_factor.recip());
    }
    checked::try_mean(&inverse_size_factors).nonfinite("local VST inverse size-factor mean")
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
    let mut inverse_column_means = Vec::with_capacity(normalization_factors.n_cols());
    for (sample, sum) in log_col_sums.into_iter().enumerate() {
        let log_mean = checked::try_div_finite_result(sum, normalization_factors.n_rows() as f64)
            .nonfinite_at(
            "local VST normalization-factor log column mean",
            Some(sample),
        )?;
        let column_geometric_mean = log_mean.exp();
        if !column_geometric_mean.is_finite() || column_geometric_mean <= 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "local VST normalization-factor geometric mean".to_string(),
                index: Some(sample),
                value: column_geometric_mean,
            });
        }
        inverse_column_means.push(column_geometric_mean.recip());
    }
    checked::try_mean(&inverse_column_means)
        .nonfinite("local VST inverse normalization-factor mean")
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
    let dispersion_count = match checked::try_mul(mean_dispersion, normalized_count)
        .nonfinite_at("mean VST dispersion count", Some(index))
    {
        Ok(value) => value,
        Err(_) => return Ok(normalized_count.log2()),
    };
    if dispersion_count.is_infinite() {
        return Ok(normalized_count.log2());
    }
    let transformed_numerator = checked::try_sub(
        checked::try_mul(2.0, dispersion_count.sqrt().asinh())
            .nonfinite_at("mean VST numerator", Some(index))?,
        mean_dispersion.ln(),
    )
    .nonfinite_at("mean VST numerator", Some(index))?;
    let transformed_numerator = checked::try_sub(transformed_numerator, 4.0_f64.ln())
        .nonfinite_at("mean VST numerator", Some(index))?;
    finite_value(
        checked::try_div_finite_result(transformed_numerator, std::f64::consts::LN_2)
            .nonfinite_at("mean VST value", Some(index))?,
        Some(index),
        "mean VST value",
    )
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
    let alpha_count = match checked::try_mul(alpha, normalized_count)
        .nonfinite_at("parametric VST alpha count", Some(index))
    {
        Ok(value) => value,
        Err(_) => return Ok(normalized_count.log2()),
    };
    if alpha_count.is_infinite() {
        return Ok(normalized_count.log2());
    }
    let root_term = checked::try_add(1.0, extra)
        .nonfinite_at("parametric VST root term", Some(index))
        .and_then(|value| {
            checked::try_add(value, alpha_count)
                .nonfinite_at("parametric VST root term", Some(index))
        })?;
    let numerator_root = checked::try_add(alpha_count.sqrt(), root_term.sqrt())
        .nonfinite_at("parametric VST numerator root", Some(index))?;
    let numerator = checked::try_sub(
        checked::try_mul(2.0, numerator_root.ln())
            .nonfinite_at("parametric VST numerator", Some(index))?,
        checked::try_mul(4.0, alpha)
            .nonfinite_at("parametric VST denominator", Some(index))?
            .ln(),
    )
    .nonfinite_at("parametric VST numerator", Some(index))?;
    finite_value(
        checked::try_div_finite_result(numerator, std::f64::consts::LN_2)
            .nonfinite_at("parametric VST value", Some(index))?,
        Some(index),
        "parametric VST value",
    )
}

// Local VST integration and FMM spline interpolation.
mod local_integral;
use local_integral::LocalVstIntegral;

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
        if (floor as u64).is_multiple_of(2) {
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
        means.push(checked::try_mean(values).nonfinite("local VST row mean")?);
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
    checked_interpolate(
        sorted[lower_idx],
        fraction,
        checked::try_sub(sorted[lower_idx + 1], sorted[lower_idx])
            .nonfinite_at("local VST quantile interpolation delta", Some(lower_idx))?,
        "local VST quantile interpolation",
    )
}

fn checked_add_assign(
    sum: &mut f64,
    term: f64,
    index: usize,
    context: &str,
) -> Result<(), DeseqError> {
    *sum = checked::try_add(*sum, term).nonfinite_at(context, Some(index))?;
    Ok(())
}

fn checked_interpolate(
    origin: f64,
    scale: f64,
    delta: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    checked::try_add(
        origin,
        checked::try_mul(scale, delta).nonfinite_at(context, Some(0))?,
    )
    .nonfinite_at(context, Some(0))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_interpolate_rejects_overflowed_product() {
        let err = checked_interpolate(0.0, f64::MAX, 2.0, "test interpolation").unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test interpolation" && index == Some(0)
        ));
    }

    #[test]
    fn checked_local_vst_arithmetic_rejects_nonfinite_subtraction_and_division() {
        let sub_err = checked::try_sub(f64::MAX, -f64::MAX)
            .nonfinite_at("test subtraction", Some(0))
            .unwrap_err();
        assert!(matches!(
            sub_err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test subtraction" && index == Some(0)
        ));

        let div_err = checked::try_div_finite_result(1.0, 0.0)
            .nonfinite_at("test division", Some(1))
            .unwrap_err();
        assert!(matches!(
            div_err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test division" && index == Some(1)
        ));
    }

    #[test]
    fn quantile_type7_keeps_large_equal_values_finite() {
        let value = quantile_type7(&[f64::MAX / 2.0, f64::MAX / 2.0], 0.5).unwrap();

        assert_eq!(value, f64::MAX / 2.0);
    }

    #[test]
    fn quantile_type7_rejects_overflowed_delta() {
        let err = quantile_type7(&[-f64::MAX, f64::MAX], 0.5).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "local VST quantile interpolation delta" && index == Some(0)
        ));
    }
}
