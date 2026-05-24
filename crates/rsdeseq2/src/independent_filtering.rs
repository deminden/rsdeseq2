use crate::errors::{invalid_dimensions, DeseqError};
use crate::multiple_testing::bh_adjust;
use crate::results::{recompute_padj, DeseqResults};

/// One filtered adjusted-p-value column.
pub type FilteredPadjColumn = Vec<Option<f64>>;

/// Filtered adjusted-p-value columns and rejection counts.
pub type FilteredPAdjustments = (Vec<FilteredPadjColumn>, Vec<usize>);

/// Options for DESeq2-style independent filtering during result assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct IndependentFilteringOptions {
    /// Whether independent filtering is enabled.
    pub enabled: bool,
    /// Target FDR threshold used to count rejections.
    pub alpha: f64,
    /// Optional theta grid. Values are filter quantiles in `[0, 1]`.
    pub theta: Option<Vec<f64>>,
}

impl Default for IndependentFilteringOptions {
    fn default() -> Self {
        Self {
            enabled: true,
            alpha: 0.1,
            theta: None,
        }
    }
}

/// Metadata from independent filtering.
#[derive(Clone, Debug, PartialEq)]
pub struct IndependentFilteringOutput {
    /// Whether independent filtering was applied.
    pub enabled: bool,
    /// Theta grid used to evaluate candidate filter cutoffs.
    pub theta: Vec<f64>,
    /// Number of adjusted p-values below `alpha` at each theta.
    pub num_rejections: Vec<usize>,
    /// Selected theta index.
    pub selected_index: Option<usize>,
    /// Selected theta value.
    pub filter_theta: Option<f64>,
    /// Selected filter threshold on the original filter scale.
    pub filter_threshold: Option<f64>,
    /// Lowess fitted rejection curve used for threshold selection.
    ///
    /// This corresponds to the `lo.fit$y` metadata that DESeq2 stores from
    /// `stats::lowess(numRej ~ theta, f=1/5)`.
    pub lowess_fit: Option<Vec<f64>>,
    /// Alpha used for rejection counting.
    pub alpha: f64,
}

/// One row of the DESeq2-style `filterNumRej` metadata table.
#[derive(Clone, Debug, PartialEq)]
pub struct IndependentFilterNumRejRow {
    /// Filter quantile evaluated by independent filtering.
    pub theta: f64,
    /// Number of adjusted p-values below the target alpha at this theta.
    pub num_rejections: usize,
}

/// One row of the DESeq2-style `lo.fit` independent-filtering metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct IndependentFilterLowessRow {
    /// Filter quantile used as the lowess x coordinate.
    pub theta: f64,
    /// Smoothed rejection count at this theta.
    pub fitted_rejections: f64,
}

/// One scalar entry from DESeq2-style independent-filtering metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct IndependentFilterMetadataEntry {
    /// Metadata key, such as `filterTheta`, `filterThreshold`, or `alpha`.
    pub name: String,
    /// Metadata value as a numeric scalar.
    pub value: f64,
}

impl IndependentFilteringOutput {
    /// Assemble the paired theta/rejection-count metadata table.
    ///
    /// DESeq2 stores this as `metadata(res)$filterNumRej` with `theta` and
    /// `numRej` columns. Rust keeps the raw vectors for direct access and
    /// exposes this paired view for wrappers and parity exporters.
    pub fn filter_num_rej(&self) -> Vec<IndependentFilterNumRejRow> {
        self.theta
            .iter()
            .copied()
            .zip(self.num_rejections.iter().copied())
            .map(|(theta, num_rejections)| IndependentFilterNumRejRow {
                theta,
                num_rejections,
            })
            .collect()
    }

    /// Assemble the paired lowess metadata table.
    ///
    /// DESeq2 stores this as `metadata(res)$lo.fit`, whose x coordinates are
    /// the theta grid and whose y coordinates are the fitted rejection counts.
    pub fn lowess_fit_table(&self) -> Vec<IndependentFilterLowessRow> {
        let Some(lowess_fit) = &self.lowess_fit else {
            return Vec::new();
        };
        self.theta
            .iter()
            .copied()
            .zip(lowess_fit.iter().copied())
            .map(|(theta, fitted_rejections)| IndependentFilterLowessRow {
                theta,
                fitted_rejections,
            })
            .collect()
    }

    /// Assemble scalar independent-filtering metadata entries.
    ///
    /// DESeq2 stores these values in `metadata(res)`. Missing optional entries
    /// are omitted, which is the disabled-filtering shape.
    pub fn scalar_metadata(&self) -> Vec<IndependentFilterMetadataEntry> {
        let mut entries = Vec::new();
        if let Some(value) = self.filter_threshold {
            entries.push(IndependentFilterMetadataEntry {
                name: "filterThreshold".to_string(),
                value,
            });
        }
        if let Some(value) = self.filter_theta {
            entries.push(IndependentFilterMetadataEntry {
                name: "filterTheta".to_string(),
                value,
            });
        }
        entries.push(IndependentFilterMetadataEntry {
            name: "alpha".to_string(),
            value: self.alpha,
        });
        entries
    }
}

/// Apply DESeq2-style independent filtering to result rows.
///
/// The current Rust path uses `baseMean` as the filter statistic, matching the
/// default DESeq2 `results()` behavior. The threshold selection follows the
/// DESeq2 logic around `filtered_p` and rejection counts, with an
/// R-`lowess`-shaped fit for the smoothed rejection curve.
pub fn apply_independent_filtering(
    results: &mut DeseqResults,
    options: &IndependentFilteringOptions,
) -> Result<IndependentFilteringOutput, DeseqError> {
    validate_alpha(options.alpha)?;
    if !options.enabled {
        recompute_padj(results);
        for row in &mut results.rows {
            row.filtered = None;
        }
        let output = IndependentFilteringOutput {
            enabled: false,
            theta: Vec::new(),
            num_rejections: Vec::new(),
            selected_index: None,
            filter_theta: None,
            filter_threshold: None,
            lowess_fit: None,
            alpha: options.alpha,
        };
        results.independent_filtering = Some(output.clone());
        return Ok(output);
    }

    let filter = results
        .rows
        .iter()
        .map(|row| row.base_mean)
        .collect::<Vec<_>>();
    let pvalues = results
        .rows
        .iter()
        .map(|row| row.pvalue)
        .collect::<Vec<_>>();
    let theta = match &options.theta {
        Some(theta) => validate_theta(theta)?,
        None => default_theta(&filter)?,
    };
    let cutoffs = theta
        .iter()
        .copied()
        .map(|value| quantile_type7(&filter, value))
        .collect::<Result<Vec<_>, _>>()?;
    let (columns, num_rejections) =
        filtered_p_adjustments(&filter, &pvalues, &cutoffs, options.alpha)?;
    let (selected_index, lowess_fit) = select_filter_index_with_lowess(&theta, &num_rejections)?;
    let selected_padj = columns.get(selected_index).ok_or_else(|| {
        invalid_dimensions(
            "independent-filter columns",
            selected_index + 1,
            columns.len(),
        )
    })?;
    let filter_threshold = cutoffs[selected_index];

    for (row, adjusted) in results.rows.iter_mut().zip(selected_padj.iter().copied()) {
        row.padj = adjusted;
        row.filtered = row.pvalue.map(|_| row.base_mean < filter_threshold);
    }

    let output = IndependentFilteringOutput {
        enabled: true,
        theta: theta.clone(),
        num_rejections,
        selected_index: Some(selected_index),
        filter_theta: Some(theta[selected_index]),
        filter_threshold: Some(filter_threshold),
        lowess_fit: Some(lowess_fit),
        alpha: options.alpha,
    };
    results.independent_filtering = Some(output.clone());
    Ok(output)
}

/// Construct the DESeq2 default theta grid from filter values.
pub fn default_theta(filter: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if filter.is_empty() {
        return Err(invalid_dimensions("independent-filter values", 1, 0));
    }
    validate_filter(filter)?;
    let zero_count = filter.iter().filter(|value| **value == 0.0).count();
    let lower = zero_count as f64 / filter.len() as f64;
    let upper = if lower < 0.95 { 0.95 } else { 1.0 };
    Ok(seq(lower, upper, 50))
}

/// Compute filtered BH-adjusted p-values for each candidate cutoff.
pub fn filtered_p_adjustments(
    filter: &[f64],
    pvalues: &[Option<f64>],
    cutoffs: &[f64],
    alpha: f64,
) -> Result<FilteredPAdjustments, DeseqError> {
    validate_alpha(alpha)?;
    validate_filter(filter)?;
    if filter.len() != pvalues.len() {
        return Err(invalid_dimensions(
            "independent-filter p-values",
            filter.len(),
            pvalues.len(),
        ));
    }

    let mut columns = Vec::with_capacity(cutoffs.len());
    let mut num_rejections = Vec::with_capacity(cutoffs.len());
    for cutoff in cutoffs {
        if !cutoff.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "independent-filter cutoff".to_string(),
                index: None,
                value: *cutoff,
            });
        }
        let mut selected = Vec::new();
        let mut selected_indices = Vec::new();
        for (idx, (filter_value, pvalue)) in filter
            .iter()
            .copied()
            .zip(pvalues.iter().copied())
            .enumerate()
        {
            if filter_value >= *cutoff {
                selected.push(pvalue);
                selected_indices.push(idx);
            }
        }
        let selected_padj = bh_adjust(&selected);
        let mut full = vec![None; filter.len()];
        for (idx, adjusted) in selected_indices.into_iter().zip(selected_padj) {
            full[idx] = adjusted;
        }
        let rejections = full
            .iter()
            .filter(|value| value.is_some_and(|adjusted| adjusted < alpha))
            .count();
        columns.push(full);
        num_rejections.push(rejections);
    }
    Ok((columns, num_rejections))
}

/// Select the independent-filter threshold index from rejection counts.
pub fn select_filter_index(theta: &[f64], num_rejections: &[usize]) -> usize {
    select_filter_index_with_lowess(theta, num_rejections)
        .map(|(idx, _)| idx)
        .unwrap_or(0)
}

/// Select the independent-filter threshold index and return the lowess fit.
pub fn select_filter_index_with_lowess(
    theta: &[f64],
    num_rejections: &[usize],
) -> Result<(usize, Vec<f64>), DeseqError> {
    if theta.len() != num_rejections.len() {
        return Err(invalid_dimensions(
            "independent-filter theta/rejection lengths",
            theta.len(),
            num_rejections.len(),
        ));
    }
    let smooth = lowess_fitted_values(theta, num_rejections, 0.2)?;
    if num_rejections.is_empty() || num_rejections.iter().copied().max().unwrap_or(0) <= 10 {
        return Ok((0, smooth));
    }
    let positive_residuals = num_rejections
        .iter()
        .copied()
        .zip(smooth.iter().copied())
        .filter_map(|(count, fitted)| (count > 0).then_some(count as f64 - fitted))
        .collect::<Vec<_>>();
    let rmse = if positive_residuals.is_empty() {
        0.0
    } else {
        (positive_residuals
            .iter()
            .map(|value| value.powi(2))
            .sum::<f64>()
            / positive_residuals.len() as f64)
            .sqrt()
    };
    let max_fit = smooth.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let threshold = max_fit - rmse;

    let selected = first_index_above(num_rejections, threshold)
        .or_else(|| first_index_above(num_rejections, 0.9 * max_fit))
        .or_else(|| first_index_above(num_rejections, 0.8 * max_fit))
        .unwrap_or(0);
    Ok((selected, smooth))
}

/// R-`lowess`-shaped fitted values for independent-filter rejection counts.
///
/// DESeq2 calls `stats::lowess(numRej ~ theta, f=1/5)` when selecting the
/// independent-filter threshold. This helper follows the same algorithmic
/// choices used by R for the DESeq2 threshold grid: floor-based nearest-neighbor
/// span, tricube distance weights, local linear fitting, three Tukey biweight
/// robustifying iterations, and the default `delta = 0.01 * diff(range(theta))`
/// interpolation shortcut.
pub fn lowess_fitted_values(
    theta: &[f64],
    num_rejections: &[usize],
    span_fraction: f64,
) -> Result<Vec<f64>, DeseqError> {
    if theta.len() != num_rejections.len() {
        return Err(invalid_dimensions(
            "lowess theta/rejection lengths",
            theta.len(),
            num_rejections.len(),
        ));
    }
    if !span_fraction.is_finite() || span_fraction <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "lowess span fraction".to_string(),
            index: None,
            value: span_fraction,
        });
    }
    let n = theta.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![num_rejections[0] as f64]);
    }
    validate_theta(theta)?;

    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|left, right| {
        theta[*left]
            .total_cmp(&theta[*right])
            .then_with(|| left.cmp(right))
    });
    let sorted_x = order.iter().map(|idx| theta[*idx]).collect::<Vec<_>>();
    let sorted_y = order
        .iter()
        .map(|idx| num_rejections[*idx] as f64)
        .collect::<Vec<_>>();

    let delta = 0.01 * (sorted_x[n - 1] - sorted_x[0]);
    let fitted_sorted = lowess_sorted(&sorted_x, &sorted_y, span_fraction, 3, delta);
    let mut fitted = vec![0.0; n];
    for (sorted_idx, original_idx) in order.into_iter().enumerate() {
        fitted[original_idx] = fitted_sorted[sorted_idx];
    }
    Ok(fitted)
}

fn validate_alpha(alpha: f64) -> Result<(), DeseqError> {
    if alpha > 0.0 && alpha < 1.0 {
        Ok(())
    } else {
        Err(DeseqError::NonFiniteValue {
            context: "independent-filter alpha".to_string(),
            index: None,
            value: alpha,
        })
    }
}

fn validate_theta(theta: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if theta.len() <= 1 {
        return Err(invalid_dimensions(
            "independent-filter theta",
            2,
            theta.len(),
        ));
    }
    let mut values = Vec::with_capacity(theta.len());
    for (idx, value) in theta.iter().copied().enumerate() {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(DeseqError::NonFiniteValue {
                context: "independent-filter theta".to_string(),
                index: Some(idx),
                value,
            });
        }
        values.push(value);
    }
    Ok(values)
}

fn validate_filter(filter: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in filter.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "independent-filter values".to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

fn quantile_type7(values: &[f64], probability: f64) -> Result<f64, DeseqError> {
    validate_filter(values)?;
    if values.is_empty() {
        return Err(invalid_dimensions("quantile values", 1, 0));
    }
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::NonFiniteValue {
            context: "quantile probability".to_string(),
            index: None,
            value: probability,
        });
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    if sorted.len() == 1 {
        return Ok(sorted[0]);
    }
    let h = 1.0 + (sorted.len() as f64 - 1.0) * probability;
    let floor = h.floor();
    let lower_idx = floor as usize - 1;
    let gamma = h - floor;
    if lower_idx + 1 >= sorted.len() {
        Ok(sorted[sorted.len() - 1])
    } else {
        Ok(sorted[lower_idx] + gamma * (sorted[lower_idx + 1] - sorted[lower_idx]))
    }
}

fn seq(start: f64, end: f64, len: usize) -> Vec<f64> {
    if len == 1 {
        return vec![start];
    }
    let step = (end - start) / (len - 1) as f64;
    (0..len).map(|idx| start + step * idx as f64).collect()
}

fn lowess_sorted(
    x: &[f64],
    y: &[f64],
    span_fraction: f64,
    robust_iterations: usize,
    delta: f64,
) -> Vec<f64> {
    let n = x.len();
    let span = ((span_fraction * n as f64 + 1e-7) as usize).clamp(2, n);
    let mut robustness_weights = vec![1.0; n];
    let mut fitted = vec![0.0; n];

    for iteration in 0..=robust_iterations {
        fitted = lowess_delta_pass(x, y, span, &robustness_weights, delta);
        if iteration == robust_iterations {
            break;
        }
        let residuals = y
            .iter()
            .copied()
            .zip(fitted.iter().copied())
            .map(|(observed, fitted)| (observed - fitted).abs())
            .collect::<Vec<_>>();
        let mean_absolute_residual = residuals.iter().sum::<f64>() / residuals.len() as f64;
        let scale = six_mad(&residuals);
        if scale < 1e-7 * mean_absolute_residual {
            break;
        }
        for (weight, residual) in robustness_weights.iter_mut().zip(residuals.iter().copied()) {
            let absolute = residual.abs();
            if absolute <= 0.001 * scale {
                *weight = 1.0;
            } else if absolute >= 0.999 * scale {
                *weight = 0.0;
            } else {
                let ratio = absolute / scale;
                *weight = (1.0 - ratio.powi(2)).powi(2);
            }
        }
    }
    fitted
}

fn lowess_delta_pass(
    x: &[f64],
    y: &[f64],
    span: usize,
    robustness_weights: &[f64],
    delta: f64,
) -> Vec<f64> {
    let n = x.len();
    let mut fitted = vec![0.0; n];
    let mut n_left = 0_usize;
    let mut n_right = span - 1;
    let mut last_estimated: Option<usize> = None;
    let mut idx = 0_usize;
    loop {
        while n_right < n - 1 {
            let left_radius = x[idx] - x[n_left];
            let next_right_radius = x[n_right + 1] - x[idx];
            if left_radius > next_right_radius {
                n_left += 1;
                n_right += 1;
            } else {
                break;
            }
        }

        let fitted_value =
            lowest_at(x, y, idx, n_left, n_right, robustness_weights).unwrap_or_else(|| y[idx]);
        fitted[idx] = fitted_value;

        if let Some(last) = last_estimated {
            if last < idx - 1 {
                let denominator = x[idx] - x[last];
                for point in last + 1..idx {
                    fitted[point] = if denominator.abs() <= f64::EPSILON {
                        fitted[last]
                    } else {
                        let fraction = (x[point] - x[last]) / denominator;
                        fraction * fitted[idx] + (1.0 - fraction) * fitted[last]
                    };
                }
            }
        }
        let mut last = idx;
        let cutoff = x[last] + delta.max(0.0);
        let mut next = last + 1;
        while next < n {
            if x[next] > cutoff {
                break;
            }
            if x[next] == x[last] {
                fitted[next] = fitted[last];
                last = next;
            }
            next += 1;
        }
        last_estimated = Some(last);
        if last >= n - 1 {
            break;
        }
        idx = (next.saturating_sub(1)).max(last + 1);
    }
    fitted
}

fn lowest_at(
    x: &[f64],
    y: &[f64],
    idx: usize,
    n_left: usize,
    n_right: usize,
    robustness_weights: &[f64],
) -> Option<f64> {
    let target_x = x[idx];
    let range = x[x.len() - 1] - x[0];
    let bandwidth = (target_x - x[n_left]).max(x[n_right] - target_x);
    let high = 0.999 * bandwidth;
    let low = 0.001 * bandwidth;
    let mut weights = Vec::new();
    let mut sum_weights = 0.0;
    let mut point = n_left;
    while point < x.len() {
        let distance = (x[point] - target_x).abs();
        if distance <= high {
            let proximity = if distance <= low {
                1.0
            } else {
                (1.0 - (distance / bandwidth).powi(3)).powi(3)
            };
            let weight = proximity * robustness_weights[point];
            weights.push((point, weight));
            sum_weights += weight;
        } else if x[point] > target_x {
            break;
        }
        point += 1;
    }
    if sum_weights <= 0.0 {
        return None;
    }

    for (_, weight) in &mut weights {
        *weight /= sum_weights;
    }

    if bandwidth > 0.0 {
        let center = weights
            .iter()
            .map(|(point, weight)| weight * x[*point])
            .sum::<f64>();
        let variance = weights
            .iter()
            .map(|(point, weight)| weight * (x[*point] - center).powi(2))
            .sum::<f64>();
        if variance.sqrt() > 0.001 * range {
            let slope_factor = (target_x - center) / variance;
            for (point, weight) in &mut weights {
                *weight *= slope_factor * (x[*point] - center) + 1.0;
            }
        }
    }

    Some(
        weights
            .iter()
            .map(|(point, weight)| weight * y[*point])
            .sum(),
    )
}

fn six_mad(residuals: &[f64]) -> f64 {
    if residuals.is_empty() {
        return 0.0;
    }
    let mut sorted = residuals.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = if sorted.len() % 2 == 0 {
        let upper = sorted.len() / 2;
        0.5 * (sorted[upper - 1] + sorted[upper])
    } else {
        sorted[sorted.len() / 2]
    };
    6.0 * median
}

fn first_index_above(values: &[usize], threshold: f64) -> Option<usize> {
    values
        .iter()
        .copied()
        .enumerate()
        .find_map(|(idx, value)| (value as f64 > threshold).then_some(idx))
}
