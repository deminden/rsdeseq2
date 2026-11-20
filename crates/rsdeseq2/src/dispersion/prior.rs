use crate::errors::{DeseqError, invalid_dimensions};
use crate::errors::{ErrorKind, Requirement, invalid_calculation, invalid_parameter};
use crate::math::{median_finite, trigamma};
use crate::numeric::checked::{self, CheckedExt};
use crate::numeric::math::ProfileMath;
use crate::numeric::special::{GammaDistribution, standard_normal_inverse_cdf};

const LOW_DF_HIST_MIN: f64 = -10.0;
const LOW_DF_HIST_MAX: f64 = 10.0;
const LOW_DF_HIST_WIDTH: f64 = 0.5;
const LOW_DF_VAR_GRID_POINTS: usize = 200;
const LOW_DF_FINE_GRID_POINTS: usize = 1000;
const LOW_DF_VAR_MAX: f64 = 8.0;
const LOW_DF_QUASI_SAMPLES: usize = 10_000;
const LOW_DF_LOESS_SPAN: f64 = 0.2;

/// Output from DESeq2-style dispersion prior variance estimation.
#[derive(Clone, Debug, PartialEq)]
pub struct DispersionPriorVarianceOutput {
    /// Final prior variance for the normal prior on `log(alpha)`.
    pub disp_prior_var: f64,
    /// Robust variance estimate of log-dispersion residuals.
    pub var_log_disp_estimates: f64,
    /// Expected sampling variance subtracted from `var_log_disp_estimates`.
    pub expected_log_dispersion_variance: f64,
    /// Residual degrees of freedom, `n_samples - n_coefficients`.
    pub residual_degrees_of_freedom: usize,
    /// Rows used for the residual variance estimate.
    pub above_min_disp: Vec<bool>,
}

/// Estimate DESeq2's dispersion prior variance.
///
/// DESeq2 computes `mad(log(dispGeneEst) - log(dispFit))^2` over rows where
/// `dispGeneEst >= 100 * minDisp`. If residual degrees of freedom are greater
/// than three, it subtracts `trigamma((m - p) / 2)` and floors the result at
/// `0.25`. If `m == p`, no sampling-variance subtraction is performed.
///
/// For residual degrees of freedom 1 through 3, DESeq2 uses a seeded
/// Monte-Carlo histogram/KL match followed by loess smoothing. This
/// implementation mirrors that shape with deterministic quasi-random samples
/// and local-linear smoothing, so repeated Rust runs are bit-stable while
/// preserving the intended low-df behavior.
pub fn estimate_dispersion_prior_variance(
    disp_gene_est: &[f64],
    disp_fit: &[f64],
    min_disp: f64,
    n_samples: usize,
    n_coefficients: usize,
) -> Result<DispersionPriorVarianceOutput, DeseqError> {
    validate_prior_variance_inputs(disp_gene_est, disp_fit, min_disp, n_samples, n_coefficients)?;
    let residual_degrees_of_freedom = n_samples - n_coefficients;
    let (residuals, above_min_disp) =
        log_dispersion_residuals_above_min(disp_gene_est, disp_fit, min_disp)?;
    if residuals.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "no data found which is greater than minDisp".to_string(),
        });
    }

    let var_log_disp_estimates = mad_squared(&residuals)?;
    let expected_log_dispersion_variance = if residual_degrees_of_freedom == 0 {
        0.0
    } else {
        trigamma(residual_degrees_of_freedom as f64 / 2.0)?
    };
    let disp_prior_var = if (1..=3).contains(&residual_degrees_of_freedom) {
        estimate_low_df_prior_variance(&residuals, residual_degrees_of_freedom)?
    } else if residual_degrees_of_freedom > 0 {
        (var_log_disp_estimates - expected_log_dispersion_variance).max(0.25)
    } else {
        var_log_disp_estimates
    };

    Ok(DispersionPriorVarianceOutput {
        disp_prior_var,
        var_log_disp_estimates,
        expected_log_dispersion_variance,
        residual_degrees_of_freedom,
        above_min_disp,
    })
}

/// Compute DESeq2's robust `mad(x)^2` estimate with R's default constant.
pub fn mad_squared(values: &[f64]) -> Result<f64, DeseqError> {
    let center = median_finite(values).ok_or_else(|| DeseqError::InvalidDispersion {
        reason: "cannot compute MAD of an empty or non-finite slice".to_string(),
    })?;
    let deviations = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .map(|value| (value - center).abs())
        .collect::<Vec<_>>();
    let mad = median_finite(&deviations).ok_or_else(|| DeseqError::InvalidDispersion {
        reason: "cannot compute MAD of an empty or non-finite slice".to_string(),
    })? * 1.4826;
    let variance = mad * mad;
    if variance.is_finite() {
        Ok(variance)
    } else {
        Err(invalid_calculation(
            ErrorKind::Dispersion,
            "MAD squared",
            None,
            variance,
            "MAD squared produced non-finite dispersion prior variance",
        ))
    }
}

/// Log-dispersion residuals for rows used by DESeq2's prior variance estimate.
pub fn log_dispersion_residuals_above_min(
    disp_gene_est: &[f64],
    disp_fit: &[f64],
    min_disp: f64,
) -> Result<(Vec<f64>, Vec<bool>), DeseqError> {
    if disp_gene_est.len() != disp_fit.len() {
        return Err(invalid_dimensions(
            "dispersion prior variance rows",
            disp_gene_est.len(),
            disp_fit.len(),
        ));
    }
    if !min_disp.is_finite() || min_disp <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "min_disp",
            Requirement::FiniteAndPositive,
        ));
    }

    let mut residuals = Vec::new();
    let mut above_min_disp = Vec::with_capacity(disp_gene_est.len());
    for (gene, (gene_est, fit)) in disp_gene_est
        .iter()
        .copied()
        .zip(disp_fit.iter().copied())
        .enumerate()
    {
        if gene_est.is_finite() && gene_est > 0.0 && fit.is_finite() && fit > 0.0 {
            let above_min = gene_est >= min_disp * 100.0;
            above_min_disp.push(above_min);
            if above_min {
                residuals.push((gene_est / fit).ln_p());
            }
        } else if gene_est.is_nan() || fit.is_nan() {
            above_min_disp.push(false);
        } else {
            let (operation, value) = if !gene_est.is_finite() || gene_est <= 0.0 {
                ("gene-wise dispersion", gene_est)
            } else {
                ("fitted dispersion", fit)
            };
            return Err(DeseqError::InvalidGeneValue {
                domain: ErrorKind::Dispersion,
                operation,
                gene,
                sample: None,
                value,
                reason: format!(
                    "dispersion prior variance row {gene} must contain positive finite dispersions"
                ),
            });
        }
    }
    Ok((residuals, above_min_disp))
}

/// Deterministic analogue of DESeq2's low-residual-df histogram/KL branch.
pub fn estimate_low_df_prior_variance(
    residuals: &[f64],
    residual_degrees_of_freedom: usize,
) -> Result<f64, DeseqError> {
    if !(1..=3).contains(&residual_degrees_of_freedom) {
        return Err(DeseqError::InvalidDimensions {
            context: "low-df dispersion prior variance residual df".to_string(),
            expected: 3,
            actual: residual_degrees_of_freedom,
        });
    }
    let observed = residuals
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > LOW_DF_HIST_MIN && *value < LOW_DF_HIST_MAX)
        .collect::<Vec<_>>();
    if observed.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "no finite low-df dispersion residuals fall inside DESeq2 histogram bounds"
                .to_string(),
        });
    }

    let obs_density = histogram_density(&observed)?;
    let (base_samples, normal_samples) =
        low_df_quasi_samples(residual_degrees_of_freedom, LOW_DF_QUASI_SAMPLES)?;
    let variance_grid = linspace(0.0, LOW_DF_VAR_MAX, LOW_DF_VAR_GRID_POINTS)?;
    let kl_divergences = variance_grid
        .iter()
        .copied()
        .map(|variance| {
            let simulated = base_samples
                .iter()
                .copied()
                .zip(normal_samples.iter().copied())
                .map(|(base, normal)| low_df_simulated_residual(base, variance, normal))
                .collect::<Result<Vec<_>, DeseqError>>()?
                .into_iter()
                .filter(|value| *value > LOW_DF_HIST_MIN && *value < LOW_DF_HIST_MAX)
                .collect::<Vec<_>>();
            let simulated_density = histogram_density(&simulated)?;
            kl_divergence(&obs_density, &simulated_density)
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;
    let fine_grid = linspace(0.0, LOW_DF_VAR_MAX, LOW_DF_FINE_GRID_POINTS)?;
    let argmin_kl = local_linear_smoothed_argmin(&variance_grid, &kl_divergences, &fine_grid)?;
    Ok(argmin_kl.max(0.25))
}

/// Estimate DESeq2's dispersion prior variance from gene-wise and trend dispersions.
///
/// This convenience wrapper returns the same detailed output as
/// [`estimate_dispersion_prior_variance`]. It exists for callers that want a
/// stage-shaped entry point named after DESeq2's dispersion-prior step.
pub fn estimate_dispersion_prior(
    disp_gene_est: &[f64],
    disp_fit: &[f64],
    min_disp: f64,
    n_samples: usize,
    n_coefficients: usize,
) -> Result<DispersionPriorVarianceOutput, DeseqError> {
    estimate_dispersion_prior_variance(disp_gene_est, disp_fit, min_disp, n_samples, n_coefficients)
}

fn validate_prior_variance_inputs(
    disp_gene_est: &[f64],
    disp_fit: &[f64],
    min_disp: f64,
    n_samples: usize,
    n_coefficients: usize,
) -> Result<(), DeseqError> {
    if disp_gene_est.len() != disp_fit.len() {
        return Err(invalid_dimensions(
            "dispersion prior variance rows",
            disp_gene_est.len(),
            disp_fit.len(),
        ));
    }
    if disp_gene_est.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion prior variance rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !min_disp.is_finite() || min_disp <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "min_disp",
            Requirement::FiniteAndPositive,
        ));
    }
    if n_samples == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion prior variance samples".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if n_coefficients == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion prior variance coefficients".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if n_samples < n_coefficients {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion prior variance residual degrees of freedom".to_string(),
            expected: n_coefficients,
            actual: n_samples,
        });
    }
    Ok(())
}

fn histogram_density(values: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if values.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "cannot compute low-df dispersion histogram for empty values".to_string(),
        });
    }
    let bins = ((LOW_DF_HIST_MAX - LOW_DF_HIST_MIN) / LOW_DF_HIST_WIDTH).round() as usize;
    let mut counts = vec![0_usize; bins];
    for value in values.iter().copied() {
        if !value.is_finite() || value <= LOW_DF_HIST_MIN || value >= LOW_DF_HIST_MAX {
            continue;
        }
        let mut bin = ((value - LOW_DF_HIST_MIN) / LOW_DF_HIST_WIDTH).floor() as usize;
        if bin >= bins {
            bin = bins - 1;
        }
        counts[bin] += 1;
    }
    let total = counts.iter().sum::<usize>();
    if total == 0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "all low-df dispersion histogram values were outside bounds".to_string(),
        });
    }
    let density_denominator = checked::try_mul(total as f64, LOW_DF_HIST_WIDTH)
        .nonfinite_at("low-df dispersion histogram density denominator", None)?;
    counts
        .into_iter()
        .map(|count| {
            checked::try_div(count as f64, density_denominator)
                .nonfinite_at("low-df dispersion histogram density", None)
        })
        .collect()
}

fn kl_divergence(observed: &[f64], simulated: &[f64]) -> Result<f64, DeseqError> {
    if observed.len() != simulated.len() {
        return Err(invalid_dimensions(
            "low-df dispersion histogram bins",
            observed.len(),
            simulated.len(),
        ));
    }
    let small = observed
        .iter()
        .copied()
        .chain(simulated.iter().copied())
        .filter(|value| *value > 0.0)
        .fold(f64::INFINITY, f64::min);
    if !small.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "low-df dispersion histograms have no positive bins".to_string(),
        });
    }
    let mut divergence = 0.0;
    for (obs, sim) in observed.iter().copied().zip(simulated.iter().copied()) {
        let obs_adjusted =
            checked::try_add(obs, small).nonfinite_at("low-df dispersion KL observed bin", None)?;
        let sim_adjusted = checked::try_add(sim, small)
            .nonfinite_at("low-df dispersion KL simulated bin", None)?;
        let ratio = checked::try_div(obs_adjusted, sim_adjusted)
            .nonfinite_at("low-df dispersion KL bin ratio", None)?;
        let log_ratio = ratio.ln_p();
        if !log_ratio.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "low-df dispersion KL log ratio".to_string(),
                index: None,
                value: log_ratio,
            });
        }
        let term =
            checked::try_mul(obs, log_ratio).nonfinite_at("low-df dispersion KL term", None)?;
        divergence =
            checked::try_add(divergence, term).nonfinite_at("low-df dispersion KL sum", None)?;
    }
    Ok(divergence)
}

fn low_df_quasi_samples(
    residual_degrees_of_freedom: usize,
    n: usize,
) -> Result<(Vec<f64>, Vec<f64>), DeseqError> {
    let chi_squared = GammaDistribution::chi_squared(residual_degrees_of_freedom as f64)
        .ok_or_else(|| DeseqError::InvalidDispersion {
            reason: format!(
                "failed to construct low-df chi-squared distribution: {residual_degrees_of_freedom} degrees of freedom"
            ),
        })?;
    let log_df = (residual_degrees_of_freedom as f64).ln_p();
    let mut base_samples = Vec::with_capacity(n);
    let mut normal_samples = Vec::with_capacity(n);
    for idx in 1..=n {
        let u_chi = clamp_probability(halton(idx, 2));
        let u_norm = clamp_probability(halton(idx, 3));
        base_samples.push(chi_squared.inverse_cdf(u_chi).ln_p() - log_df);
        normal_samples.push(standard_normal_inverse_cdf(u_norm));
    }
    Ok((base_samples, normal_samples))
}

fn halton(mut index: usize, base: usize) -> f64 {
    let mut fraction = 1.0;
    let mut result = 0.0;
    while index > 0 {
        fraction /= base as f64;
        result += fraction * (index % base) as f64;
        index /= base;
    }
    result
}

fn clamp_probability(value: f64) -> f64 {
    value.clamp(1e-12, 1.0 - 1e-12)
}

fn linspace(start: f64, end: f64, len: usize) -> Result<Vec<f64>, DeseqError> {
    if len == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "low-df prior variance grid".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !start.is_finite() || !end.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "low-df prior variance grid endpoint".to_string(),
            index: None,
            value: if start.is_finite() { end } else { start },
        });
    }
    if len == 1 {
        return Ok(vec![start]);
    }
    let span =
        checked::try_sub(end, start).nonfinite_at("low-df prior variance grid span", None)?;
    let step = checked::try_div(span, len as f64 - 1.0)
        .nonfinite_at("low-df prior variance grid step", None)?;
    (0..len)
        .map(|idx| {
            let offset = checked::try_mul(step, idx as f64)
                .nonfinite_at("low-df prior variance grid offset", None)?;
            checked::try_add(start, offset)
                .nonfinite_at("low-df prior variance grid value", None)
                .map_err(|_| DeseqError::NonFiniteValue {
                    context: "low-df prior variance grid value".to_string(),
                    index: Some(idx),
                    value: start + offset,
                })
        })
        .collect()
}

fn local_linear_smoothed_argmin(x: &[f64], y: &[f64], fine_x: &[f64]) -> Result<f64, DeseqError> {
    if x.len() != y.len() {
        return Err(invalid_dimensions("low-df loess grid", x.len(), y.len()));
    }
    if x.is_empty() || fine_x.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "low-df loess grid".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let span_points = ((x.len() as f64 * LOW_DF_LOESS_SPAN).ceil() as usize).clamp(2, x.len());
    let mut best_x = fine_x[0];
    let mut best_y = f64::INFINITY;
    for target in fine_x.iter().copied() {
        let mut distances = x
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, value)| ((value - target).abs(), idx))
            .collect::<Vec<_>>();
        distances.sort_by(|left, right| left.0.total_cmp(&right.0));
        let max_distance = distances[span_points - 1].0.max(f64::EPSILON);
        let mut sw = 0.0;
        let mut sx = 0.0;
        let mut sy = 0.0;
        let mut sxx = 0.0;
        let mut sxy = 0.0;
        for (_, idx) in distances.iter().copied().take(span_points) {
            let distance = checked::try_sub(x[idx], target)
                .nonfinite_at("low-df prior variance smoother distance", None)?
                .abs();
            let scaled = checked::try_div(distance, max_distance)
                .nonfinite_at("low-df prior variance smoother scaled distance", None)?
                .min(1.0);
            let weight = tricube_weight(scaled)?;
            sw = checked::try_add(sw, weight)
                .nonfinite_at("low-df prior variance smoother weights", None)?;
            sx = checked::try_add(
                sx,
                checked::try_mul(weight, x[idx])
                    .nonfinite_at("low-df prior variance smoother x", None)?,
            )
            .nonfinite_at("low-df prior variance smoother x", None)?;
            sy = checked::try_add(
                sy,
                checked::try_mul(weight, y[idx])
                    .nonfinite_at("low-df prior variance smoother y", None)?,
            )
            .nonfinite_at("low-df prior variance smoother y", None)?;
            let weighted_x = checked::try_mul(weight, x[idx])
                .nonfinite_at("low-df prior variance smoother xx", None)?;
            sxx = checked::try_add(
                sxx,
                checked::try_mul(weighted_x, x[idx])
                    .nonfinite_at("low-df prior variance smoother xx", None)?,
            )
            .nonfinite_at("low-df prior variance smoother xx", None)?;
            sxy = checked::try_add(
                sxy,
                checked::try_mul(weighted_x, y[idx])
                    .nonfinite_at("low-df prior variance smoother xy", None)?,
            )
            .nonfinite_at("low-df prior variance smoother xy", None)?;
        }
        if sw <= 0.0 {
            continue;
        }
        let denominator = checked::try_sub(
            checked::try_mul(sw, sxx)
                .nonfinite_at("low-df prior variance smoother denominator", None)?,
            checked::try_mul(sx, sx)
                .nonfinite_at("low-df prior variance smoother denominator", None)?,
        )
        .nonfinite_at("low-df prior variance smoother denominator", None)?;
        let predicted = if denominator.abs() > f64::EPSILON {
            let numerator = checked::try_sub(
                checked::try_mul(sw, sxy)
                    .nonfinite_at("low-df prior variance smoother slope", None)?,
                checked::try_mul(sx, sy)
                    .nonfinite_at("low-df prior variance smoother slope", None)?,
            )
            .nonfinite_at("low-df prior variance smoother slope", None)?;
            let slope = checked::try_div(numerator, denominator)
                .nonfinite_at("low-df prior variance smoother slope", None)?;
            let intercept = checked::try_div(
                checked::try_sub(
                    sy,
                    checked::try_mul(slope, sx)
                        .nonfinite_at("low-df prior variance smoother intercept", None)?,
                )
                .nonfinite_at("low-df prior variance smoother intercept", None)?,
                sw,
            )
            .nonfinite_at("low-df prior variance smoother intercept", None)?;
            checked::try_add(
                intercept,
                checked::try_mul(slope, target)
                    .nonfinite_at("low-df prior variance smoother prediction", None)?,
            )
            .nonfinite_at("low-df prior variance smoother prediction", None)?
        } else {
            checked::try_div(sy, sw).nonfinite_at("low-df prior variance smoother mean", None)?
        };
        if predicted < best_y {
            best_y = predicted;
            best_x = target;
        }
    }
    if !best_y.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "low-df prior variance smoothing failed".to_string(),
        });
    }
    Ok(best_x)
}

fn low_df_simulated_residual(base: f64, variance: f64, normal: f64) -> Result<f64, DeseqError> {
    checked::try_add(
        base,
        checked::try_mul(variance.sqrt(), normal)
            .nonfinite_at("low-df simulated residual scale", None)?,
    )
    .nonfinite_at("low-df simulated residual", None)
}

fn tricube_weight(scaled: f64) -> Result<f64, DeseqError> {
    let scaled_square = checked::try_mul(scaled, scaled)
        .nonfinite_at("low-df prior variance smoother scaled square", None)?;
    let scaled_cubed = checked::try_mul(scaled_square, scaled)
        .nonfinite_at("low-df prior variance smoother scaled cube", None)?;
    let one_minus_scaled_cubed = checked::try_sub(1.0, scaled_cubed)
        .nonfinite_at("low-df prior variance smoother tricube complement", None)?;
    let weight_square = checked::try_mul(one_minus_scaled_cubed, one_minus_scaled_cubed)
        .nonfinite_at("low-df prior variance smoother tricube weight", None)?;
    checked::try_mul(weight_square, one_minus_scaled_cubed)
        .nonfinite_at("low-df prior variance smoother tricube weight", None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_df_smoother_rejects_overflowed_weighted_accumulation() {
        let err =
            local_linear_smoothed_argmin(&[2e154, 2.000000000000001e154], &[1.0, 2.0], &[2e154])
                .unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, .. }
                if context == "low-df prior variance smoother xx"
        ));
    }

    #[test]
    fn low_df_smoother_rejects_overflowed_response_accumulation() {
        let x = (0..20).map(|idx| idx as f64 * 0.1).collect::<Vec<_>>();
        let y = vec![1e308; x.len()];
        let err = local_linear_smoothed_argmin(&x, &y, &[0.0]).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, .. }
                if context == "low-df prior variance smoother y"
        ));
    }

    #[test]
    fn low_df_simulated_residual_rejects_overflowed_scale() {
        let err = low_df_simulated_residual(0.0, f64::MAX, f64::MAX).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, .. }
                if context == "low-df simulated residual scale"
        ));
    }

    #[test]
    fn low_df_linspace_rejects_nonfinite_grid_arithmetic() {
        let endpoint_err = linspace(f64::INFINITY, 1.0, 3).unwrap_err();
        assert!(matches!(
            endpoint_err,
            DeseqError::NonFiniteValue { context, .. }
                if context == "low-df prior variance grid endpoint"
        ));

        let span_err = linspace(-f64::MAX, f64::MAX, 3).unwrap_err();
        assert!(matches!(
            span_err,
            DeseqError::NonFiniteValue { context, .. }
                if context == "low-df prior variance grid span"
        ));

        assert_eq!(linspace(2.0, 2.0, 3).unwrap(), vec![2.0, 2.0, 2.0]);
    }

    #[test]
    fn checked_div_rejects_zero_and_nonfinite_prior_arithmetic() {
        let zero_err = checked::try_div(1.0, 0.0)
            .nonfinite_at("test low-df division", None)
            .unwrap_err();
        assert!(matches!(
            zero_err,
            DeseqError::NonFiniteValue { context, .. } if context == "test low-df division"
        ));

        let nonfinite_err = checked::try_div(f64::NAN, 1.0)
            .nonfinite_at("test low-df division", None)
            .unwrap_err();
        assert!(matches!(
            nonfinite_err,
            DeseqError::NonFiniteValue { context, .. } if context == "test low-df division"
        ));
    }

    #[test]
    fn tricube_weight_matches_boundary_values() {
        assert_eq!(tricube_weight(0.0).unwrap(), 1.0);
        assert_eq!(tricube_weight(1.0).unwrap(), 0.0);
    }

    #[test]
    fn kl_divergence_rejects_overflowed_bin_adjustment() {
        let err = kl_divergence(&[f64::MAX], &[1.0]).unwrap_err();

        assert!(matches!(err, DeseqError::NonFiniteValue { .. }));
    }

    #[test]
    fn kl_divergence_keeps_large_finite_bins_finite() {
        let divergence = kl_divergence(&[1e100, 2e100], &[2e100, 1e100]).unwrap();

        assert!(divergence.is_finite());
    }
}
