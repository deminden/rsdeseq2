//! Complete location-scale shrinkage preparation and summaries in Rust.
//! The default public adapter uses this engine; advanced option gates remain.
use super::location_missing::{NA, excluded, is_na, log_marginal};
use super::{
    AdaptiveOutputSelection, LocationLikelihood, MixtureComponents, MixtureFit, MixtureFitOptions,
    MixtureGrid, MixtureGridOptions,
};
use crate::DeseqError;

/// Batch controls. The default matches the wrapper's normal-mixture shrink
/// configuration; it is distinct from the original ashr function's defaults.
#[derive(Clone, Debug)]
pub struct LocationShrinkageOptions {
    pub grid: MixtureGridOptions,
    pub supplied_grid: Option<MixtureGrid>,
    pub fixed_prior: bool,
    pub likelihood: LocationLikelihood,
    pub alpha: f64,
    pub prune_threshold: f64,
    pub observation_weights: Option<Vec<f64>>,
    pub fit_options: MixtureFitOptions,
    /// Raw posterior CDF cutoffs, evaluated with the unpruned fitted prior.
    pub thresholds: Vec<f64>,
}
impl Default for LocationShrinkageOptions {
    fn default() -> Self {
        Self {
            grid: MixtureGridOptions::default(),
            supplied_grid: None,
            fixed_prior: false,
            likelihood: LocationLikelihood::Normal,
            alpha: 0.0,
            prune_threshold: 1e-10,
            observation_weights: None,
            fit_options: MixtureFitOptions {
                likelihood_offset: 1e-6,
                active_iterations: 20,
                positive_coordinate_convergence: true,
                sqp: super::MixtureSqpOptions {
                    ashr_compatibility: true,
                    automatic_active_iterations: true,
                    search_direction_threshold: Some(1e-14),
                    ..Default::default()
                },
                ..Default::default()
            },
            thresholds: vec![],
        }
    }
}
/// Full batch summaries, fitted prior and explicit diagnostics. Missing summaries
/// retain their floating-point missing values rather than being rounded in R.
#[derive(Clone, Debug)]
pub struct LocationShrinkage {
    pub fitted_grid: MixtureGrid,
    pub fitted_weights: Vec<f64>,
    pub fit: Option<MixtureFit>,
    pub posterior_mean: Vec<f64>,
    pub posterior_sd: Vec<f64>,
    pub negative_probability: Vec<f64>,
    pub positive_probability: Vec<f64>,
    pub local_false_discovery_rate: Vec<f64>,
    pub local_false_sign_rate: Vec<f64>,
    pub svalues: Vec<f64>,
    pub qvalues: Vec<f64>,
    /// Observation-major raw posterior CDF values, in requested cutoff order.
    pub cdf: Vec<f64>,
    pub log_likelihood: f64,
    pub warnings: Vec<String>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("location shrinkage: {reason}"),
    }
}

fn zero_clip(value: f64) -> f64 {
    if value < 0.0 { 0.0 } else { value }
}

/// Ordered cumulative local-error averages, including the original ashr
/// missing-entry replacement recycling. This compatibility operation is kept
/// separate from `svalues`, whose missing values remain missing by contract.
fn legacy_qvalues(rates: &[f64], warnings: &mut Vec<String>) -> Vec<f64> {
    let mut order = (0..rates.len()).collect::<Vec<_>>();
    order.sort_by(|&i, &j| match (rates[i].is_nan(), rates[j].is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => rates[i].partial_cmp(&rates[j]).unwrap(),
    });
    let count = order.iter().take_while(|&&i| !rates[i].is_nan()).count();
    if count == 0 {
        return vec![f64::from_bits(0x7ff0_0000_0000_07a2); rates.len()];
    }
    let mut cumulative = 0.0;
    let averages = order[..count]
        .iter()
        .enumerate()
        .map(|(rank, &i)| {
            cumulative += rates[i];
            cumulative / (rank + 1) as f64
        })
        .collect::<Vec<_>>();
    if !rates.len().is_multiple_of(count) {
        warnings.push("number of items to replace is not a multiple of replacement length".into());
    }
    let mut result = vec![0.0; rates.len()];
    for (rank, &i) in order.iter().enumerate() {
        result[i] = averages[rank % count];
    }
    result
}

/// Prepare, fit and summarize a complete observation batch without an R runtime.
/// Excluded observations use prior summaries; raw requested CDFs retain their
/// separate missing/zero/infinite-SE policy. Supplied grids permit fixed priors.
pub fn shrink_location_batch(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
) -> Result<LocationShrinkage, DeseqError> {
    shrink_location_batch_with_outputs(
        estimates,
        standard_errors,
        options,
        AdaptiveOutputSelection::level(2.0)?,
    )
}
pub(super) fn shrink_location_batch_with_outputs(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
    outputs: AdaptiveOutputSelection,
) -> Result<LocationShrinkage, DeseqError> {
    let super::location_fit::FittedLocation {
        grid,
        fitted_weights,
        fit,
        x,
        se,
        scale,
        mut warnings,
    } = super::location_fit::fit_location(estimates, standard_errors, options)?;
    let n = estimates.len();
    if matches!(&options.likelihood,LocationLikelihood::Builtin(b) if matches!(b.name(),"logF"|"pois"|"binom"))
        && x.iter().chain(&se).any(|v| v.is_nan())
    {
        return Err(invalid(
            "this likelihood requires nonmissing observation bounds",
        ));
    }
    let rows = 0..n;
    let raw_prior = grid
        .prior_with_likelihood(&fitted_weights, &options.likelihood)?
        .orient_for(&x, &se, rows.clone())?;
    if let MixtureComponents::TruncatedNormal { lower, upper, .. } = &grid.components
        && options.thresholds.iter().any(|&cutoff| {
            !lower
                .iter()
                .zip(upper)
                .any(|(&a, &b)| a <= cutoff && cutoff <= b)
        })
    {
        return Err(invalid(
            "truncated-normal posterior CDF requires a component containing each cutoff",
        ));
    }

    let keep = fitted_weights
        .iter()
        .enumerate()
        .filter_map(|(i, &w)| (w > options.prune_threshold).then_some(i))
        .collect::<Vec<_>>();
    if keep.is_empty() {
        return Err(invalid("pruning removed every prior component"));
    }
    let mut weights = keep.iter().map(|&i| fitted_weights[i]).collect::<Vec<_>>();
    let total: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= total;
    }
    let summary_grid = MixtureGrid {
        components: grid.components.select(&keep),
        initial_weights: weights.clone(),
        dirichlet: vec![1.0; keep.len()],
        warnings: vec![],
    };
    let summary_prior = summary_grid
        .prior_with_likelihood(&weights, &options.likelihood)?
        .orient_for(&x, &se, rows.clone())?;
    let (prior_mean, prior_sd, zero) = summary_grid.components.moments(&weights);
    let prior_negative = summary_grid.components.prior_cdf(&weights, 0.0)?;
    let prior_zero: f64 = weights
        .iter()
        .zip(&zero)
        .filter(|(_, z)| **z)
        .map(|(w, _)| *w)
        .sum();
    let mut means = vec![prior_mean; n];
    let mut sds = vec![prior_sd; n];
    let mut negative = vec![prior_negative; n];
    let mut lfdr = vec![prior_zero; n];
    let mut cdfs = Vec::with_capacity(
        n.checked_mul(options.thresholds.len())
            .ok_or_else(|| invalid("CDF matrix shape overflows"))?,
    );
    let mut log_likelihood = 0.0;
    let missing_rows = x
        .iter()
        .zip(&se)
        .filter(|(x, s)| x.is_nan() || s.is_nan())
        .count();
    let mut missing_row = 0;
    for i in 0..n {
        let included = se[i] > 0.0 && se[i].is_finite() && x[i].is_finite();
        if included {
            let posterior = summary_prior.posterior_at(x[i], se[i], i, n)?;
            means[i] = posterior.mean();
            sds[i] = posterior.standard_deviation();
            lfdr[i] = zero_clip(
                posterior
                    .weights()
                    .iter()
                    .zip(&zero)
                    .filter(|(_, z)| **z)
                    .map(|(w, _)| *w)
                    .sum(),
            );
            negative[i] = zero_clip(posterior.cdf(0.0) - lfdr[i]);
            cdfs.extend(raw_prior.cdf_values_at(x[i], se[i], &options.thresholds, i, n)?);
            let ll = raw_prior.log_likelihoods_at(x[i], se[i], i, n)?;
            log_likelihood += log_marginal(&ll, &fitted_weights)
                - options.alpha
                    * standard_errors[if standard_errors.len() == 1 { 0 } else { i }].ln();
        } else {
            let (raw_cdfs, ll) = excluded(
                &grid.components,
                &fitted_weights,
                x[i],
                se[i],
                &options.thresholds,
                missing_row,
                missing_rows,
            );
            if x[i].is_nan() || se[i].is_nan() {
                missing_row += 1;
            }
            cdfs.extend(raw_cdfs);
            log_likelihood += ll
                - options.alpha
                    * standard_errors[if standard_errors.len() == 1 { 0 } else { i }].ln();
        }
        means[i] *= scale[i];
        sds[i] *= scale[i];
    }
    if x.iter().chain(&se).any(|v| is_na(*v)) {
        log_likelihood = NA;
    }
    let positive = negative
        .iter()
        .zip(&lfdr)
        .map(|(n, z)| zero_clip(1.0 - n - z))
        .collect::<Vec<_>>();
    let lfsr = negative
        .iter()
        .zip(&lfdr)
        .map(|(&n, &z)| zero_clip(if n > 0.5 * (1.0 - z) { 1.0 - n } else { n + z }))
        .collect::<Vec<_>>();
    let mut sign_warnings = vec![];
    let mut null_warnings = vec![];
    let svalues = if outputs.svalues {
        legacy_qvalues(&lfsr, &mut sign_warnings)
    } else {
        vec![]
    };
    let qvalues = if outputs.qvalues {
        legacy_qvalues(&lfdr, &mut null_warnings)
    } else {
        vec![]
    };
    if let Some(f) = super::likelihood_diagnostics::WarningLikelihood::new(&options.likelihood) {
        let raw = super::likelihood_diagnostics::Diagnostics::new(f, &grid.components, &x, &se);
        if outputs.flash_data {
            warnings.extend_from_slice(&raw.likelihood);
            warnings.extend_from_slice(&raw.mean);
            warnings.extend_from_slice(&raw.second);
        }
        if outputs.log_likelihood {
            warnings.extend_from_slice(&raw.likelihood);
        }
        if outputs.log_likelihood_ratio {
            warnings.extend_from_slice(&raw.likelihood);
            warnings.extend_from_slice(&raw.null);
        }
        let summary =
            super::likelihood_diagnostics::Diagnostics::new(f, &summary_grid.components, &x, &se);
        summary.append_selected_summary(&mut warnings, outputs, &sign_warnings, &null_warnings);
        for &cutoff in &options.thresholds {
            warnings.extend_from_slice(&raw.likelihood);
            warnings.extend(super::likelihood_diagnostics::cdf(
                f,
                &grid.components,
                &x,
                &se,
                cutoff,
            ));
        }
    } else {
        warnings.extend(sign_warnings);
        warnings.extend(null_warnings);
    }
    Ok(LocationShrinkage {
        fitted_grid: grid,
        fitted_weights,
        fit,
        posterior_mean: means,
        posterior_sd: sds,
        negative_probability: negative,
        positive_probability: positive,
        local_false_discovery_rate: lfdr,
        local_false_sign_rate: lfsr,
        svalues,
        qvalues,
        cdf: cdfs,
        log_likelihood,
        warnings,
    })
}
