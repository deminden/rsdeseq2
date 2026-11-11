//! Shared location preprocessing and likelihood-only mixture fitting.
use super::{
    LocationLikelihood, LocationShrinkageOptions, MixtureComponents, MixtureFit, MixtureGrid,
    prepare_mixture_grid,
};
use crate::DeseqError;

pub(super) struct FittedLocation {
    pub grid: MixtureGrid,
    pub fitted_weights: Vec<f64>,
    pub fit: Option<MixtureFit>,
    pub x: Vec<f64>,
    pub se: Vec<f64>,
    pub scale: Vec<f64>,
    pub warnings: Vec<String>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("location shrinkage: {reason}"),
    }
}
pub(super) fn fit_location(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
) -> Result<FittedLocation, DeseqError> {
    if standard_errors.len() != estimates.len() && standard_errors.len() != 1 {
        return Err(invalid("SEs must have length one or match estimates"));
    }
    if !options.alpha.is_finite()
        || !options.prune_threshold.is_finite()
        || options.prune_threshold < 0.0
        || options.thresholds.iter().any(|v| v.is_nan())
    {
        return Err(invalid(
            "finite alpha/nonnegative pruning threshold and nonmissing CDF cutoffs are required",
        ));
    }
    if !estimates.iter().enumerate().any(|(i, x)| {
        !x.is_nan() && !standard_errors[if standard_errors.len() == 1 { 0 } else { i }].is_nan()
    }) {
        return Err(invalid("all input observations are missing"));
    }
    if options.fixed_prior && options.supplied_grid.is_none() {
        return Err(invalid("fixed_prior requires a supplied mixture grid"));
    }
    let n = estimates.len();
    let mut x = Vec::with_capacity(n);
    let mut se = Vec::with_capacity(n);
    let mut scale = Vec::with_capacity(n);
    let mut included = vec![];
    for (i, &estimate) in estimates.iter().enumerate() {
        let s = standard_errors[if standard_errors.len() == 1 { 0 } else { i }];
        if s < 0.0 {
            return Err(invalid("standard errors must be nonnegative"));
        }
        let factor = s.powf(options.alpha);
        let value = estimate / factor;
        let error = s.powf(1.0 - options.alpha);
        if !value.is_nan() && !error.is_nan() && error != 0.0 && error != f64::INFINITY {
            if !value.is_finite() || !error.is_finite() {
                return Err(invalid("nonexcluded estimates and SEs must be finite"));
            }
            included.push(i);
        }
        x.push(value);
        se.push(error);
        scale.push(factor);
    }
    let grid = match &options.supplied_grid {
        Some(g) => g.clone(),
        None => prepare_mixture_grid(&x, &se, &options.grid)?,
    };
    if let MixtureComponents::TruncatedNormal { sds, .. } = &grid.components
        && sds.contains(&0.0)
        && x.iter().zip(&se).any(|(x, s)| x.is_nan() || s.is_nan())
    {
        return Err(invalid(
            "missing observations with zero-SD truncated components are undefined",
        ));
    }
    if grid.dirichlet.len() != grid.initial_weights.len()
        || grid.dirichlet.iter().any(|v| {
            !v.is_finite()
                || *v <= 0.0
                || (!options.fit_options.optimizer.is_variational() && *v < 1.0)
        })
    {
        return Err(invalid(
            "finite Dirichlet parameters >= 1 must match the prior components",
        ));
    }
    let mut warnings = grid.warnings.clone();
    let initial_prior = grid.prior_with_likelihood(&grid.initial_weights, &options.likelihood)?;
    if let LocationLikelihood::Restricted(f) = &options.likelihood {
        warnings.extend(super::likelihood_diagnostics::likelihood(
            f,
            &grid.components,
            &x,
            &se,
        ));
    }
    let fit = if options.fixed_prior {
        None
    } else {
        let fit_x = included.iter().map(|&i| x[i]).collect::<Vec<_>>();
        let fit_se = included.iter().map(|&i| se[i]).collect::<Vec<_>>();
        let weights = if let Some(w) = &options.observation_weights {
            if w.len() != n {
                return Err(invalid(
                    "observation weights must match the complete observation vector",
                ));
            }
            Some(included.iter().map(|&i| w[i]).collect::<Vec<_>>())
        } else {
            None
        };
        // Use the likelihood selected above for every component family.
        let k = grid.initial_weights.len();
        fit_x
            .len()
            .checked_mul(k)
            .ok_or_else(|| invalid("likelihood matrix shape overflows"))?;
        let mut ll = vec![];
        for (&value, &error) in fit_x.iter().zip(&fit_se) {
            ll.extend(initial_prior.log_likelihoods(value, error)?);
        }
        let f = super::fit_mixture_proportions(
            crate::matrix::MatrixView::new(
                fit_x.len(),
                k,
                &ll,
                crate::matrix::MatrixOrder::RowMajor,
            )?,
            weights.as_deref(),
            &grid.dirichlet,
            Some(&grid.initial_weights),
            options.fit_options,
        )?;
        if !f.converged {
            warnings.push("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.".into());
        }
        Some(f)
    };
    let fitted_weights = fit
        .as_ref()
        .map(|f| f.weights.clone())
        .unwrap_or_else(|| grid.initial_weights.clone());
    Ok(FittedLocation {
        grid,
        fitted_weights,
        fit,
        x,
        se,
        scale,
        warnings,
    })
}
