//! Dispersion estimation stages and DESeq2-shaped low-level primitives.

pub mod gene_est;
pub mod map;
pub mod prior;
pub mod trend;

pub use gene_est::{
    DispersionLineSearchOutput, DispersionPrior, GeneWiseDispersionFitMethod,
    GeneWiseDispersionInput, GeneWiseDispersionOptions, GeneWiseDispersionOutput,
    WeightedDispersionFitInput, cox_reid_adjustment, cox_reid_adjustment_derivative,
    cox_reid_adjustment_derivative_weighted, cox_reid_adjustment_second_derivative,
    cox_reid_adjustment_second_derivative_weighted, cox_reid_adjustment_weighted,
    dispersion_log_posterior, dispersion_log_posterior_derivative,
    dispersion_log_posterior_derivative_with_prior,
    dispersion_log_posterior_derivative_with_prior_and_weights,
    dispersion_log_posterior_second_derivative,
    dispersion_log_posterior_second_derivative_with_prior,
    dispersion_log_posterior_second_derivative_with_prior_and_weights,
    dispersion_log_posterior_with_prior, dispersion_log_posterior_with_prior_and_weights,
    dispersion_nb_log_likelihood_kernel, dispersion_nb_log_likelihood_kernel_derivative,
    dispersion_nb_log_likelihood_kernel_derivative_weighted,
    dispersion_nb_log_likelihood_kernel_second_derivative,
    dispersion_nb_log_likelihood_kernel_second_derivative_weighted,
    dispersion_nb_log_likelihood_kernel_weighted, dispersion_prior_derivative,
    dispersion_prior_log_density, dispersion_prior_second_derivative,
    estimate_gene_wise_dispersions_glm_mu, estimate_gene_wise_dispersions_linear_mu,
    fit_dispersion_grid, fit_dispersion_grid_no_cr, fit_dispersion_grid_no_cr_with_prior,
    fit_dispersion_grid_no_cr_with_prior_and_weights, fit_dispersion_grid_with_prior,
    fit_dispersion_grid_with_prior_and_weights, fit_dispersion_line_search,
    fit_dispersion_line_search_no_cr, fit_dispersion_line_search_no_cr_with_prior,
    fit_dispersion_line_search_no_cr_with_prior_and_weights, fit_dispersion_line_search_with_prior,
    fit_dispersion_line_search_with_prior_and_weights, initial_dispersion_estimates,
    linear_model_mu, moments_dispersion_estimates,
    moments_dispersion_estimates_with_normalization_factors, rough_dispersion_estimates,
};

pub use prior::{
    DispersionPriorVarianceOutput, estimate_dispersion_prior, estimate_dispersion_prior_variance,
    estimate_low_df_prior_variance, log_dispersion_residuals_above_min, mad_squared,
};

pub use map::{
    MapDispersionInput, MapDispersionOptions, MapDispersionOutput, estimate_map_dispersions,
    map_dispersion_initial_value, map_dispersion_outlier,
};

pub use trend::{
    DispersionTrendFit, LocalDispersionTrend, LocalDispersionTrendFit, LocalDispersionTrendOptions,
    MeanDispersionTrend, MeanDispersionTrendFit, MeanDispersionTrendOptions,
    ParametricDispersionTrend, ParametricDispersionTrendFit, ParametricDispersionTrendOptions,
    fit_dispersion_trend, fit_local_dispersion_trend, fit_mean_dispersion_trend,
    fit_parametric_dispersion_trend, local_trend_use_for_fit, mean_trend_use_for_fit,
    mean_trend_use_for_mean, parametric_trend_use_for_fit,
};
