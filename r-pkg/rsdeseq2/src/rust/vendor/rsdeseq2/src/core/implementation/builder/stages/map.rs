//! Builder map stages.
use super::*;

impl DeseqBuilder {
    /// Run linear-mu gene-wise, selected trend, prior variance, and MAP dispersion stages.
    ///
    /// This fills final `dispersion` values using the builder's `fit_type`.
    /// Parametric, local, and mean trends are implemented. It follows the implemented
    /// subset of DESeq2's
    /// `estimateDispersionsMAP(type="DESeq2")`: no observation weights and
    /// deterministic prior-variance estimation, including the low-df
    /// histogram/KL branch.
    pub fn fit_map_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_map_dispersions_linear_mu(counts, design));
        }

        let fit = self.fit_dispersion_trend_linear_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    /// Run linear-mu gene-wise, parametric trend, prior variance, and MAP dispersion stages.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_map_dispersions_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_map_dispersions_linear_mu_parametric(counts, design));
        }

        let fit = self.fit_parametric_dispersion_trend_linear_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    /// Run GLM-mu gene-wise, selected trend, prior variance, and MAP dispersion stages.
    ///
    /// This fills final `dispersion` values using the builder's `fit_type`.
    /// Parametric, local, and mean trends are implemented. Builder-supplied
    /// observation weights flow through the GLM-mu gene-wise, MAP, and native
    /// Wald stages after DESeq2-style preprocessing.
    pub fn fit_map_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_map_dispersions_glm_mu(counts, design));
        }

        let fit = self.fit_dispersion_trend_glm_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    /// Run GLM-mu gene-wise, parametric trend, prior variance, and MAP dispersion stages.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_map_dispersions_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_map_dispersions_glm_mu_parametric(counts, design));
        }

        let fit = self.fit_parametric_dispersion_trend_glm_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    pub(super) fn fit_map_dispersions_glm_mu_with_dispersion_function(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        trend: &DispersionTrendFit,
        disp_prior_var: f64,
        var_log_disp_estimates: f64,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        let fit = self.attach_existing_dispersion_trend(fit, trend)?;
        self.attach_map_dispersions_with_prior_values(
            counts,
            design,
            fit,
            Some(disp_prior_var),
            Some(var_log_disp_estimates),
        )
    }

    pub(super) fn attach_map_dispersions(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        self.attach_map_dispersions_with_prior_values(counts, design, fit, None, None)
    }

    pub(super) fn attach_map_dispersions_with_prior_values(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        mut fit: DeseqFit,
        supplied_disp_prior_var: Option<f64>,
        supplied_var_log_disp_estimates: Option<f64>,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before MAP fitting".to_string(),
                })?;
        let disp_fit = fit
            .disp_fit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted dispersion trend is required before MAP fitting".to_string(),
            })?;
        let mu = fit
            .mu
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted means are required before MAP fitting".to_string(),
            })?;
        // Replacement refits reuse the original trend/prior scalars; ordinary
        // MAP fitting estimates them from the current gene-wise/trend state.
        let prior_variance = match (supplied_disp_prior_var, supplied_var_log_disp_estimates) {
            (Some(disp_prior_var), Some(var_log_disp_estimates)) => {
                crate::dispersion::DispersionPriorVarianceOutput {
                    disp_prior_var,
                    var_log_disp_estimates,
                    expected_log_dispersion_variance: f64::NAN,
                    residual_degrees_of_freedom: design.n_samples() - design.n_coefficients(),
                    above_min_disp: Vec::new(),
                }
            }
            (None, None) => estimate_dispersion_prior_variance(
                disp_gene_est,
                disp_fit,
                self.gene_wise_dispersion_options.min_disp,
                design.n_samples(),
                design.n_coefficients(),
            )?,
            _ => {
                return Err(DeseqError::InvalidOptions {
                    reason:
                        "dispersion prior variance and log-dispersion variance must be supplied together"
                            .to_string(),
                })
            }
        };
        let map = estimate_map_dispersions(
            MapDispersionInput {
                counts,
                design,
                mu,
                disp_gene_est,
                disp_fit,
                all_zero: &fit.all_zero,
                observation_weights: fit.observation_weights.as_ref(),
                disp_prior_var: prior_variance.disp_prior_var,
                var_log_disp_estimates: prior_variance.var_log_disp_estimates,
            },
            MapDispersionOptions::from(self.gene_wise_dispersion_options),
        )?;
        fit.disp_prior_var = Some(prior_variance.disp_prior_var);
        fit.var_log_disp_estimates = Some(prior_variance.var_log_disp_estimates);
        fit.disp_map = Some(map.disp_map);
        fit.disp_iter = Some(map.disp_iter);
        fit.disp_outlier = Some(map.disp_outlier);
        fit.dispersion_converged = Some(map.converged);
        fit.dispersion = Some(map.dispersion);
        Ok(fit)
    }
}
