//! Builder gene wise stages.
use super::*;

impl DeseqBuilder {
    /// Leaf workflow behind [`DeseqBuilder::fit_gene_wise_dispersions_linear_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_gene_wise_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_gene_wise_dispersions_linear_mu(counts, design));
        }

        self.ensure_no_observation_weights("native linear-mu dispersion estimation")?;
        design.validate_full_rank("linear-mu dispersion")?;
        let stages = self.normalization_stages(counts)?;
        let dispersion = estimate_gene_wise_dispersions_linear_mu(
            GeneWiseDispersionInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                normalized_counts: &stages.normalized,
                base_mean: &stages.base_mean,
                base_var: &stages.base_var,
                all_zero: &stages.all_zero,
                observation_weights: None,
            },
            self.gene_wise_dispersion_options,
        )?;
        let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_gene_wise_dispersions_glm_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_gene_wise_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_gene_wise_dispersions_glm_mu(counts, design));
        }

        design.validate_full_rank("GLM-mu dispersion")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let estimator = if self.automatic_dispersion_means {
            estimate_gene_wise_dispersions_auto
        } else {
            estimate_gene_wise_dispersions_glm_mu
        };
        let dispersion = estimator(
            GeneWiseDispersionInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                normalized_counts: &stages.normalized,
                base_mean: &stages.base_mean,
                base_var: &stages.base_var,
                all_zero: &stages.all_zero,
                observation_weights: stages.observation_weights.as_ref(),
            },
            self.gene_wise_dispersion_options,
            self.irls_options.clone(),
        )?;
        let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
    }
}
