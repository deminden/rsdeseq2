//! Builder gene wise stages.
use super::*;
use crate::dispersion::gene_est::{
    estimate_gene_wise_dispersions_for_workflow,
    estimate_gene_wise_dispersions_linear_mu_with_moments,
};

impl DeseqBuilder {
    pub(super) fn leaf_gene_wise_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.in_pool(|| {
            self.ensure_no_observation_weights("native linear-mu dispersion estimation")?;
            design.validate_full_rank("linear-mu dispersion")?;
            let stages = self.normalization_stages(counts)?;
            let dispersion = estimate_gene_wise_dispersions_linear_mu_with_moments(
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
                self.formula_sample_factors,
            )?;
            let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
            fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
            fit.disp_gene_iter = Some(dispersion.disp_iter);
            fit.dispersion_converged = Some(dispersion.converged);
            fit.mu = Some(dispersion.mu);
            Ok(fit)
        })
    }

    pub(super) fn leaf_gene_wise_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.in_pool(|| {
            design.validate_full_rank("GLM-mu dispersion")?;
            let stages = self.normalization_stages_for_design(counts, design)?;
            let dispersion = estimate_gene_wise_dispersions_for_workflow(
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
                self.automatic_dispersion_means,
                self.formula_sample_factors,
            )?;
            let mut fit = self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
            fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
            fit.disp_gene_iter = Some(dispersion.disp_iter);
            fit.dispersion_converged = Some(dispersion.converged);
            fit.mu = Some(dispersion.mu);
            Ok(fit)
        })
    }
}
