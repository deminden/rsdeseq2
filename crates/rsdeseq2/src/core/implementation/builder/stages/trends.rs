//! Builder trends stages.
use super::*;

impl DeseqBuilder {
    /// Leaf workflow behind [`DeseqBuilder::fit_parametric_dispersion_trend_linear_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_parametric_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.leaf_parametric_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_mean_dispersion_trend_linear_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_mean_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_mean_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_parametric_dispersion_trend_glm_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_parametric_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_parametric_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_mean_dispersion_trend_glm_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_mean_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_mean_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_local_dispersion_trend_linear_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_local_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_local_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_local_dispersion_trend_glm_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_local_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_local_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.leaf_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    pub(super) fn attach_parametric_dispersion_trend(
        &self,
        fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        self.attach_estimated_dispersion_trend(fit, FitType::Parametric)
    }

    pub(super) fn attach_mean_dispersion_trend(
        &self,
        fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        self.attach_estimated_dispersion_trend(fit, FitType::Mean)
    }

    pub(super) fn attach_local_dispersion_trend(
        &self,
        fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        self.attach_estimated_dispersion_trend(fit, FitType::Local)
    }

    /// DESeq2 `estimateDispersionsFit`: shared row rule, parametric-to-local
    /// fallback and local fit on `useForFit` rows only.
    fn attach_estimated_dispersion_trend(
        &self,
        mut fit: DeseqFit,
        fit_type: FitType,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend = estimate_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            fit_type,
            self.gene_wise_dispersion_options.min_disp,
        )?;
        fit.disp_fit = Some(trend.disp_fit().to_vec());
        fit.dispersion_trend = Some(trend);
        Ok(fit)
    }

    pub(super) fn attach_existing_dispersion_trend(
        &self,
        mut fit: DeseqFit,
        trend: &DispersionTrendFit,
    ) -> Result<DeseqFit, DeseqError> {
        fit.disp_fit = Some(trend.evaluate_many_allow_missing(&fit.base_mean)?);
        fit.dispersion_trend = Some(trend.clone());
        Ok(fit)
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_dispersion_trend_linear_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_dispersion_trend_linear_mu(counts, design));
        }

        match self.fit_type {
            FitType::Parametric => self.leaf_parametric_dispersion_trend_linear_mu(counts, design),
            FitType::Mean => self.leaf_mean_dispersion_trend_linear_mu(counts, design),
            FitType::Local => self.leaf_local_dispersion_trend_linear_mu(counts, design),
        }
    }

    /// Leaf workflow behind [`DeseqBuilder::fit_dispersion_trend_glm_mu`] and [`DeseqBuilder::run`].
    pub(super) fn leaf_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.leaf_dispersion_trend_glm_mu(counts, design));
        }

        match self.fit_type {
            FitType::Parametric => self.leaf_parametric_dispersion_trend_glm_mu(counts, design),
            FitType::Mean => self.leaf_mean_dispersion_trend_glm_mu(counts, design),
            FitType::Local => self.leaf_local_dispersion_trend_glm_mu(counts, design),
        }
    }
}
