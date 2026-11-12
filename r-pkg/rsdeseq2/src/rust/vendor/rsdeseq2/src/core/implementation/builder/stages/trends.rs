//! Builder trends stages.
use super::*;

impl DeseqBuilder {
    /// Run linear-mu gene-wise dispersion estimation and fit the parametric trend.
    ///
    /// This mirrors the implemented subset of `estimateDispersionsGeneEst`
    /// followed by `estimateDispersionsFit(fitType="parametric")`. It fills
    /// `dispGeneEst`, `dispFit`, and the linear fitted mean matrix, but it does
    /// not yet estimate prior variance or MAP dispersions.
    pub fn fit_parametric_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_parametric_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the mean trend.
    ///
    /// This mirrors the implemented subset of `estimateDispersionsGeneEst`
    /// followed by `estimateDispersionsFit(fitType="mean")`. It fills
    /// `dispGeneEst`, a constant `dispFit` for non-all-zero rows, and the
    /// linear fitted mean matrix, but it does not yet estimate prior variance
    /// or MAP dispersions for mean-trend fits.
    pub fn fit_mean_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_mean_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the parametric trend.
    ///
    /// This mirrors the current non-`linearMu` gene-wise branch
    /// followed by `estimateDispersionsFit(fitType="parametric")`.
    pub fn fit_parametric_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_parametric_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the mean trend.
    pub fn fit_mean_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_mean_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_local_dispersion_trend_linear_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_local_dispersion_trend_glm_mu(counts, design));
        }

        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    pub(super) fn attach_parametric_dispersion_trend(
        &self,
        mut fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_parametric_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            ParametricDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..ParametricDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Parametric(trend_fit));
        Ok(fit)
    }

    pub(super) fn attach_mean_dispersion_trend(
        &self,
        mut fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_mean_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            MeanDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..MeanDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Mean(trend_fit));
        Ok(fit)
    }

    pub(super) fn attach_local_dispersion_trend(
        &self,
        mut fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_local_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            LocalDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..LocalDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Local(trend_fit));
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

    /// Run the implemented linear-mu dispersion trend path selected by `fit_type`.
    ///
    /// `Parametric`, `Local`, and `Mean` are currently implemented.
    /// `GlmGamPoi` returns `UnsupportedFeature` until a parity implementation
    /// is added.
    pub fn fit_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_dispersion_trend_linear_mu(counts, design));
        }

        match self.fit_type {
            FitType::Parametric => self.fit_parametric_dispersion_trend_linear_mu(counts, design),
            FitType::Mean => self.fit_mean_dispersion_trend_linear_mu(counts, design),
            FitType::Local => self.fit_local_dispersion_trend_linear_mu(counts, design),
            FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
                feature: "linear-mu glmGamPoi dispersion trend fitting".to_string(),
            }),
        }
    }

    /// Run the implemented GLM-mu dispersion trend path selected by `fit_type`.
    ///
    /// `Parametric`, `Local`, and `Mean` are currently implemented.
    /// `GlmGamPoi` returns `UnsupportedFeature` until a parity implementation
    /// is added.
    pub fn fit_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_dispersion_trend_glm_mu(counts, design));
        }

        match self.fit_type {
            FitType::Parametric => self.fit_parametric_dispersion_trend_glm_mu(counts, design),
            FitType::Mean => self.fit_mean_dispersion_trend_glm_mu(counts, design),
            FitType::Local => self.fit_local_dispersion_trend_glm_mu(counts, design),
            FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
                feature: "GLM-mu glmGamPoi dispersion trend fitting".to_string(),
            }),
        }
    }
}
