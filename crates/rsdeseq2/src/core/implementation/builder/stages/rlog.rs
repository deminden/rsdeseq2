//! Builder rlog stages.
use super::*;

impl DeseqBuilder {
    /// Apply the implemented GLM-mu rlog workflow in one builder call.
    ///
    /// This fits GLM-mu gene-wise dispersions and the selected trend, then
    /// applies the fit-state rlog transform with default IRLS options. As in
    /// DESeq2's `rlogData`, sample effects are fit with the trend values before
    /// MAP shrinkage to final per-gene dispersions.
    pub fn rlog_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.rlog_glm_mu(counts, design));
        }

        self.rlog_glm_mu_with_fit(counts, design)
            .map(|output| output.rlog)
    }

    /// Apply the implemented GLM-mu rlog workflow with explicit rlog IRLS options.
    pub fn rlog_glm_mu_with_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.rlog_glm_mu_with_options(counts, design, rlog_irls_options));
        }

        self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply the implemented GLM-mu rlog workflow and retain the dispersion fit state.
    pub fn rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.rlog_glm_mu_with_fit(counts, design));
        }

        self.rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Apply the implemented GLM-mu rlog workflow with explicit options and retained fit state.
    pub fn rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)
            });
        }

        let fit = self.fit_dispersion_trend_glm_mu(counts, design)?;
        let rlog = fit.regularized_log_transform_with_options(counts, rlog_irls_options)?;
        Ok(RlogGlmMuOutput {
            rlog,
            fit,
            design_mode: RlogDesignMode::DesignAware,
        })
    }

    /// Learn rlog intercepts and immediately run a frozen-intercept rlog reuse pass.
    pub fn frozen_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.frozen_rlog_glm_mu_with_fit(counts, design));
        }

        self.frozen_rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Learn rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.frozen_rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)
            });
        }

        let source = self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)?;
        let frozen_rlog = source.fit.frozen_rlog(
            counts,
            &source.rlog.intercept,
            source.rlog.sample_prior_variance,
        )?;
        Ok(FrozenRlogGlmMuOutput {
            source_rlog: source.rlog,
            frozen_rlog,
            fit: source.fit,
            design_mode: source.design_mode,
        })
    }

    /// Apply the implemented GLM-mu rlog workflow with an intercept-only design.
    ///
    /// This mirrors the implemented part of DESeq2's `blind=TRUE` rlog shape:
    /// dispersion fitting and rlog prior estimation ignore sample groups by
    /// using a one-column all-ones design.
    pub fn blind_rlog_glm_mu(&self, counts: &CountMatrix) -> Result<RlogOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.blind_rlog_glm_mu(counts));
        }

        self.blind_rlog_glm_mu_with_fit(counts)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog with explicit rlog IRLS options.
    pub fn blind_rlog_glm_mu_with_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.blind_rlog_glm_mu_with_options(counts, rlog_irls_options));
        }

        self.blind_rlog_glm_mu_with_fit_and_options(counts, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog and retain the intercept-only dispersion fit state.
    pub fn blind_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.blind_rlog_glm_mu_with_fit(counts));
        }

        self.blind_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Apply blind GLM-mu rlog with explicit options and retained fit state.
    pub fn blind_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.blind_rlog_glm_mu_with_fit_and_options(counts, rlog_irls_options)
            });
        }

        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        let mut output =
            self.rlog_glm_mu_with_fit_and_options(counts, &design, rlog_irls_options)?;
        output.design_mode = RlogDesignMode::Blind;
        Ok(output)
    }

    /// Learn blind rlog intercepts and run a frozen-intercept rlog reuse pass.
    pub fn blind_frozen_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.blind_frozen_rlog_glm_mu_with_fit(counts));
        }

        self.blind_frozen_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Learn blind rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn blind_frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| {
                self.blind_frozen_rlog_glm_mu_with_fit_and_options(counts, rlog_irls_options)
            });
        }

        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        let mut output =
            self.frozen_rlog_glm_mu_with_fit_and_options(counts, &design, rlog_irls_options)?;
        output.design_mode = RlogDesignMode::Blind;
        Ok(output)
    }
}
