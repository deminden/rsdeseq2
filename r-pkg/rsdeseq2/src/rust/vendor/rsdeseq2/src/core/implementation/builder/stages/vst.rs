//! Builder vst stages.
use super::*;

impl DeseqBuilder {
    /// Fit the selected GLM-mu dispersion trend on a fast-VST subset.
    ///
    /// This is a building block for the high-level fast `vst()` workflow:
    /// size factors and normalization factors are derived from the full
    /// dataset, the deterministic fast-VST row subset is selected from the
    /// full-data `baseMean`, and the dispersion trend is fit on the subset
    /// count matrix. The returned subset keeps the original row indices and
    /// aligned normalized counts/factors for inspection. Selection uses
    /// ordinary means of normalized counts, matching DESeq2 1.52.0 even when
    /// the trend fit propagates observation weights.
    pub fn fit_fast_vst_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<(DeseqFit, FastVstSubset), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_fast_vst_dispersion_trend_glm_mu(counts, design, nsub));
        }

        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "fast VST subset size must be positive".to_string(),
            });
        }
        design.validate_full_rank("fast VST GLM-mu dispersion trend")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let selection_means = base_mean(&stages.normalized)?;
        let subset = build_fast_vst_subset(
            counts,
            &stages.normalized,
            &selection_means,
            nsub,
            stages.normalization_factors.as_ref(),
            stages.observation_weights.as_ref(),
        )?;
        let mut subset_builder = self.clone();
        subset_builder.size_factor_options.supplied_size_factors = Some(stages.size_factors);
        subset_builder.normalization_factors = subset.normalization_factors.clone();
        subset_builder.observation_weights = subset.observation_weights.clone();
        let fit = subset_builder.fit_dispersion_trend_glm_mu(&subset.counts, design)?;
        Ok((fit, subset))
    }

    /// Fit the selected GLM-mu dispersion trend on the default fast-VST subset.
    ///
    /// Uses DESeq2's default `nsub=1000`.
    pub fn fit_default_fast_vst_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, FastVstSubset), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_default_fast_vst_dispersion_trend_glm_mu(counts, design));
        }

        self.fit_fast_vst_dispersion_trend_glm_mu(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply a fast-VST transform using a GLM-mu trend fit on the fast subset.
    ///
    /// This mirrors the implemented part of DESeq2's fast `vst()` workflow:
    /// the dispersion trend is estimated on the deterministic subset, then the
    /// fitted trend is applied to the full normalized count matrix. The subset
    /// fit and row-aligned subset bundle are returned for diagnostics.
    pub fn fast_vst_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<FastVstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fast_vst_glm_mu(counts, design, nsub));
        }

        let (subset_fit, subset) =
            self.fit_fast_vst_dispersion_trend_glm_mu(counts, design, nsub)?;
        let trend_fit =
            subset_fit
                .dispersion_trend
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "a fitted fast-VST dispersion trend is required".to_string(),
                })?;
        // The subset supplies only the trend; the transform is still evaluated
        // against the full count matrix and the caller's full-data offsets.
        let full_fit = self.base_fit(
            counts,
            Some(design.clone()),
            BaseFitInput {
                size_factors: subset_fit.size_factors.clone(),
                normalization_factors: self.normalization_factors.clone(),
                observation_weights: None,
                weights_fail: None,
                weights_design_rank: None,
                base_mean: vec![f64::NAN; counts.n_genes()],
                base_var: vec![f64::NAN; counts.n_genes()],
                all_zero: counts.all_zero_flags(),
            },
        );
        let transformed = full_fit.variance_stabilizing_transform_with_trend(counts, trend_fit)?;
        Ok(FastVstGlmMuOutput {
            transformed,
            subset_fit,
            subset,
        })
    }

    /// Apply fast VST using DESeq2's default `nsub=1000`.
    pub fn default_fast_vst_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<FastVstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.default_fast_vst_glm_mu(counts, design));
        }

        self.fast_vst_glm_mu(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply Rust's automatic GLM-mu VST path.
    ///
    /// When at least `nsub` rows have `baseMean > 5`, the dispersion trend is
    /// fit on the deterministic fast-VST subset and applied to the full count
    /// matrix. Otherwise, the selected GLM-mu dispersion trend is fit on the
    /// full count matrix before transforming the full matrix. DESeq2 1.52.0
    /// instead stops when too few rows qualify. Both weighted and unweighted
    /// Rust paths select rows using ordinary normalized-count row means, as
    /// DESeq2 does; the automatic builder helper retains its documented
    /// full-data fallback for too few eligible rows.
    pub fn vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.vst_glm_mu_auto(counts, design, nsub));
        }

        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "automatic VST subset size must be positive".to_string(),
            });
        }
        let stages = self.normalization_stages_for_design(counts, design)?;
        let selection_means = base_mean(&stages.normalized)?;
        let eligible_rows = count_fast_vst_eligible_rows(&selection_means)?;
        // Rust uses the deterministic subset when enough rows qualify and the
        // full dataset otherwise. DESeq2's wrapper errors in the latter case.
        if eligible_rows >= nsub {
            let output = self.fast_vst_glm_mu(counts, design, nsub)?;
            return Ok(VstGlmMuOutput {
                transformed: output.transformed,
                trend_fit: output.subset_fit,
                trend_source: VstTrendSource::FastSubset {
                    nsub,
                    eligible_rows,
                },
                fast_subset: Some(output.subset),
            });
        }

        let trend_fit = self.fit_dispersion_trend_glm_mu(counts, design)?;
        let transformed = trend_fit.variance_stabilizing_transform(counts)?;
        Ok(VstGlmMuOutput {
            transformed,
            trend_fit,
            trend_source: VstTrendSource::FullData {
                nsub,
                eligible_rows,
                reason: VstFullDataReason::InsufficientEligibleRows,
            },
            fast_subset: None,
        })
    }

    /// Apply automatic GLM-mu VST using DESeq2's default `nsub=1000`.
    pub fn default_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.default_vst_glm_mu_auto(counts, design));
        }

        self.vst_glm_mu_auto(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply automatic GLM-mu VST with an intercept-only design.
    ///
    /// This mirrors the implemented part of DESeq2's `blind=TRUE` VST shape:
    /// the transform ignores sample groups by fitting the selected dispersion
    /// trend with a one-column all-ones design, then uses the same automatic
    /// fast-subset/full-data decision as [`Self::vst_glm_mu_auto`].
    pub fn blind_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        nsub: usize,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.blind_vst_glm_mu_auto(counts, nsub));
        }

        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        self.vst_glm_mu_auto(counts, &design, nsub)
    }

    /// Apply blind automatic GLM-mu VST using DESeq2's default `nsub=1000`.
    pub fn default_blind_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.default_blind_vst_glm_mu_auto(counts));
        }

        self.blind_vst_glm_mu_auto(counts, DEFAULT_FAST_VST_NSUB)
    }
}
