impl DeseqBuilder {
    /// Run only the implemented initial normalization stages.
    pub fn fit_size_factors_and_base_means(
        &self,
        counts: &CountMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let stages = self.normalization_stages(counts)?;
        Ok(Self::base_fit(counts, None, stages.into_base_fit_input()))
    }

    /// Run initial normalization stages with design-aware observation-weight checks.
    ///
    /// This is useful for parity checks against DESeq2's early metadata when a
    /// `weights` assay is present: raw weights are used for `baseMean` and
    /// `baseVar`, weights are row-normalized for fitting checks, design/rank
    /// failures are recorded in `weights_fail`, and failed rows are marked in
    /// `all_zero` for downstream skipping.
    pub fn fit_size_factors_and_base_means_with_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let stages = self.normalization_stages_for_design(counts, design)?;
        Ok(Self::base_fit(
            counts,
            Some(design.clone()),
            stages.into_base_fit_input(),
        ))
    }

    /// Run the current linear-mu gene-wise dispersion estimator.
    ///
    /// This is a narrow Phase 3 stage for designs where DESeq2's
    /// `linearMu=TRUE` branch is appropriate. It estimates size factors,
    /// normalized counts, base row metadata, linear fitted means, and
    /// fixed-mean gene-wise dispersions. Cox-Reid correction, iterative GLM
    /// mean refits, trend fitting, and MAP shrinkage remain future stages.
    pub fn fit_gene_wise_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
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
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
    }

    /// Run the current GLM-mu gene-wise dispersion estimator.
    ///
    /// This is the first non-`linearMu` foundation for
    /// `estimateDispersionsGeneEst`: rough/moments starts are followed by
    /// fixed-dispersion NB GLM mean fitting and fixed-mean dispersion
    /// optimization. Builder-supplied observation weights are preprocessed in
    /// the same design-aware stage used by fixed-dispersion Wald/LRT paths and
    /// then passed into both the fixed-dispersion mean fit and fixed-mean
    /// dispersion objective.
    pub fn fit_gene_wise_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        design.validate_full_rank("GLM-mu dispersion")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let dispersion = estimate_gene_wise_dispersions_glm_mu(
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
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
    }

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
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the mean trend.
    pub fn fit_mean_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    fn attach_parametric_dispersion_trend(
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

    fn attach_mean_dispersion_trend(&self, mut fit: DeseqFit) -> Result<DeseqFit, DeseqError> {
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

    fn attach_local_dispersion_trend(&self, mut fit: DeseqFit) -> Result<DeseqFit, DeseqError> {
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

    fn attach_existing_dispersion_trend(
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
        match self.fit_type {
            FitType::Parametric => self.fit_parametric_dispersion_trend_glm_mu(counts, design),
            FitType::Mean => self.fit_mean_dispersion_trend_glm_mu(counts, design),
            FitType::Local => self.fit_local_dispersion_trend_glm_mu(counts, design),
            FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
                feature: "GLM-mu glmGamPoi dispersion trend fitting".to_string(),
            }),
        }
    }

    /// Fit the selected GLM-mu dispersion trend on DESeq2's fast-VST subset.
    ///
    /// This is a building block for the high-level fast `vst()` workflow:
    /// size factors and normalization factors are derived from the full
    /// dataset, the deterministic fast-VST row subset is selected from the
    /// full-data `baseMean`, and the dispersion trend is fit on the subset
    /// count matrix. The returned subset keeps the original row indices and
    /// aligned normalized counts/factors for inspection.
    pub fn fit_fast_vst_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<(DeseqFit, FastVstSubset), DeseqError> {
        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "fast VST subset size must be positive".to_string(),
            });
        }
        design.validate_full_rank("fast VST GLM-mu dispersion trend")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let subset = build_fast_vst_subset(
            counts,
            &stages.normalized,
            &stages.base_mean,
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
        let full_fit = Self::base_fit(
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
        self.fast_vst_glm_mu(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply the implemented GLM-mu VST path with DESeq2-like fast-path selection.
    ///
    /// When at least `nsub` rows have `baseMean > 5`, the dispersion trend is
    /// fit on the deterministic fast-VST subset and applied to the full count
    /// matrix. Otherwise, the selected GLM-mu dispersion trend is fit on the
    /// full count matrix before transforming the full matrix.
    pub fn vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "automatic VST subset size must be positive".to_string(),
            });
        }
        let base_fit = self.fit_size_factors_and_base_means_with_design(counts, design)?;
        let eligible_rows = base_fit.fast_vst_eligible_count()?;
        // Match DESeq2's fast path gate: enough high-baseMean rows use the
        // deterministic subset, otherwise the trend is fit on all rows.
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
        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        self.vst_glm_mu_auto(counts, &design, nsub)
    }

    /// Apply blind automatic GLM-mu VST using DESeq2's default `nsub=1000`.
    pub fn default_blind_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        self.blind_vst_glm_mu_auto(counts, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply the implemented GLM-mu rlog workflow in one builder call.
    ///
    /// This fits GLM-mu gene-wise, trend, prior, and MAP dispersion stages
    /// using the builder's current options, then applies the fit-state rlog
    /// transform with default IRLS options.
    pub fn rlog_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogOutput, DeseqError> {
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
        self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply the implemented GLM-mu rlog workflow and retain the dispersion fit state.
    pub fn rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        self.rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Apply the implemented GLM-mu rlog workflow with explicit options and retained fit state.
    pub fn rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
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
        self.frozen_rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Learn rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
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
        self.blind_rlog_glm_mu_with_fit(counts)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog with explicit rlog IRLS options.
    pub fn blind_rlog_glm_mu_with_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        self.blind_rlog_glm_mu_with_fit_and_options(counts, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog and retain the intercept-only dispersion fit state.
    pub fn blind_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        self.blind_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Apply blind GLM-mu rlog with explicit options and retained fit state.
    pub fn blind_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
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
        self.blind_frozen_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Learn blind rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn blind_frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        let mut output =
            self.frozen_rlog_glm_mu_with_fit_and_options(counts, &design, rlog_irls_options)?;
        output.design_mode = RlogDesignMode::Blind;
        Ok(output)
    }

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
        let fit = self.fit_parametric_dispersion_trend_glm_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    fn fit_map_dispersions_glm_mu_with_dispersion_function(
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

    fn attach_map_dispersions(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        self.attach_map_dispersions_with_prior_values(counts, design, fit, None, None)
    }

    fn attach_map_dispersions_with_prior_values(
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
