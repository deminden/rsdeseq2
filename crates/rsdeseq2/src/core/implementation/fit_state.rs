/// Inspectable fit state for all implemented and future DESeq2 stages.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqFit {
    /// Count matrix summary.
    pub counts_summary: CountsSummary,
    /// Optional design matrix information.
    pub design: Option<DesignMatrix>,
    /// Optional reduced design matrix information for likelihood-ratio tests.
    pub reduced_design: Option<DesignMatrix>,
    /// Optional formula/model-frame metadata used to construct the design.
    ///
    /// This lets already-fitted object routes keep factor levels, references,
    /// and R-cleaned alias behavior without asking callers to repeat the same
    /// metadata for later result or replacement/refit requests.
    pub model_frame: Option<FormulaModelFrame>,
    /// Estimated or supplied size factors.
    pub size_factors: Vec<f64>,
    /// Optional gene/sample normalization factors.
    pub normalization_factors: Option<RowMajorMatrix<f64>>,
    /// Optional DESeq2-style observation weights used by builder stages.
    ///
    /// When a design is available these are row-normalized weights from
    /// `getAndCheckWeights`-style preprocessing. `baseMean` and `baseVar` still
    /// use the raw caller-supplied weights first, matching DESeq2's
    /// `getBaseMeansAndVariances` ordering.
    pub observation_weights: Option<RowMajorMatrix<f64>>,
    /// Rows that failed DESeq2-style observation-weight design checks.
    pub weights_fail: Option<Vec<bool>>,
    /// Rank of the unweighted design during observation-weight preprocessing.
    pub weights_design_rank: Option<usize>,
    /// Per-gene base means.
    pub base_mean: Vec<f64>,
    /// Per-gene sample variance of normalized counts, matching DESeq2 `baseVar`.
    pub base_var: Vec<f64>,
    /// Per-gene all-zero flags, matching DESeq2 `allZero`.
    ///
    /// For design-aware observation-weight stages, rows with `weightsFail` are
    /// also marked true, matching DESeq2's downstream skip behavior.
    pub all_zero: Vec<bool>,
    /// Gene-wise dispersion estimates.
    pub disp_gene_est: Option<Vec<f64>>,
    /// Gene-wise dispersion iteration counts, matching DESeq2 `dispGeneIter`.
    pub disp_gene_iter: Option<Vec<usize>>,
    /// Fitted dispersion trend values.
    pub disp_fit: Option<Vec<f64>>,
    /// Fitted dispersion trend object used by implemented VST dispatch.
    pub dispersion_trend: Option<DispersionTrendFit>,
    /// MAP dispersion estimates before outlier replacement.
    pub disp_map: Option<Vec<f64>>,
    /// Final dispersion estimates.
    pub dispersion: Option<Vec<f64>>,
    /// MAP dispersion iteration counts.
    pub disp_iter: Option<Vec<usize>>,
    /// MAP dispersion outlier flags.
    pub disp_outlier: Option<Vec<bool>>,
    /// Dispersion prior variance.
    pub disp_prior_var: Option<f64>,
    /// Robust variance of log-dispersion estimates used for MAP outlier checks.
    pub var_log_disp_estimates: Option<f64>,
    /// Dispersion convergence flags.
    pub dispersion_converged: Option<Vec<bool>>,
    /// GLM beta estimates.
    pub beta: Option<RowMajorMatrix<f64>>,
    /// GLM beta standard errors.
    pub beta_se: Option<RowMajorMatrix<f64>>,
    /// Fallback optimizer starting beta values on log2 scale.
    pub beta_optim_start: Option<RowMajorMatrix<f64>>,
    /// Per-gene GLM beta covariance matrices on log2 scale.
    ///
    /// Stored as genes x `(n_coefficients * n_coefficients)`, with each gene
    /// row containing a row-major coefficient covariance matrix.
    pub beta_covariance: Option<RowMajorMatrix<f64>>,
    /// GLM beta convergence flags.
    pub beta_converged: Option<Vec<bool>>,
    /// GLM beta iteration counts.
    pub beta_iter: Option<Vec<usize>>,
    /// Rust fallback-optimizer iterations for rows routed after IRLS.
    pub beta_optim_iter: Option<Vec<f64>>,
    /// Rust fallback-optimizer objective at the starting parameter vector.
    pub beta_optim_start_objective: Option<Vec<f64>>,
    /// Final Rust fallback-optimizer objective for rows routed after IRLS.
    pub beta_optim_objective: Option<Vec<f64>>,
    /// Projected gradient norm at the final Rust fallback-optimizer parameters.
    pub beta_optim_gradient_norm: Option<Vec<f64>>,
    /// Per-gene fitted-model log likelihoods from the full GLM.
    pub log_like: Option<Vec<f64>>,
    /// Per-gene full-model deviance, matching DESeq2's `-2 * logLike` field.
    pub full_deviance: Option<Vec<f64>>,
    /// Per-gene reduced-model log likelihoods for LRT pipelines.
    pub reduced_log_like: Option<Vec<f64>>,
    /// Per-gene reduced-model beta convergence flags for LRT pipelines.
    pub reduced_beta_converged: Option<Vec<bool>>,
    /// Per-gene reduced-model beta iteration counts for LRT pipelines.
    pub reduced_beta_iter: Option<Vec<usize>>,
    /// Reduced-model fitted mean matrix for LRT pipelines.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub reduced_mu: Option<RowMajorMatrix<f64>>,
    /// Reduced-model hat diagonal matrix for LRT pipelines.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub reduced_hat_diagonal: Option<RowMajorMatrix<f64>>,
    /// Fitted mean matrix.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub mu: Option<RowMajorMatrix<f64>>,
    /// Cook's distance matrix.
    pub cooks: Option<RowMajorMatrix<f64>>,
    /// Per-gene maximum Cook's distance over eligible samples.
    pub max_cooks: Option<Vec<Option<f64>>>,
    /// Hat diagonal matrix.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub hat_diagonal: Option<RowMajorMatrix<f64>>,
    /// Wald output.
    pub wald: Option<WaldOutput>,
    /// LRT output.
    pub lrt: Option<LrtOutput>,
}

/// Builder for implemented and future DESeq2 workflow stages.
#[derive(Clone, Debug)]
pub struct DeseqBuilder {
    fit_type: FitType,
    test: TestType,
    size_factor_options: SizeFactorOptions,
    normalization_factors: Option<RowMajorMatrix<f64>>,
    observation_weights: Option<RowMajorMatrix<f64>>,
    observation_weight_options: ObservationWeightOptions,
    execution_mode: ExecutionMode,
    threads: Option<usize>,
    reduced_design: Option<DesignMatrix>,
    model_frame: Option<FormulaModelFrame>,
    irls_options: IrlsOptions,
    gene_wise_dispersion_options: GeneWiseDispersionOptions,
    wald_test_options: WaldTestOptions,
    cooks_cutoff: CooksCutoff,
    independent_filtering_options: IndependentFilteringOptions,
}

impl DeseqFit {
    /// Reconstruct DESeq2-style normalized counts for the count matrix used to create this fit.
    pub fn normalized_counts(
        &self,
        counts: &CountMatrix,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "normalized counts")?;
        match &self.normalization_factors {
            Some(factors) => normalized_counts_with_factors(counts, factors),
            None => normalized_counts(counts, &self.size_factors),
        }
    }

    /// Apply DESeq2's `normTransform` to the count matrix used to create this fit.
    pub fn norm_transform(&self, counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
        let normalized = self.normalized_counts(counts)?;
        norm_transform(&normalized)
    }

    /// Apply VST using this fit's stored dispersion trend.
    ///
    /// This is a lower-level analogue of DESeq2's
    /// `getVarianceStabilizedData`: the fit must already include a fitted
    /// dispersion trend from one of the implemented trend stages.
    pub fn variance_stabilizing_transform(
        &self,
        counts: &CountMatrix,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "variance-stabilizing transform")?;
        let trend_fit =
            self.dispersion_trend
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "a fitted dispersion trend is required before VST".to_string(),
                })?;
        self.variance_stabilizing_transform_with_trend(counts, trend_fit)
    }

    /// Apply VST using a caller-supplied fitted dispersion trend.
    ///
    /// This supports fast VST, where normalization and
    /// full-count reconstruction come from the full fit but the trend may have
    /// been estimated on a deterministic subset.
    pub fn variance_stabilizing_transform_with_trend(
        &self,
        counts: &CountMatrix,
        trend_fit: &DispersionTrendFit,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "variance-stabilizing transform")?;
        let normalized = self.normalized_counts(counts)?;
        match &self.normalization_factors {
            Some(factors) => {
                vst_with_dispersion_trend_and_normalization_factors(&normalized, trend_fit, factors)
            }
            None => vst_with_dispersion_trend_and_size_factors(
                &normalized,
                trend_fit,
                &self.size_factors,
            ),
        }
    }

    /// Short alias for [`DeseqFit::variance_stabilizing_transform`].
    pub fn vst(&self, counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
        self.variance_stabilizing_transform(counts)
    }

    /// Formula/model-frame metadata stored with this fit, if available.
    pub fn current_model_frame(&self) -> Option<&FormulaModelFrame> {
        self.model_frame.as_ref()
    }

    /// Assemble `results(contrast=...)` rows from this already-fitted object.
    ///
    /// The stored test state chooses Wald or LRT behavior. Character triplet
    /// contrasts use retained formula/model-frame metadata when available.
    pub fn results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        match (self.wald.is_some(), self.lrt.is_some()) {
            (true, false) => self.wald_results_contrast_request(counts, contrast, gene_names),
            (false, true) => self.lrt_results_contrast_request(counts, contrast, gene_names),
            (true, true) => Err(DeseqError::InvalidOptions {
                reason: "fitted-object contrast results require either Wald or LRT output, not both"
                    .to_string(),
            }),
            (false, false) => Err(DeseqError::InvalidOptions {
                reason: "fitted-object contrast results require stored Wald or LRT output"
                    .to_string(),
            }),
        }
    }

    /// Assemble Wald `results(contrast=...)` rows from this already-fitted object.
    ///
    /// Character triplet contrasts use stored formula/model-frame metadata for
    /// factor levels and references. List and numeric contrasts use the design
    /// coefficient names carried by the fitted GLM state.
    pub fn wald_results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "Wald results contrast request")?;
        let fit = self.as_nbinom_glm_fit()?;
        let numeric_contrast = match contrast {
            ResultsContrast::Character { .. } => {
                let model_frame =
                    self.model_frame
                        .as_ref()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "character results contrast requires stored formula model-frame metadata; use wald_results_factor_level_contrast when supplying explicit sample levels".to_string(),
                        })?;
                let factor_contrast = factor_level_contrast_from_model_frame(contrast, model_frame)?
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "character results contrast does not match stored formula model-frame metadata".to_string(),
                    })?;
                return self.wald_results_factor_level_contrast(counts, factor_contrast, gene_names);
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                contrast.as_contrast_spec()
            }
        };
        let numeric = resolve_contrast(&fit.model_matrix, &numeric_contrast)?;
        let mut output =
            wald_test_contrast_with_options(&fit, &numeric, &WaldTestOptions::default())?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, &fit.model_matrix, &numeric)?;
        apply_contrast_all_zero_to_wald_contrast(
            &mut output,
            &contrast_all_zero,
            &self.all_zero,
        )?;
        let mut results = build_wald_contrast_results(
            &self.base_mean,
            &fit,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        results.apply_wald_test_options(&WaldTestOptions::default());
        self.attach_fit_diagnostics_to_results(&mut results)?;
        results.set_resolved_contrast_metadata(
            numeric_contrast.result_name(),
            numeric_contrast.comparison(),
            &numeric,
        );
        Ok(results)
    }

    /// Assemble Wald `results(contrast=...)` rows for an explicit factor-level contrast.
    ///
    /// This is the fitted-object route for callers that already know the
    /// per-sample factor levels. Formula-built fits can instead use
    /// [`DeseqFit::wald_results_contrast_request`] with stored model-frame
    /// metadata.
    pub fn wald_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        contrast: FactorLevelContrast<'_>,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "Wald factor-level results contrast")?;
        let fit = self.as_nbinom_glm_fit()?;
        let contrast_spec = factor_level_contrast_spec(contrast);
        let numeric = resolve_contrast(&fit.model_matrix, &contrast_spec)?;
        let contrast_all_zero = contrast_all_zero_factor_levels(
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let mut output =
            wald_test_contrast_with_options(&fit, &numeric, &WaldTestOptions::default())?;
        apply_contrast_all_zero_to_wald_contrast(
            &mut output,
            &contrast_all_zero,
            &self.all_zero,
        )?;
        let mut results = build_wald_contrast_results(
            &self.base_mean,
            &fit,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        results.apply_wald_test_options(&WaldTestOptions::default());
        self.attach_fit_diagnostics_to_results(&mut results)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric);
        Ok(results)
    }

    /// Assemble LRT `results(contrast=...)` rows from this already-fitted object.
    ///
    /// The LRT statistic and p-value come from the stored full-vs-reduced test;
    /// the requested contrast only changes the reported effect-size columns.
    pub fn lrt_results_contrast_request(
        &self,
        counts: &CountMatrix,
        contrast: &ResultsContrast,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "LRT results contrast request")?;
        let fit = self.as_nbinom_glm_fit()?;
        let lrt = self.lrt.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored LRT output is required for LRT contrast results".to_string(),
        })?;
        let numeric_contrast = match contrast {
            ResultsContrast::Character { .. } => {
                let model_frame =
                    self.model_frame
                        .as_ref()
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "character results contrast requires stored formula model-frame metadata; use lrt_results_factor_level_contrast when supplying explicit sample levels".to_string(),
                        })?;
                let factor_contrast = factor_level_contrast_from_model_frame(contrast, model_frame)?
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "character results contrast does not match stored formula model-frame metadata".to_string(),
                    })?;
                return self.lrt_results_factor_level_contrast(counts, factor_contrast, gene_names);
            }
            ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                contrast.as_contrast_spec()
            }
        };
        let numeric = resolve_contrast(&fit.model_matrix, &numeric_contrast)?;
        let output = wald_test_contrast_with_options(&fit, &numeric, &WaldTestOptions::default())?;
        let mut results = build_lrt_contrast_results(
            &self.base_mean,
            &fit,
            lrt,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        self.attach_fit_diagnostics_to_results(&mut results)?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, &fit.model_matrix, &numeric)?;
        apply_contrast_all_zero_to_lrt_results(&mut results, &contrast_all_zero, &self.all_zero)?;
        results.set_resolved_contrast_metadata(
            numeric_contrast.result_name(),
            numeric_contrast.comparison(),
            &numeric,
        );
        Ok(results)
    }

    /// Assemble LRT `results(contrast=...)` rows for an explicit factor-level contrast.
    ///
    /// The stored full-vs-reduced LRT supplies the statistic and p-value; the
    /// explicit factor contrast supplies reported effect-size columns and
    /// character-contrast all-zero handling.
    pub fn lrt_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        contrast: FactorLevelContrast<'_>,
        gene_names: Option<&[String]>,
    ) -> Result<DeseqResults, DeseqError> {
        validate_fit_counts_shape(self, counts, "LRT factor-level results contrast")?;
        let fit = self.as_nbinom_glm_fit()?;
        let lrt = self.lrt.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored LRT output is required for LRT contrast results".to_string(),
        })?;
        let contrast_spec = factor_level_contrast_spec(contrast);
        let numeric = resolve_contrast(&fit.model_matrix, &contrast_spec)?;
        let contrast_all_zero = contrast_all_zero_factor_levels(
            counts,
            contrast.sample_levels,
            contrast.numerator,
            contrast.denominator,
        )?;
        let output = wald_test_contrast_with_options(&fit, &numeric, &WaldTestOptions::default())?;
        let mut results = build_lrt_contrast_results(
            &self.base_mean,
            &fit,
            lrt,
            &output,
            gene_names,
            self.dispersion.as_deref(),
        )?;
        self.attach_fit_diagnostics_to_results(&mut results)?;
        apply_contrast_all_zero_to_lrt_results(&mut results, &contrast_all_zero, &self.all_zero)?;
        let (result_name, comparison) = factor_level_result_metadata(contrast);
        results.set_resolved_contrast_metadata(result_name, comparison, &numeric);
        Ok(results)
    }

    /// Apply the implemented regularized-log sample-effect transform.
    ///
    /// The fit must already contain fitted trend dispersions (`dispFit`) and
    /// final gene-wise dispersions (`dispersion`). The rlog sample-effect prior
    /// is estimated from this fit's normalized counts, `baseMean`, and
    /// `dispFit`, then the sample-effect ridge GLM is fit with this fit's
    /// size factors or normalization factors.
    pub fn regularized_log_transform(
        &self,
        counts: &CountMatrix,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_options(counts, IrlsOptions::default())
    }

    /// Apply the implemented regularized-log transform with explicit IRLS options.
    pub fn regularized_log_transform_with_options(
        &self,
        counts: &CountMatrix,
        options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        validate_fit_counts_shape(self, counts, "regularized-log transform")?;
        if self.all_zero.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog allZero rows",
                counts.n_genes(),
                self.all_zero.len(),
            ));
        }
        let disp_fit = self
            .disp_fit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted dispersion trend values are required before rlog".to_string(),
            })?;
        let dispersions =
            self.dispersion
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "final dispersions are required before rlog".to_string(),
                })?;
        let nonzero_rows = nonzero_gene_indices(&self.all_zero);
        if nonzero_rows.len() != counts.n_genes() {
            if nonzero_rows.is_empty() {
                let transformed =
                    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 0.0)?;
                let offset_mode = match &self.normalization_factors {
                    Some(_) => RlogOffsetMode::NormalizationFactors,
                    None => RlogOffsetMode::SizeFactors,
                };
                return Ok(RlogOutput {
                    transformed,
                    intercept: vec![0.0; counts.n_genes()],
                    sample_prior_variance: 0.0,
                    offset_mode,
                });
            }
            let compact_counts = compact_counts(counts, &nonzero_rows)?;
            let compact_base_mean = compact_f64_values(&self.base_mean, &nonzero_rows)?;
            let compact_disp_fit = compact_f64_values(disp_fit, &nonzero_rows)?;
            let compact_dispersions = compact_f64_values(dispersions, &nonzero_rows)?;
            let compact_output = match &self.normalization_factors {
                Some(factors) => {
                    let compact_factors = compact_matrix_rows(factors, &nonzero_rows)?;
                    rlog_with_estimated_prior_and_normalization_factors(
                        &compact_counts,
                        &compact_factors,
                        &compact_base_mean,
                        &compact_disp_fit,
                        &compact_dispersions,
                        options,
                    )?
                }
                None => rlog_with_estimated_prior_and_size_factors(
                    &compact_counts,
                    &self.size_factors,
                    &compact_base_mean,
                    &compact_disp_fit,
                    &compact_dispersions,
                    options,
                )?,
            };
            return expand_rlog_output_with_all_zero_rows(
                compact_output,
                &self.all_zero,
                counts.n_samples(),
            );
        }
        match &self.normalization_factors {
            Some(factors) => rlog_with_estimated_prior_and_normalization_factors(
                counts,
                factors,
                &self.base_mean,
                disp_fit,
                dispersions,
                options,
            ),
            None => rlog_with_estimated_prior_and_size_factors(
                counts,
                &self.size_factors,
                &self.base_mean,
                disp_fit,
                dispersions,
                options,
            ),
        }
    }

    /// Short alias for [`DeseqFit::regularized_log_transform`].
    pub fn rlog(&self, counts: &CountMatrix) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform(counts)
    }

    /// Apply rlog with supplied frozen intercepts and sample-effect prior variance.
    ///
    /// This uses the fit state's final dispersions and size-factor or
    /// normalization-factor offsets, while treating `frozen_intercept` as the
    /// fixed log2 intercept vector. It is the fit-level building block for
    /// frozen rlog reuse.
    pub fn regularized_log_transform_with_frozen_intercept(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_frozen_intercept_and_options(
            counts,
            frozen_intercept,
            sample_prior_variance,
            IrlsOptions::default(),
        )
    }

    /// Apply frozen-intercept rlog with explicit IRLS options.
    pub fn regularized_log_transform_with_frozen_intercept_and_options(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
        options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        validate_fit_counts_shape(self, counts, "frozen regularized-log transform")?;
        if frozen_intercept.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog frozen intercept rows",
                counts.n_genes(),
                frozen_intercept.len(),
            ));
        }
        if self.all_zero.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog allZero rows",
                counts.n_genes(),
                self.all_zero.len(),
            ));
        }
        let dispersions =
            self.dispersion
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "final dispersions are required before frozen rlog".to_string(),
                })?;
        let nonzero_rows = nonzero_gene_indices(&self.all_zero);
        if nonzero_rows.len() != counts.n_genes() {
            if nonzero_rows.is_empty() {
                let transformed =
                    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 0.0)?;
                let offset_mode = match &self.normalization_factors {
                    Some(_) => RlogOffsetMode::NormalizationFactors,
                    None => RlogOffsetMode::SizeFactors,
                };
                return Ok(RlogOutput {
                    transformed,
                    intercept: frozen_intercept.to_vec(),
                    sample_prior_variance,
                    offset_mode,
                });
            }
            let compact_counts = compact_counts(counts, &nonzero_rows)?;
            let compact_dispersions = compact_f64_values(dispersions, &nonzero_rows)?;
            let compact_intercept = compact_f64_values(frozen_intercept, &nonzero_rows)?;
            let compact_output = match &self.normalization_factors {
                Some(factors) => {
                    let compact_factors = compact_matrix_rows(factors, &nonzero_rows)?;
                    rlog_frozen_with_normalization_factors(
                        &compact_counts,
                        &compact_factors,
                        &compact_dispersions,
                        sample_prior_variance,
                        &compact_intercept,
                        options,
                    )?
                }
                None => rlog_frozen_with_size_factors(
                    &compact_counts,
                    &self.size_factors,
                    &compact_dispersions,
                    sample_prior_variance,
                    &compact_intercept,
                    options,
                )?,
            };
            return expand_frozen_rlog_output_with_all_zero_rows(
                compact_output,
                &self.all_zero,
                counts.n_samples(),
                frozen_intercept,
            );
        }
        match &self.normalization_factors {
            Some(factors) => rlog_frozen_with_normalization_factors(
                counts,
                factors,
                dispersions,
                sample_prior_variance,
                frozen_intercept,
                options,
            ),
            None => rlog_frozen_with_size_factors(
                counts,
                &self.size_factors,
                dispersions,
                sample_prior_variance,
                frozen_intercept,
                options,
            ),
        }
    }

    /// Short alias for [`DeseqFit::regularized_log_transform_with_frozen_intercept`].
    pub fn frozen_rlog(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_frozen_intercept(
            counts,
            frozen_intercept,
            sample_prior_variance,
        )
    }

    /// Number of rows eligible under this fit's fast-VST selection statistic.
///
    /// This uses the stored `baseMean` vector and the `baseMean > 5` rule. For
    /// weighted fits the stored values are weighted, unlike DESeq2 1.52.0's
    /// ordinary normalized-count row means.
    pub fn fast_vst_eligible_count(&self) -> Result<usize, DeseqError> {
        count_fast_vst_eligible_rows(&self.base_mean)
    }

    /// Build the row-aligned subset used by this fit's fast-VST trend fit.
    ///
    /// The returned bundle includes raw counts, normalized counts, optional
    /// normalization factors, optional observation weights, and original row
    /// indices, all selected from this fit's `baseMean` vector.
    pub fn fast_vst_subset(
        &self,
        counts: &CountMatrix,
        nsub: usize,
    ) -> Result<FastVstSubset, DeseqError> {
        validate_fit_counts_shape(self, counts, "fast VST subset")?;
        let normalized = self.normalized_counts(counts)?;
        build_fast_vst_subset(
            counts,
            &normalized,
            &self.base_mean,
            nsub,
            self.normalization_factors.as_ref(),
            self.observation_weights.as_ref(),
        )
    }

    fn as_nbinom_glm_fit(&self) -> Result<NbinomGlmFit, DeseqError> {
        let design = self.design.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored design is required for fitted-object results".to_string(),
        })?;
        let beta = self.beta.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored beta matrix is required for fitted-object results".to_string(),
        })?;
        let beta_se = self
            .beta_se
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "stored beta standard-error matrix is required for fitted-object results"
                    .to_string(),
            })?;
        let beta_optim_start =
            self.beta_optim_start
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta optimizer starts are required for fitted-object results"
                        .to_string(),
                })?;
        let beta_converged =
            self.beta_converged
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta convergence flags are required for fitted-object results"
                        .to_string(),
                })?;
        let beta_iter = self.beta_iter.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored beta iteration counts are required for fitted-object results".to_string(),
        })?;
        let beta_optim_iter =
            self.beta_optim_iter
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta optimizer iterations are required for fitted-object results"
                        .to_string(),
                })?;
        let beta_optim_start_objective =
            self.beta_optim_start_objective
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta optimizer start objectives are required for fitted-object results"
                        .to_string(),
                })?;
        let beta_optim_objective =
            self.beta_optim_objective
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta optimizer objectives are required for fitted-object results"
                        .to_string(),
                })?;
        let beta_optim_gradient_norm =
            self.beta_optim_gradient_norm
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored beta optimizer gradient norms are required for fitted-object results"
                        .to_string(),
                })?;
        let log_like = self.log_like.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored log likelihoods are required for fitted-object results".to_string(),
        })?;
        let mu = self.mu.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "stored fitted means are required for fitted-object results".to_string(),
        })?;
        let hat_diagonal =
            self.hat_diagonal
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "stored hat diagonals are required for fitted-object results".to_string(),
                })?;
        Ok(NbinomGlmFit {
            log_like: log_like.clone(),
            beta_converged: beta_converged.clone(),
            beta: beta.clone(),
            beta_se: beta_se.clone(),
            beta_optim_start: beta_optim_start.clone(),
            beta_covariance: self.beta_covariance.clone(),
            mu: mu.clone(),
            beta_iter: beta_iter.clone(),
            beta_optim_iter: beta_optim_iter.clone(),
            beta_optim_start_objective: beta_optim_start_objective.clone(),
            beta_optim_objective: beta_optim_objective.clone(),
            beta_optim_gradient_norm: beta_optim_gradient_norm.clone(),
            model_matrix: design.clone(),
            n_terms: design.n_coefficients(),
            hat_diagonal: hat_diagonal.clone(),
        })
    }

    fn attach_fit_diagnostics_to_results(
        &self,
        results: &mut DeseqResults,
    ) -> Result<(), DeseqError> {
        if results.rows.len() != self.all_zero.len() {
            return Err(invalid_dimensions(
                "fitted-object result rows",
                self.all_zero.len(),
                results.rows.len(),
            ));
        }
        if let Some(max_cooks) = &self.max_cooks {
            if max_cooks.len() != results.rows.len() {
                return Err(invalid_dimensions(
                    "fitted-object maxCooks rows",
                    results.rows.len(),
                    max_cooks.len(),
                ));
            }
            for (row, max_cook) in results.rows.iter_mut().zip(max_cooks.iter().copied()) {
                row.max_cooks = max_cook;
            }
        }
        for (row, all_zero) in results.rows.iter_mut().zip(&self.all_zero) {
            if *all_zero {
                row.converged = None;
                row.max_cooks = None;
            }
        }
        Ok(())
    }
}

fn factor_level_contrast_spec(contrast: FactorLevelContrast<'_>) -> ContrastSpec {
    match contrast.reference {
        Some(reference) => ContrastSpec::factor_level_with_reference(
            contrast.factor,
            contrast.numerator,
            contrast.denominator,
            reference,
        ),
        None => ContrastSpec::factor_level(contrast.factor, contrast.numerator, contrast.denominator),
    }
}

fn validate_fit_counts_shape(
    fit: &DeseqFit,
    counts: &CountMatrix,
    context: &str,
) -> Result<(), DeseqError> {
    let expected = fit.counts_summary.n_genes * fit.counts_summary.n_samples;
    let actual = counts.n_genes() * counts.n_samples();
    if fit.counts_summary.n_genes != counts.n_genes()
        || fit.counts_summary.n_samples != counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected,
            actual,
        });
    }
    Ok(())
}
