impl DeseqBuilder {
    fn normalization_stages(
        &self,
        counts: &CountMatrix,
    ) -> Result<NormalizationStages, DeseqError> {
        self.normalization_stages_inner(counts, None)
    }

    fn normalization_stages_for_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<NormalizationStages, DeseqError> {
        self.normalization_stages_inner(counts, Some(design))
    }

    fn normalization_stages_inner(
        &self,
        counts: &CountMatrix,
        design: Option<&DesignMatrix>,
    ) -> Result<NormalizationStages, DeseqError> {
        // DESeq2 gives assay-shaped normalization factors precedence over
        // size-factor estimation; supplied size factors are kept only for APIs
        // that still need sample-scale metadata.
        let (size_factors, normalization_factors, normalized) = match &self.normalization_factors {
            Some(factors) => {
                validate_normalization_factors(counts, factors)?;
                let size_factors = match &self.size_factor_options.supplied_size_factors {
                    Some(size_factors) => {
                        normalized_counts(counts, size_factors)?;
                        size_factors.clone()
                    }
                    None => vec![1.0; counts.n_samples()],
                };
                let factors = factors.clone();
                let normalized = normalized_counts_with_factors(counts, &factors)?;
                (size_factors, Some(factors), normalized)
            }
            None => {
                let control_gene_indices = self
                    .size_factor_options
                    .control_genes
                    .as_ref()
                    .map(|control_genes| control_genes.to_indices(counts.n_genes()))
                    .transpose()?;
                let size_factors = match &self.size_factor_options.supplied_size_factors {
                    Some(size_factors) => size_factors.clone(),
                    None => estimate_size_factors_with_options(
                        counts,
                        self.size_factor_options.method,
                        self.size_factor_options.geo_means.as_deref(),
                        control_gene_indices.as_deref(),
                    )?,
                };
                let normalized = normalized_counts(counts, &size_factors)?;
                (size_factors, None, normalized)
            }
        };
        let weighted_metadata = self.weighted_base_metadata(counts, design, &normalized)?;
        let raw_all_zero = counts.all_zero_flags();
        let all_zero = match &weighted_metadata.weights_fail {
            Some(weights_fail) => raw_all_zero
                .iter()
                .copied()
                .zip(weights_fail.iter().copied())
                .map(|(all_zero, weights_fail)| all_zero || weights_fail)
                .collect(),
            None => raw_all_zero,
        };
        Ok(NormalizationStages {
            size_factors,
            base_mean: weighted_metadata.base_mean,
            base_var: weighted_metadata.base_var,
            all_zero,
            normalized,
            normalization_factors,
            observation_weights: weighted_metadata.observation_weights,
            weights_fail: weighted_metadata.weights_fail,
            weights_design_rank: weighted_metadata.weights_design_rank,
        })
    }

    fn weighted_base_metadata(
        &self,
        counts: &CountMatrix,
        design: Option<&DesignMatrix>,
        normalized: &RowMajorMatrix<f64>,
    ) -> Result<WeightedBaseMetadata, DeseqError> {
        let Some(weights) = &self.observation_weights else {
            return Ok(WeightedBaseMetadata {
                base_mean: base_mean(normalized)?,
                base_var: base_variance(normalized)?,
                observation_weights: None,
                weights_fail: None,
                weights_design_rank: None,
            });
        };
        validate_observation_weights_for_counts(counts, weights)?;
        match design {
            Some(design) => {
                // Design-aware paths use checked, row-normalized weights for
                // fitting, but baseMean/baseVar intentionally use raw weights
                // to match DESeq2's metadata ordering.
                let checked = preprocess_observation_weights_with_options(
                    weights,
                    design,
                    self.observation_weight_options,
                )?;
                let base_mean = base_mean_with_weights(normalized, weights)?;
                let base_var = base_variance_with_weights(normalized, weights)?;
                Ok(WeightedBaseMetadata {
                    base_mean,
                    base_var,
                    observation_weights: Some(checked.weights),
                    weights_fail: Some(checked.weights_fail),
                    weights_design_rank: Some(checked.design_rank),
                })
            }
            None => Ok(WeightedBaseMetadata {
                base_mean: base_mean_with_weights(normalized, weights)?,
                base_var: base_variance_with_weights(normalized, weights)?,
                observation_weights: Some(weights.clone()),
                weights_fail: None,
                weights_design_rank: None,
            }),
        }
    }

    fn ensure_no_observation_weights(&self, feature: &str) -> Result<(), DeseqError> {
        if self.observation_weights.is_some() {
            return Err(DeseqError::UnsupportedFeature {
                feature: format!("{feature} with observation weights"),
            });
        }
        Ok(())
    }

    fn fixed_dispersion_wald_components(
        &self,
        input: WaldPipelineInput<'_>,
    ) -> Result<WaldPipelineOutput, DeseqError> {
        if input.coefficient >= input.design.n_coefficients() {
            return Err(DeseqError::InvalidDimensions {
                context: "pipeline Wald coefficient index".to_string(),
                expected: input.design.n_coefficients().saturating_sub(1),
                actual: input.coefficient,
            });
        }
        let FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        } = self.fixed_dispersion_glm_components(FixedDispersionGlmInput {
            counts: input.counts,
            design: input.design,
            size_factors: input.size_factors,
            normalization_factors: input.normalization_factors,
            observation_weights: input.observation_weights,
            all_zero: input.all_zero,
            dispersions: input.dispersions,
        })?;
        let mut wald = wald_test_coefficient_with_options(
            &glm_fit,
            input.coefficient,
            &self.wald_test_options,
        )?;
        mask_wald_degrees_of_freedom_for_all_zero_rows(&mut wald, input.all_zero)?;
        let cooks = calculate_cooks_distance(
            input.counts,
            input.normalized,
            &glm_fit.mu,
            &glm_fit.hat_diagonal,
            input.design,
        )?;
        let mut results = build_wald_results_from_wald(
            input.base_mean,
            &glm_fit,
            input.coefficient,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
            &wald,
        )?;
        results.apply_wald_test_options(&self.wald_test_options);
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.design.n_samples(),
            input.design.n_coefficients(),
        )?;
        if let Some(contrast) =
            self.model_frame_factor_level_contrast_for_coefficient(input.design, input.coefficient)
        {
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                input.counts,
                &cooks.cooks,
                contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(WaldPipelineOutput {
            glm_fit,
            wald,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_wald_contrast_components(
        &self,
        input: FixedDispersionGlmInput<'_>,
        normalized: &RowMajorMatrix<f64>,
        base_mean: &[f64],
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
    ) -> Result<WaldContrastPipelineOutput, DeseqError> {
        let FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        } = self.fixed_dispersion_glm_components(input)?;
        let mut wald_contrast =
            wald_test_contrast_with_options(&glm_fit, contrast, &self.wald_test_options)?;
        mask_wald_degrees_of_freedom_for_all_zero_rows(&mut wald_contrast.wald, input.all_zero)?;
        let contrast_all_zero = match contrast_all_zero_override {
            Some(flags) => {
                if flags.len() != input.counts.n_genes() {
                    return Err(invalid_dimensions(
                        "contrastAllZero rows",
                        input.counts.n_genes(),
                        flags.len(),
                    ));
                }
                flags.to_vec()
            }
            None => contrast_all_zero_numeric(input.counts, input.design, contrast)?,
        };
        apply_contrast_all_zero_to_wald_contrast(
            &mut wald_contrast,
            &contrast_all_zero,
            input.all_zero,
        )?;
        let cooks = calculate_cooks_distance(
            input.counts,
            normalized,
            &glm_fit.mu,
            &glm_fit.hat_diagonal,
            input.design,
        )?;
        let mut results = build_wald_contrast_results(
            base_mean,
            &glm_fit,
            &wald_contrast,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
        )?;
        results.apply_wald_test_options(&self.wald_test_options);
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.design.n_samples(),
            input.design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(WaldContrastPipelineOutput {
            glm_fit,
            wald_contrast,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_lrt_components(
        &self,
        input: LrtPipelineInput<'_>,
    ) -> Result<LrtPipelineOutput, DeseqError> {
        validate_lrt_pipeline_input(&input)?;
        let nonzero_gene_indices = nonzero_gene_indices(input.all_zero);
        let (full_fit, reduced_fit) = if nonzero_gene_indices.is_empty() {
            (
                all_zero_glm_fit(input.counts, input.full_design)?,
                all_zero_glm_fit(input.counts, input.reduced_design)?,
            )
        } else {
            let compact_counts = compact_counts(input.counts, &nonzero_gene_indices)?;
            let compact_dispersions = nonzero_gene_indices
                .iter()
                .map(|gene| input.dispersions[*gene])
                .collect::<Vec<_>>();
            let compact_normalization_factors = input
                .normalization_factors
                .map(|factors| compact_matrix_rows(factors, &nonzero_gene_indices))
                .transpose()?;
            let compact_weights = input
                .observation_weights
                .map(|weights| compact_matrix_rows(weights, &nonzero_gene_indices))
                .transpose()?;
            let compact_full_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.full_design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            let compact_reduced_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.reduced_design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            (
                expand_glm_fit(compact_full_fit, input.all_zero)?,
                expand_glm_fit(compact_reduced_fit, input.all_zero)?,
            )
        };

        let expanded_dispersions =
            mask_all_zero_values_with_nan_rows(input.dispersions, input.all_zero)?;
        let lrt = lrt_test(&full_fit, &reduced_fit)?;
        let cooks = calculate_cooks_distance(
            input.counts,
            input.normalized,
            &full_fit.mu,
            &full_fit.hat_diagonal,
            input.full_design,
        )?;
        let mut results = build_lrt_results(
            input.base_mean,
            &full_fit,
            &lrt,
            input.coefficient,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
        )?;
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.full_design.n_samples(),
            input.full_design.n_coefficients(),
        )?;
        if let Some(contrast) = self
            .model_frame_factor_level_contrast_for_coefficient(input.full_design, input.coefficient)
        {
            apply_cooks_cutoff_for_factor_level_metadata(
                &mut results,
                cooks_cutoff,
                input.counts,
                &cooks.cooks,
                contrast,
            )?;
        } else {
            apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        }
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(LrtPipelineOutput {
            full_fit,
            reduced_fit,
            lrt,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_glm_components(
        &self,
        input: FixedDispersionGlmInput<'_>,
    ) -> Result<FixedDispersionGlmOutput, DeseqError> {
        input.design.validate_full_rank("GLM")?;
        if input.dispersions.len() != input.counts.n_genes() {
            return Err(invalid_dimensions(
                "pipeline dispersions",
                input.counts.n_genes(),
                input.dispersions.len(),
            ));
        }

        let nonzero_gene_indices = nonzero_gene_indices(input.all_zero);
        let glm_fit = if nonzero_gene_indices.is_empty() {
            all_zero_glm_fit(input.counts, input.design)?
        } else {
            let compact_counts = compact_counts(input.counts, &nonzero_gene_indices)?;
            let compact_dispersions = nonzero_gene_indices
                .iter()
                .map(|gene| input.dispersions[*gene])
                .collect::<Vec<_>>();
            let compact_normalization_factors = input
                .normalization_factors
                .map(|factors| compact_matrix_rows(factors, &nonzero_gene_indices))
                .transpose()?;
            let compact_weights = input
                .observation_weights
                .map(|weights| compact_matrix_rows(weights, &nonzero_gene_indices))
                .transpose()?;
            let compact_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            expand_glm_fit(compact_fit, input.all_zero)?
        };
        let expanded_dispersions =
            mask_all_zero_values_with_nan_rows(input.dispersions, input.all_zero)?;
        Ok(FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        })
    }

    fn base_fit(
        &self,
        counts: &CountMatrix,
        design: Option<DesignMatrix>,
        input: BaseFitInput,
    ) -> DeseqFit {
        DeseqFit {
            counts_summary: counts.summary(),
            design,
            reduced_design: None,
            model_frame: self.model_frame.clone(),
            size_factors: input.size_factors,
            normalization_factors: input.normalization_factors,
            observation_weights: input.observation_weights,
            weights_fail: input.weights_fail,
            weights_design_rank: input.weights_design_rank,
            base_mean: input.base_mean,
            base_var: input.base_var,
            all_zero: input.all_zero,
            disp_gene_est: None,
            disp_gene_iter: None,
            disp_fit: None,
            dispersion_trend: None,
            disp_map: None,
            dispersion: None,
            disp_iter: None,
            disp_outlier: None,
            disp_prior_var: None,
            var_log_disp_estimates: None,
            dispersion_converged: None,
            beta: None,
            beta_se: None,
            beta_optim_start: None,
            beta_covariance: None,
            beta_converged: None,
            beta_iter: None,
            beta_optim_iter: None,
            beta_optim_start_objective: None,
            beta_optim_objective: None,
            beta_optim_gradient_norm: None,
            log_like: None,
            full_deviance: None,
            reduced_log_like: None,
            reduced_beta_converged: None,
            reduced_beta_iter: None,
            reduced_mu: None,
            reduced_hat_diagonal: None,
            mu: None,
            cooks: None,
            max_cooks: None,
            hat_diagonal: None,
            wald: None,
            lrt: None,
        }
    }
}
