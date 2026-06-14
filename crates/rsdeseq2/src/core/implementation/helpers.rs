fn default_results_coefficient(design: &DesignMatrix) -> Result<usize, DeseqError> {
    design
        .n_coefficients()
        .checked_sub(1)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "default results coefficient".to_string(),
            expected: 1,
            actual: 0,
        })
}

fn validate_observation_weights_for_counts(
    counts: &CountMatrix,
    weights: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
        return Err(DeseqError::InvalidDimensions {
            context: "observation weights".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: weights.len(),
        });
    }
    Ok(())
}

fn validate_lrt_pipeline_input(input: &LrtPipelineInput<'_>) -> Result<(), DeseqError> {
    if input.dispersions.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "pipeline dispersions",
            input.counts.n_genes(),
            input.dispersions.len(),
        ));
    }
    if input.full_design.n_samples() != input.reduced_design.n_samples() {
        return Err(invalid_dimensions(
            "LRT reduced design samples",
            input.full_design.n_samples(),
            input.reduced_design.n_samples(),
        ));
    }
    if input.full_design.n_coefficients() <= input.reduced_design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT full/reduced coefficients".to_string(),
            expected: input.reduced_design.n_coefficients() + 1,
            actual: input.full_design.n_coefficients(),
        });
    }
    if input.coefficient >= input.full_design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT result coefficient index".to_string(),
            expected: input.full_design.n_coefficients().saturating_sub(1),
            actual: input.coefficient,
        });
    }
    if input.all_zero.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "LRT all-zero rows",
            input.counts.n_genes(),
            input.all_zero.len(),
        ));
    }
    if input.base_mean.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "LRT baseMean rows",
            input.counts.n_genes(),
            input.base_mean.len(),
        ));
    }
    if input.normalized.n_rows() != input.counts.n_genes()
        || input.normalized.n_cols() != input.counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT normalized counts".to_string(),
            expected: input.counts.n_genes() * input.counts.n_samples(),
            actual: input.normalized.len(),
        });
    }
    input.full_design.validate_full_rank("LRT full")?;
    input.reduced_design.validate_full_rank("LRT reduced")?;
    Ok(())
}

fn is_intercept_only_design(design: &DesignMatrix) -> bool {
    design.n_coefficients() == 1
        && design
            .matrix()
            .as_slice()
            .iter()
            .all(|value| (*value - 1.0).abs() <= f64::EPSILON)
}

fn fit_fixed_dispersion_model(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    weights: Option<&RowMajorMatrix<f64>>,
    dispersions: &[f64],
    irls_options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    if is_intercept_only_design(design)
        && irls_options.uses_intercept_shortcut_for_coefficients(design.n_coefficients())?
    {
        match normalization_factors {
            Some(factors) => fit_intercept_only_fixed_dispersion_with_normalization_factors(
                counts,
                factors,
                dispersions,
                weights,
            ),
            None => fit_intercept_only_fixed_dispersion_with_weights(
                counts,
                size_factors,
                dispersions,
                weights,
            ),
        }
    } else {
        match normalization_factors {
            Some(factors) => fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
                counts,
                design,
                factors,
                dispersions,
                weights,
                irls_options,
            ),
            None => fit_fixed_dispersion_irls_with_weights(
                counts,
                design,
                size_factors,
                dispersions,
                weights,
                irls_options,
            ),
        }
    }
}

fn merge_replacement_refit_results(
    original_results: &DeseqResults,
    refit_results: Option<&DeseqResults>,
    refit_plan: &CooksRefitPlan,
) -> Result<DeseqResults, DeseqError> {
    if original_results.rows.len() != refit_plan.replacement.replace.len() {
        return Err(invalid_dimensions(
            "replacement-refit result rows",
            refit_plan.replacement.replace.len(),
            original_results.rows.len(),
        ));
    }
    if let Some(refit_results) = refit_results {
        if refit_results.rows.len() != original_results.rows.len() {
            return Err(invalid_dimensions(
                "replacement-refit refit result rows",
                original_results.rows.len(),
                refit_results.rows.len(),
            ));
        }
    }
    if refit_plan.replaced_base_mean.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "replacement-refit baseMean rows",
            original_results.rows.len(),
            refit_plan.replaced_base_mean.len(),
        ));
    }
    if refit_plan.post_refit_max_cooks.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "replacement-refit maxCooks rows",
            original_results.rows.len(),
            refit_plan.post_refit_max_cooks.len(),
        ));
    }

    // DESeq2 replacement refits preserve the original table shape: only rows
    // selected for refit get statistical values from replacement counts, while
    // every row receives replacement-aware baseMean/maxCooks metadata.
    let mut merged = original_results.clone();
    for (gene, row) in merged.rows.iter_mut().enumerate() {
        row.base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            row.max_cooks = refit_plan.post_refit_max_cooks[gene];
            row.cooks_outlier = None;
            row.filtered = None;
        }
    }

    if let Some(refit_results) = refit_results {
        for gene in refit_plan.refit_rows.iter().copied() {
            merged.rows[gene] = refit_results.rows[gene].clone();
            merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
            merged.rows[gene].cooks_outlier = None;
            merged.rows[gene].filtered = None;
        }
    }

    for gene in refit_plan.new_all_zero_rows.iter().copied() {
        clear_replacement_all_zero_result(&mut merged.rows[gene]);
        merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
        }
    }

    merged.independent_filtering = None;
    Ok(merged)
}

fn replacement_refit_plan_from_original(
    counts: &CountMatrix,
    design: &DesignMatrix,
    original_fit: &DeseqFit,
    replacement_options: &CooksReplacementOptions,
) -> Result<CooksRefitPlan, DeseqError> {
    let original_cooks = original_fit
        .cooks
        .as_ref()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "Cook's distances are required before replacement refit".to_string(),
        })?;
    let normalized = match original_fit.normalization_factors.as_ref() {
        Some(normalization_factors) => {
            normalized_counts_with_factors(counts, normalization_factors)?
        }
        None => normalized_counts(counts, &original_fit.size_factors)?,
    };
    prepare_cooks_replacement_refit(
        counts,
        &normalized,
        &original_fit.size_factors,
        original_fit.normalization_factors.as_ref(),
        original_cooks,
        design,
        replacement_options,
    )
}

fn replacement_dispersion_inputs(
    original_fit: &DeseqFit,
) -> Result<(&DispersionTrendFit, f64, f64), DeseqError> {
    let trend =
        original_fit
            .dispersion_trend
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "original dispersion function is required before replacement refit"
                    .to_string(),
            })?;
    let disp_prior_var =
        original_fit
            .disp_prior_var
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "original dispersion prior variance is required before replacement refit"
                    .to_string(),
            })?;
    let var_log_disp_estimates =
        original_fit
            .var_log_disp_estimates
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "original log-dispersion variance is required before replacement refit"
                    .to_string(),
            })?;
    Ok((trend, disp_prior_var, var_log_disp_estimates))
}

fn apply_contrast_metadata_to_replacement_output(
    output: &mut CooksReplacementWaldOutput,
    result_name: String,
    comparison: String,
    contrast: Option<&[f64]>,
) {
    output.original_results.set_resolved_contrast_metadata(
        result_name.clone(),
        comparison.clone(),
        contrast.unwrap_or(&[]),
    );
    if contrast.is_none() {
        output.original_results.metadata.contrast = None;
    }
    if let Some(refit_results) = &mut output.refit_results {
        refit_results.set_resolved_contrast_metadata(
            result_name.clone(),
            comparison.clone(),
            contrast.unwrap_or(&[]),
        );
        if contrast.is_none() {
            refit_results.metadata.contrast = None;
        }
    }
    output
        .results
        .set_resolved_contrast_metadata(result_name, comparison, contrast.unwrap_or(&[]));
    if contrast.is_none() {
        output.results.metadata.contrast = None;
    }
}

fn apply_lrt_contrast_metadata_to_replacement_output(
    output: &mut CooksReplacementLrtOutput,
    result_name: String,
    comparison: String,
    contrast: Option<&[f64]>,
) {
    output.original_results.set_resolved_contrast_metadata(
        result_name.clone(),
        comparison.clone(),
        contrast.unwrap_or(&[]),
    );
    if contrast.is_none() {
        output.original_results.metadata.contrast = None;
    }
    if let Some(refit_results) = &mut output.refit_results {
        refit_results.set_resolved_contrast_metadata(
            result_name.clone(),
            comparison.clone(),
            contrast.unwrap_or(&[]),
        );
        if contrast.is_none() {
            refit_results.metadata.contrast = None;
        }
    }
    output
        .results
        .set_resolved_contrast_metadata(result_name, comparison, contrast.unwrap_or(&[]));
    if contrast.is_none() {
        output.results.metadata.contrast = None;
    }
}

fn factor_level_result_metadata(contrast: FactorLevelContrast<'_>) -> (String, String) {
    let spec = match contrast.reference {
        Some(reference) => ContrastSpec::factor_level_with_reference(
            contrast.factor,
            contrast.numerator,
            contrast.denominator,
            reference,
        ),
        None => ContrastSpec::factor_level(
            contrast.factor,
            contrast.numerator,
            contrast.denominator,
        ),
    };
    (spec.result_name(), spec.comparison())
}

fn apply_cooks_cutoff_for_factor_level_metadata(
    results: &mut DeseqResults,
    cutoff: Option<f64>,
    counts: &CountMatrix,
    cooks: &RowMajorMatrix<f64>,
    contrast: FactorLevelContrast<'_>,
) -> Result<(), DeseqError> {
    if factor_level_contrast_is_single_two_level_condition(contrast) {
        apply_cooks_cutoff_with_low_count_heuristic(results, cutoff, counts, cooks)
    } else {
        apply_cooks_cutoff(results, cutoff)
    }
}

fn factor_level_contrast_is_single_two_level_condition(contrast: FactorLevelContrast<'_>) -> bool {
    let mut saw_numerator = false;
    let mut saw_denominator = false;
    for level in contrast.sample_levels {
        if level == contrast.numerator {
            saw_numerator = true;
        } else if level == contrast.denominator {
            saw_denominator = true;
        } else {
            return false;
        }
    }
    saw_numerator && saw_denominator
}

fn clear_replacement_all_zero_result(row: &mut DeseqResultRow) {
    row.log2_fold_change = Some(0.0);
    row.lfc_se = Some(0.0);
    row.stat = Some(0.0);
    row.pvalue = Some(1.0);
    row.padj = None;
    row.dispersion = None;
    row.converged = None;
    row.cooks_outlier = None;
    row.filtered = None;
}

fn compact_counts(counts: &CountMatrix, gene_indices: &[usize]) -> Result<CountMatrix, DeseqError> {
    counts.select_rows(gene_indices)
}

fn compact_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    row_indices: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * matrix.n_cols());
    for row in row_indices {
        values.extend_from_slice(matrix.row(*row)?);
    }
    RowMajorMatrix::from_row_major(row_indices.len(), matrix.n_cols(), values)
}

fn compact_f64_values(values: &[f64], row_indices: &[usize]) -> Result<Vec<f64>, DeseqError> {
    let mut compact = Vec::with_capacity(row_indices.len());
    for row in row_indices {
        let Some(value) = values.get(*row) else {
            return Err(invalid_dimensions(
                "compact vector rows",
                row + 1,
                values.len(),
            ));
        };
        compact.push(*value);
    }
    Ok(compact)
}

fn expand_rlog_output_with_all_zero_rows(
    compact_output: RlogOutput,
    all_zero: &[bool],
    n_samples: usize,
) -> Result<RlogOutput, DeseqError> {
    if compact_output.transformed.n_cols() != n_samples {
        return Err(invalid_dimensions(
            "expanded rlog columns",
            n_samples,
            compact_output.transformed.n_cols(),
        ));
    }
    let mut values = vec![0.0; all_zero.len() * n_samples];
    let mut compact_row = 0_usize;
    for (gene, is_zero) in all_zero.iter().copied().enumerate() {
        if is_zero {
            continue;
        }
        let src = compact_output.transformed.row(compact_row)?;
        let start = gene * n_samples;
        values[start..start + n_samples].copy_from_slice(src);
        compact_row += 1;
    }
    if compact_row != compact_output.transformed.n_rows() {
        return Err(invalid_dimensions(
            "expanded rlog non-zero rows",
            compact_row,
            compact_output.transformed.n_rows(),
        ));
    }
    Ok(RlogOutput {
        transformed: RowMajorMatrix::from_row_major(all_zero.len(), n_samples, values)?,
        intercept: expand_rlog_intercepts_with_all_zero_rows(&compact_output, all_zero)?,
        sample_prior_variance: compact_output.sample_prior_variance,
        offset_mode: compact_output.offset_mode,
    })
}

fn expand_rlog_intercepts_with_all_zero_rows(
    compact_output: &RlogOutput,
    all_zero: &[bool],
) -> Result<Vec<f64>, DeseqError> {
    let expected_nonzero = all_zero.iter().filter(|is_zero| !**is_zero).count();
    if compact_output.intercept.len() != expected_nonzero {
        return Err(invalid_dimensions(
            "expanded rlog intercepts",
            expected_nonzero,
            compact_output.intercept.len(),
        ));
    }
    let mut compact_row = 0usize;
    let mut intercept = Vec::with_capacity(all_zero.len());
    for is_zero in all_zero {
        if *is_zero {
            intercept.push(0.0);
        } else {
            intercept.push(compact_output.intercept[compact_row]);
            compact_row += 1;
        }
    }
    Ok(intercept)
}

fn expand_frozen_rlog_output_with_all_zero_rows(
    compact_output: RlogOutput,
    all_zero: &[bool],
    n_samples: usize,
    frozen_intercept: &[f64],
) -> Result<RlogOutput, DeseqError> {
    if frozen_intercept.len() != all_zero.len() {
        return Err(invalid_dimensions(
            "expanded frozen rlog intercepts",
            all_zero.len(),
            frozen_intercept.len(),
        ));
    }
    if compact_output.transformed.n_cols() != n_samples {
        return Err(invalid_dimensions(
            "expanded frozen rlog columns",
            n_samples,
            compact_output.transformed.n_cols(),
        ));
    }
    let mut values = Vec::with_capacity(all_zero.len() * n_samples);
    let mut compact_row = 0usize;
    for (gene, is_zero) in all_zero.iter().enumerate() {
        if *is_zero {
            values.extend(std::iter::repeat_n(frozen_intercept[gene], n_samples));
        } else {
            let src = compact_output.transformed.row(compact_row)?;
            values.extend_from_slice(src);
            compact_row += 1;
        }
    }
    if compact_row != compact_output.transformed.n_rows() {
        return Err(invalid_dimensions(
            "expanded frozen rlog non-zero rows",
            compact_row,
            compact_output.transformed.n_rows(),
        ));
    }
    Ok(RlogOutput {
        transformed: RowMajorMatrix::from_row_major(all_zero.len(), n_samples, values)?,
        intercept: frozen_intercept.to_vec(),
        sample_prior_variance: compact_output.sample_prior_variance,
        offset_mode: compact_output.offset_mode,
    })
}

fn all_zero_glm_fit(
    counts: &CountMatrix,
    design: &DesignMatrix,
) -> Result<NbinomGlmFit, DeseqError> {
    Ok(NbinomGlmFit {
        log_like: vec![f64::NAN; counts.n_genes()],
        beta_converged: vec![false; counts.n_genes()],
        beta: RowMajorMatrix::from_elem(counts.n_genes(), design.n_coefficients(), f64::NAN)?,
        beta_se: RowMajorMatrix::from_elem(counts.n_genes(), design.n_coefficients(), f64::NAN)?,
        beta_optim_start: RowMajorMatrix::from_elem(
            counts.n_genes(),
            design.n_coefficients(),
            f64::NAN,
        )?,
        beta_covariance: Some(RowMajorMatrix::from_elem(
            counts.n_genes(),
            design.n_coefficients() * design.n_coefficients(),
            f64::NAN,
        )?),
        mu: RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), f64::NAN)?,
        beta_iter: vec![0; counts.n_genes()],
        beta_optim_iter: vec![f64::NAN; counts.n_genes()],
        beta_optim_start_objective: vec![f64::NAN; counts.n_genes()],
        beta_optim_objective: vec![f64::NAN; counts.n_genes()],
        beta_optim_gradient_norm: vec![f64::NAN; counts.n_genes()],
        model_matrix: design.clone(),
        n_terms: design.n_coefficients(),
        hat_diagonal: RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), f64::NAN)?,
    })
}

fn expand_glm_fit(
    compact_fit: NbinomGlmFit,
    all_zero: &[bool],
) -> Result<NbinomGlmFit, DeseqError> {
    Ok(NbinomGlmFit {
        log_like: expand_gene_values_with_nan_rows(&compact_fit.log_like, all_zero)?,
        beta_converged: expand_gene_values_with_fill_rows(
            &compact_fit.beta_converged,
            all_zero,
            false,
        )?,
        beta: expand_matrix_with_nan_rows(&compact_fit.beta, all_zero)?,
        beta_se: expand_matrix_with_nan_rows(&compact_fit.beta_se, all_zero)?,
        beta_optim_start: expand_matrix_with_nan_rows(&compact_fit.beta_optim_start, all_zero)?,
        beta_covariance: compact_fit
            .beta_covariance
            .as_ref()
            .map(|matrix| expand_matrix_with_nan_rows(matrix, all_zero))
            .transpose()?,
        mu: expand_matrix_with_nan_rows(&compact_fit.mu, all_zero)?,
        beta_iter: expand_gene_values_with_fill_rows(&compact_fit.beta_iter, all_zero, 0)?,
        beta_optim_iter: expand_gene_values_with_nan_rows(&compact_fit.beta_optim_iter, all_zero)?,
        beta_optim_start_objective: expand_gene_values_with_nan_rows(
            &compact_fit.beta_optim_start_objective,
            all_zero,
        )?,
        beta_optim_objective: expand_gene_values_with_nan_rows(
            &compact_fit.beta_optim_objective,
            all_zero,
        )?,
        beta_optim_gradient_norm: expand_gene_values_with_nan_rows(
            &compact_fit.beta_optim_gradient_norm,
            all_zero,
        )?,
        model_matrix: compact_fit.model_matrix,
        n_terms: compact_fit.n_terms,
        hat_diagonal: expand_matrix_with_nan_rows(&compact_fit.hat_diagonal, all_zero)?,
    })
}

fn mask_wald_degrees_of_freedom_for_all_zero_rows(
    wald: &mut WaldOutput,
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let Some(degrees_of_freedom) = &mut wald.degrees_of_freedom else {
        return Ok(());
    };
    if degrees_of_freedom.len() != all_zero.len() {
        return Err(invalid_dimensions(
            "Wald degrees of freedom all-zero mask",
            all_zero.len(),
            degrees_of_freedom.len(),
        ));
    }
    for (df, is_all_zero) in degrees_of_freedom.iter_mut().zip(all_zero.iter().copied()) {
        if is_all_zero {
            *df = None;
        }
    }
    Ok(())
}

fn apply_contrast_all_zero_to_wald_contrast(
    contrast: &mut WaldContrastOutput,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let n_genes = contrast.log2_fold_change.len();
    if contrast_all_zero.len() != n_genes {
        return Err(invalid_dimensions(
            "contrastAllZero rows",
            n_genes,
            contrast_all_zero.len(),
        ));
    }
    if all_zero.len() != n_genes {
        return Err(invalid_dimensions("allZero rows", n_genes, all_zero.len()));
    }
    for gene in 0..n_genes {
        if contrast_all_zero[gene] && !all_zero[gene] {
            contrast.log2_fold_change[gene] = Some(0.0);
            contrast.wald.stat[gene] = Some(0.0);
            contrast.wald.pvalue[gene] = Some(1.0);
        }
    }
    Ok(())
}

fn apply_contrast_all_zero_to_lrt_results(
    results: &mut DeseqResults,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let n_genes = results.rows.len();
    if contrast_all_zero.len() != n_genes {
        return Err(invalid_dimensions(
            "contrastAllZero rows",
            n_genes,
            contrast_all_zero.len(),
        ));
    }
    if all_zero.len() != n_genes {
        return Err(invalid_dimensions("allZero rows", n_genes, all_zero.len()));
    }
    for gene in 0..n_genes {
        if contrast_all_zero[gene] && !all_zero[gene] {
            results.rows[gene].log2_fold_change = Some(0.0);
        }
    }
    Ok(())
}

fn attach_glm_fit(fit: &mut DeseqFit, glm_fit: NbinomGlmFit) {
    let full_deviance = glm_fit
        .log_like
        .iter()
        .map(|log_like| full_deviance_from_log_like(*log_like))
        .collect();
    fit.beta = Some(glm_fit.beta);
    fit.beta_se = Some(glm_fit.beta_se);
    fit.beta_optim_start = Some(glm_fit.beta_optim_start);
    fit.beta_covariance = glm_fit.beta_covariance;
    fit.beta_converged = Some(glm_fit.beta_converged);
    fit.beta_iter = Some(glm_fit.beta_iter);
    fit.beta_optim_iter = Some(glm_fit.beta_optim_iter);
    fit.beta_optim_start_objective = Some(glm_fit.beta_optim_start_objective);
    fit.beta_optim_objective = Some(glm_fit.beta_optim_objective);
    fit.beta_optim_gradient_norm = Some(glm_fit.beta_optim_gradient_norm);
    fit.log_like = Some(glm_fit.log_like);
    fit.full_deviance = Some(full_deviance);
    fit.mu = Some(glm_fit.mu);
    fit.hat_diagonal = Some(glm_fit.hat_diagonal);
}

fn full_deviance_from_log_like(log_like: f64) -> f64 {
    checked_product2(-2.0, log_like).unwrap_or(f64::NAN)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let deviance = left * right;
    if left.is_finite() && right.is_finite() && deviance.is_finite() {
        Some(deviance)
    } else {
        None
    }
}

fn validate_pipeline_wald_coefficient(
    design: &DesignMatrix,
    coefficient: usize,
) -> Result<(), DeseqError> {
    if coefficient >= design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "pipeline Wald coefficient index".to_string(),
            expected: design.n_coefficients().saturating_sub(1),
            actual: coefficient,
        });
    }
    Ok(())
}
