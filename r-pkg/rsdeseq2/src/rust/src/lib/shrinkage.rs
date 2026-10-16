// Coefficient shrinkage export.

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_lfc_shrink_normal(
    n_genes: usize,
    n_samples: usize,
    counts_values: *const f64,
    size_factors_values: *const f64,
    normalization_factors_values: *const f64,
    observation_weights_values: *const f64,
    n_coefficients: usize,
    design_values: *const f64,
    dispersions_values: *const f64,
    base_mean_values: *const f64,
    disp_fit_values: *const f64,
    coefficient: usize,
    intercept_index: isize,
    lfc_threshold: f64,
    beta_output: *mut f64,
    beta_se_output: *mut f64,
    prior_variance_output: *mut f64,
    stat_output: *mut f64,
    pvalue_output: *mut f64,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts_values)? };
        let design_values =
            unsafe { row_major_from_column_major(n_samples, n_coefficients, design_values)? };
        let mut coefficient_names = (0..n_coefficients)
            .map(|index| format!("coefficient{index}"))
            .collect::<Vec<_>>();
        if intercept_index >= 0 {
            let index = intercept_index as usize;
            if index >= n_coefficients {
                return Err("intercept index is outside the design matrix".into());
            }
            coefficient_names[index] = "Intercept".to_string();
        }
        if coefficient >= n_coefficients {
            return Err("coefficient index is outside the design matrix".into());
        }
        let design = DesignMatrix::from_row_major(
            n_samples,
            n_coefficients,
            design_values.into_values(),
            Some(coefficient_names),
        )
        .map_err(|error| error.to_string())?;
        let dispersions = unsafe { input_slice(dispersions_values, n_genes, "dispersions")? };
        let base_mean = unsafe { input_slice(base_mean_values, n_genes, "base means")? };
        let disp_fit = unsafe { input_slice(disp_fit_values, n_genes, "dispersion fit")? };
        let weights = if observation_weights_values.is_null() {
            None
        } else {
            let weights = unsafe {
                row_major_from_column_major(n_genes, n_samples, observation_weights_values)?
            };
            Some(
                preprocess_observation_weights(&weights, &design)
                    .map_err(|e| e.to_string())?
                    .weights,
            )
        };
        let fitted = if !normalization_factors_values.is_null() {
            let factors = unsafe {
                row_major_from_column_major(n_genes, n_samples, normalization_factors_values)?
            };
            fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
                &counts,
                &design,
                BetaPriorNormalizationFactorWeightInput {
                    normalization_factors: &factors,
                    weights: weights.as_ref(),
                },
                dispersions,
                base_mean,
                disp_fit,
                BetaPriorRefitOptions::default(),
            )
        } else {
            let size_factors =
                unsafe { input_slice(size_factors_values, n_samples, "size factors")? };
            if let Some(weights) = weights.as_ref() {
                fit_glms_with_estimated_beta_prior_variance_and_weights(
                    &counts,
                    &design,
                    BetaPriorSizeFactorWeightInput {
                        size_factors,
                        weights: Some(weights),
                    },
                    dispersions,
                    base_mean,
                    disp_fit,
                    BetaPriorRefitOptions::default(),
                )
            } else {
                fit_glms_with_estimated_beta_prior_variance(
                    &counts,
                    &design,
                    size_factors,
                    dispersions,
                    base_mean,
                    disp_fit,
                    BetaPriorRefitOptions::default(),
                )
            }
        }
        .map_err(|error| error.to_string())?;
        if !lfc_threshold.is_finite() || lfc_threshold < 0.0 {
            return Err("lfcThreshold must be finite and non-negative".into());
        }
        let wald = wald_test_coefficient_with_options(
            &fitted.prior_fit,
            coefficient,
            &WaldTestOptions::normal()
                .with_lfc_threshold(lfc_threshold, WaldAlternative::GreaterAbs),
        )
        .map_err(|error| error.to_string())?;
        let beta = unsafe { output_slice(beta_output, n_genes, "shrunken coefficients")? };
        let beta_se = unsafe {
            output_slice(
                beta_se_output,
                n_genes,
                "shrunken coefficient standard errors",
            )?
        };
        for gene in 0..n_genes {
            beta[gene] = *fitted
                .prior_fit
                .beta
                .get(gene, coefficient)
                .ok_or("missing shrunken coefficient")?;
            beta_se[gene] = *fitted
                .prior_fit
                .beta_se
                .get(gene, coefficient)
                .ok_or("missing shrunken coefficient SE")?;
        }
        unsafe {
            output_slice(
                prior_variance_output,
                n_coefficients,
                "beta prior variances",
            )?
        }
        .copy_from_slice(&fitted.beta_prior_variance);
        let stats = unsafe { output_slice(stat_output, n_genes, "shrunken Wald statistics")? };
        let pvalues = unsafe { output_slice(pvalue_output, n_genes, "shrunken Wald p-values")? };
        for gene in 0..n_genes {
            stats[gene] = wald.stat[gene].unwrap_or(f64::NAN);
            pvalues[gene] = wald.pvalue[gene].unwrap_or(f64::NAN);
        }
        Ok(())
    })
}
