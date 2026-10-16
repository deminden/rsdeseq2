// Expanded-model contrast shrinkage export.

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_lfc_shrink_normal_contrast(
    n_genes: usize,
    n_samples: usize,
    counts_values: *const f64,
    size_factors_values: *const f64,
    normalization_factors_values: *const f64,
    observation_weights_values: *const f64,
    n_coefficients: usize,
    expanded_design_values: *const f64,
    n_standard_coefficients: usize,
    standard_beta_values: *const f64,
    standard_factor_groups_values: *const i32,
    expanded_prior_map_values: *const i32,
    dispersions_values: *const f64,
    base_mean_values: *const f64,
    disp_fit_values: *const f64,
    factor_groups_values: *const i32,
    contrast_values: *const f64,
    intercept_index: isize,
    standard_intercept_index: isize,
    lfc_threshold: f64,
    beta_output: *mut f64,
    beta_se_output: *mut f64,
    prior_variance_output: *mut f64,
    stat_output: *mut f64,
    pvalue_output: *mut f64,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts_values)? };
        let design_values = unsafe {
            row_major_from_column_major(n_samples, n_coefficients, expanded_design_values)?
        };
        let mut coefficient_names = (0..n_coefficients)
            .map(|index| format!("coefficient{index}"))
            .collect::<Vec<_>>();
        if intercept_index >= 0 {
            let index = intercept_index as usize;
            if index >= n_coefficients {
                return Err("intercept index is outside the expanded design matrix".into());
            }
            coefficient_names[index] = "Intercept".to_string();
        }
        let design = DesignMatrix::from_row_major(
            n_samples,
            n_coefficients,
            design_values.into_values(),
            Some(coefficient_names.clone()),
        )
        .map_err(|error| error.to_string())?;
        let dispersions = unsafe { input_slice(dispersions_values, n_genes, "dispersions")? };
        let base_mean = unsafe { input_slice(base_mean_values, n_genes, "base means")? };
        let disp_fit = unsafe { input_slice(disp_fit_values, n_genes, "dispersion fit")? };
        let factor_groups = unsafe {
            input_slice(
                factor_groups_values,
                n_coefficients,
                "factor coefficient groups",
            )?
        };
        let standard_factor_groups = unsafe {
            input_slice(
                standard_factor_groups_values,
                n_standard_coefficients,
                "standard-design factor coefficient groups",
            )?
        };
        let expanded_prior_map = unsafe {
            input_slice(
                expanded_prior_map_values,
                n_coefficients,
                "expanded prior map",
            )?
        };
        if standard_factor_groups.iter().any(|group| *group < 0)
            || expanded_prior_map.iter().any(|column| *column < 0)
        {
            return Err("factor groups and prior-map indices must be non-negative".into());
        }
        let contrast = unsafe { input_slice(contrast_values, n_coefficients, "contrast")? };
        if factor_groups.iter().any(|group| *group < 0) {
            return Err("factor coefficient groups must be non-negative".into());
        }
        if !lfc_threshold.is_finite() || lfc_threshold < 0.0 {
            return Err("lfcThreshold must be finite and non-negative".into());
        }
        let normalization_factors = if !normalization_factors_values.is_null() {
            unsafe {
                row_major_from_column_major(n_genes, n_samples, normalization_factors_values)?
            }
        } else {
            let size_factors =
                unsafe { input_slice(size_factors_values, n_samples, "size factors")? };
            normalization_factors_from_size_factors(&counts, size_factors)
                .map_err(|error| error.to_string())?
        };
        let weights = if observation_weights_values.is_null() {
            None
        } else {
            let weights = unsafe {
                row_major_from_column_major(n_genes, n_samples, observation_weights_values)?
            };
            Some(
                preprocess_observation_weights(&weights, &design)
                    .map_err(|error| error.to_string())?
                    .weights,
            )
        };

        let standard_beta = unsafe {
            row_major_from_column_major(n_genes, n_standard_coefficients, standard_beta_values)?
        };
        let mut augmented_groups = standard_factor_groups
            .iter()
            .map(|group| *group as usize)
            .collect::<Vec<_>>();
        let mut contrast_pairs = Vec::new();
        let max_group = augmented_groups.iter().copied().max().unwrap_or(0);
        for group in 1..=max_group {
            let columns = augmented_groups
                .iter()
                .enumerate()
                .filter_map(|(column, candidate)| (*candidate == group).then_some(column))
                .collect::<Vec<_>>();
            if columns.is_empty() {
                return Err(format!(
                    "standard-design factor coefficient group {group} is empty"
                ));
            }
            for left in 0..columns.len() {
                for right in left + 1..columns.len() {
                    contrast_pairs.push((group, columns[left], columns[right]));
                    augmented_groups.push(group);
                }
            }
        }
        let augmented_columns = n_standard_coefficients + contrast_pairs.len();
        let mut augmented_values = Vec::with_capacity(n_genes * augmented_columns);
        for gene in 0..n_genes {
            let beta_row = standard_beta.row(gene).map_err(|error| error.to_string())?;
            augmented_values.extend_from_slice(beta_row);
            for (_, left, right) in &contrast_pairs {
                augmented_values.push(beta_row[*left] - beta_row[*right]);
            }
        }
        let augmented_beta =
            RowMajorMatrix::from_row_major(n_genes, augmented_columns, augmented_values)
                .map_err(|error| error.to_string())?;
        let mut standard_names = (0..n_standard_coefficients)
            .map(|index| format!("standardCoefficient{index}"))
            .collect::<Vec<_>>();
        if standard_intercept_index >= 0 {
            let index = standard_intercept_index as usize;
            if index >= n_standard_coefficients {
                return Err("intercept index is outside the standard design matrix".into());
            }
            standard_names[index] = "Intercept".to_string();
        }
        let augmented_names = standard_names
            .iter()
            .cloned()
            .chain((0..contrast_pairs.len()).map(|index| format!("factorContrast{index}")))
            .collect::<Vec<_>>();
        let estimated_priors = estimate_beta_prior_variance(
            &augmented_beta,
            base_mean,
            disp_fit,
            Some(&augmented_names),
            BetaPriorVarianceOptions::default(),
        )
        .map_err(|error| error.to_string())?;
        let mut prior_variance = vec![f64::NAN; n_coefficients];
        for (expanded_column, standard_column) in expanded_prior_map.iter().enumerate() {
            if *standard_column > 0 {
                let standard_column = (*standard_column - 1) as usize;
                if standard_column >= n_standard_coefficients {
                    return Err("expanded prior map refers outside the standard design".into());
                }
                prior_variance[expanded_column] = estimated_priors[standard_column];
            }
        }
        for group in 1..=max_group {
            let members = augmented_groups
                .iter()
                .enumerate()
                .filter_map(|(column, candidate)| (*candidate == group).then_some(column))
                .collect::<Vec<_>>();
            if members.is_empty()
                || !factor_groups
                    .iter()
                    .any(|candidate| *candidate as usize == group)
            {
                return Err(format!(
                    "expanded-design factor coefficient group {group} is empty"
                ));
            }
            let mean = members
                .iter()
                .map(|column| estimated_priors[*column])
                .sum::<f64>()
                / members.len() as f64;
            for (column, candidate) in factor_groups.iter().enumerate() {
                if *candidate as usize == group {
                    prior_variance[column] = mean;
                }
            }
        }
        if prior_variance
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(
                "expanded model prior variances could not be mapped to every coefficient".into(),
            );
        }
        let prior_fit = fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            &design,
            &normalization_factors,
            dispersions,
            weights.as_ref(),
            &prior_variance,
            Default::default(),
        )
        .map_err(|error| error.to_string())?;
        let wald = wald_test_contrast_with_options(
            &prior_fit,
            contrast,
            &WaldTestOptions::normal()
                .with_lfc_threshold(lfc_threshold, WaldAlternative::GreaterAbs),
        )
        .map_err(|error| error.to_string())?;
        unsafe { output_slice(beta_output, n_genes, "shrunken contrast estimates")? }
            .iter_mut()
            .zip(wald.log2_fold_change.iter())
            .for_each(|(out, value)| *out = value.unwrap_or(f64::NAN));
        unsafe { output_slice(beta_se_output, n_genes, "shrunken contrast standard errors")? }
            .iter_mut()
            .zip(wald.lfc_se.iter())
            .for_each(|(out, value)| *out = value.unwrap_or(f64::NAN));
        unsafe {
            output_slice(
                prior_variance_output,
                n_coefficients,
                "beta prior variances",
            )?
        }
        .copy_from_slice(&prior_variance);
        unsafe { output_slice(stat_output, n_genes, "shrunken Wald statistics")? }
            .iter_mut()
            .zip(wald.wald.stat.iter())
            .for_each(|(out, value)| *out = value.unwrap_or(f64::NAN));
        unsafe { output_slice(pvalue_output, n_genes, "shrunken Wald p-values")? }
            .iter_mut()
            .zip(wald.wald.pvalue.iter())
            .for_each(|(out, value)| *out = value.unwrap_or(f64::NAN));
        Ok(())
    })
}
