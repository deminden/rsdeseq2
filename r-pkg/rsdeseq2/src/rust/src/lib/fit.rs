// Wald/LRT fitting export, replacement and fitted output buffers.

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return. The
/// replacement output is a distinct writable owner-pointer slot. On success it
/// receives NULL or an owned result which C must drop on success or R unwinding.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_fit(
    n_genes: usize,
    n_samples: usize,
    counts_values: *const f64,
    size_factors_values: *const f64,
    normalization_factors_values: *const f64,
    observation_weights_values: *const f64,
    n_coefficients: usize,
    design_values: *const f64,
    reduced_coefficients: usize,
    reduced_values: *const f64,
    fit_type: i32,
    test_type: i32,
    min_replicates_for_replace: usize,
    cooks_replacement_cutoff: f64,
    threads: usize,
    execution_mode: i32,
    result_output: *mut f64,
    beta_output: *mut f64,
    beta_se_output: *mut f64,
    covariance_output: *mut f64,
    mu_output: *mut f64,
    cooks_output: *mut f64,
    replacement_output: *mut *mut replacement_counts::ReplacementCounts,
    replace_output: *mut u8,
    replaceable_samples_output: *mut u8,
    refit_rows_output: *mut u8,
    new_all_zero_rows_output: *mut u8,
    n_refit_output: *mut usize,
    size_factors_output: *mut f64,
    dispersion_output: *mut f64,
    disp_fit_output: *mut f64,
    base_var_output: *mut f64,
    all_zero_output: *mut i32,
    lrt_df_output: *mut i32,
) -> i32 {
    let mut pending = None;
    let status = run_backend(|| {
        if replacement_output.is_null() { return Err("replacement owner output pointer is null".into()); }
        unsafe { *replacement_output = ptr::null_mut(); }
        let lrt_df = unsafe { output_slice(lrt_df_output, 1, "LRT degrees of freedom")? };
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts_values)? };
        let design_matrix =
            unsafe { row_major_from_column_major(n_samples, n_coefficients, design_values)? };
        let design = DesignMatrix::from_row_major(
            n_samples,
            n_coefficients,
            design_matrix.into_values(),
            None,
        )
        .map_err(|error| error.to_string())?;
        let selected_fit_type = match fit_type {
            0 => FitType::Parametric,
            1 => FitType::Local,
            2 => FitType::Mean,
            _ => return Err("fitType must be parametric, local, or mean".into()),
        };
        let selected_test_type = match test_type {
            0 => TestType::Wald,
            1 => TestType::Lrt,
            _ => return Err("test must be Wald or LRT".into()),
        };
        let selected_execution_mode = match execution_mode {
            0 => ExecutionMode::Strict,
            1 => ExecutionMode::Fast,
            2 => ExecutionMode::Diagnostic,
            _ => return Err("parityMode must be strict, fast, or diagnostic".into()),
        };
        let mut builder = DeseqBuilder::new()
            .automatic_dispersion_means()
            .fit_type(selected_fit_type)
            .test(selected_test_type)
            .execution_mode(selected_execution_mode)
            .threads(threads.max(1));
        if !size_factors_values.is_null() {
            builder = builder.size_factors(
                unsafe { input_slice(size_factors_values, n_samples, "size factors")? }.to_vec(),
            );
        }
        if !normalization_factors_values.is_null() {
            let factors = unsafe {
                row_major_from_column_major(n_genes, n_samples, normalization_factors_values)?
            };
            builder = builder.normalization_factors(factors);
        }
        if !observation_weights_values.is_null() {
            let weights = unsafe {
                row_major_from_column_major(n_genes, n_samples, observation_weights_values)?
            };
            builder = builder.observation_weights(weights);
        }
        if selected_test_type == TestType::Lrt {
            let reduced_matrix = unsafe {
                row_major_from_column_major(n_samples, reduced_coefficients, reduced_values)?
            };
            let reduced = DesignMatrix::from_row_major(
                n_samples,
                reduced_coefficients,
                reduced_matrix.into_values(),
                None,
            )
            .map_err(|error| error.to_string())?;
            builder = builder.reduced_design(reduced);
        }
        let (fit, results, replacement) = if min_replicates_for_replace == 0 {
            let result = builder.fit_with_results(&counts, &design)
                .map_err(|error| error.to_string())?;
            (result.0, result.1, None)
        } else {
            let mut replacement_options = CooksReplacementOptions::new(cooks_replacement_cutoff);
            replacement_options.min_replicates = min_replicates_for_replace;
            let output = builder
                .fit_with_test_results_with_cooks_replacement(
                    &counts,
                    &design,
                    &replacement_options,
                )
                .map_err(|error| error.to_string())?;
            let (merged_fit, results, plan) = output.into_merged_fit()
                .map_err(|error| error.to_string())?;
            (merged_fit, results, Some(plan))
        };
        lrt_df[0] = fit.lrt.as_ref().map(|lrt| i32::try_from(lrt.degrees_of_freedom)
            .map_err(|_| "LRT degrees of freedom exceed the R integer range"))
            .transpose()?.unwrap_or(0);
        const COLUMNS: usize = 13;
        let output = unsafe {
            output_slice(
                result_output,
                n_genes
                    .checked_mul(COLUMNS)
                    .ok_or("result dimensions overflow")?,
                "result table",
            )?
        };
        let mut newly_zero = vec![false; n_genes];
        if let Some(replacement) = &replacement {
            for &gene in &replacement.new_all_zero_rows {
                newly_zero[gene] = true;
            }
        }
        for (gene, row) in results.rows.iter().enumerate() {
            let raw_fit = &fit;
            let (raw_stat, raw_pvalue) = if newly_zero[gene] {
                (Some(0.0), Some(1.0))
            } else if fit.all_zero[gene] {
                (None, None)
            } else if let Some(lrt) = &raw_fit.lrt {
                (lrt.deviance[gene], lrt.pvalue[gene])
            } else if let Some(wald) = &raw_fit.wald {
                (wald.stat[gene], wald.pvalue[gene])
            } else {
                return Err("fit is missing raw test statistics".into());
            };
            let values = [
                row.base_mean,
                row.log2_fold_change.unwrap_or(f64::NAN),
                row.lfc_se.unwrap_or(f64::NAN),
                row.stat.unwrap_or(f64::NAN),
                row.pvalue.unwrap_or(f64::NAN),
                row.padj.unwrap_or(f64::NAN),
                row.dispersion
                    .unwrap_or_else(|| fit.dispersion.as_ref().map_or(f64::NAN, |d| d[gene])),
                row.converged
                    .map_or(f64::NAN, |value| f64::from(value as u8)),
                row.max_cooks.unwrap_or(f64::NAN),
                row.cooks_outlier
                    .map_or(f64::NAN, |value| f64::from(value as u8)),
                row.filtered
                    .map_or(f64::NAN, |value| f64::from(value as u8)),
                raw_pvalue.unwrap_or(f64::NAN),
                raw_stat.unwrap_or(f64::NAN),
            ];
            for (column, value) in values.into_iter().enumerate() {
                output[gene + column * n_genes] = value;
            }
        }
        let beta = fit
            .beta
            .as_ref()
            .ok_or("Rust fit did not retain coefficient estimates")?;
        let beta_se = fit
            .beta_se
            .as_ref()
            .ok_or("Rust fit did not retain coefficient standard errors")?;
        let covariance = fit
            .beta_covariance
            .as_ref()
            .ok_or("Rust fit did not retain coefficient covariance")?;
        let mu = fit
            .mu
            .as_ref()
            .ok_or("Rust fit did not retain fitted means")?;
        let cooks = fit
            .cooks
            .as_ref()
            .ok_or("Rust fit did not retain Cook's distances")?;
        let beta_out = unsafe {
            output_slice(
                beta_output,
                n_genes
                    .checked_mul(n_coefficients)
                    .ok_or("beta dimensions overflow")?,
                "beta matrix",
            )?
        };
        let beta_se_out = unsafe {
            output_slice(
                beta_se_output,
                n_genes
                    .checked_mul(n_coefficients)
                    .ok_or("beta-SE dimensions overflow")?,
                "beta-SE matrix",
            )?
        };
        let covariance_len = n_genes
            .checked_mul(n_coefficients)
            .and_then(|value| value.checked_mul(n_coefficients))
            .ok_or("covariance dimensions overflow")?;
        let covariance_out =
            unsafe { output_slice(covariance_output, covariance_len, "beta covariance")? };
        for gene in 0..n_genes {
            let beta_row = beta.row(gene).map_err(|error| error.to_string())?;
            let beta_se_row = beta_se.row(gene).map_err(|error| error.to_string())?;
            let covariance_row = covariance.row(gene).map_err(|error| error.to_string())?;
            for coefficient in 0..n_coefficients {
                beta_out[gene + coefficient * n_genes] = beta_row[coefficient];
                beta_se_out[gene + coefficient * n_genes] = beta_se_row[coefficient];
                for other in 0..n_coefficients {
                    let r_offset = coefficient * n_coefficients + other;
                    let out_offset =
                        gene + coefficient * n_genes + other * n_genes * n_coefficients;
                    covariance_out[out_offset] = covariance_row[r_offset];
                }
            }
        }
        unsafe { copy_matrix_to_r(mu, mu_output)? };
        unsafe { copy_matrix_to_r(cooks, cooks_output)? };
        let replace = unsafe { output_slice(replace_output, n_genes, "replacement flags")? };
        let replaceable_samples =
            unsafe { output_slice(replaceable_samples_output, n_samples, "replaceable samples")? };
        let refit_rows = unsafe { output_slice(refit_rows_output, n_genes, "refit rows")? };
        let new_all_zero_rows =
            unsafe { output_slice(new_all_zero_rows_output, n_genes, "new all-zero rows")? };
        let n_refit = unsafe { output_slice(n_refit_output, 1, "replacement refit count")? };
        if let Some(output) = &replacement {
            let plan = output;
            let has_replaceable_samples = plan.has_replaceable_samples();
            n_refit[0] = if has_replaceable_samples {
                plan.n_refit
            } else {
                0
            };
            for (gene, flag) in replace.iter_mut().enumerate() {
                *flag = if has_replaceable_samples {
                    match plan.replacement.replace[gene] {
                        Some(value) => u8::from(value),
                        None => 2,
                    }
                } else {
                    2
                };
            }
            for (sample, flag) in replaceable_samples.iter_mut().enumerate() {
                *flag = u8::from(plan.replacement.replaceable_samples[sample]);
            }
            refit_rows.fill(0);
            new_all_zero_rows.fill(0);
            if has_replaceable_samples {
                for &row in &plan.refit_rows {
                    refit_rows[row] = 1;
                }
                for &row in &plan.new_all_zero_rows {
                    new_all_zero_rows[row] = 1;
                }
            }
        } else {
            n_refit[0] = 0;
            replace.fill(2);
            replaceable_samples.fill(0);
            refit_rows.fill(0);
            new_all_zero_rows.fill(0);

        }
        unsafe { output_slice(size_factors_output, n_samples, "size factors")? }
            .copy_from_slice(&fit.size_factors);
        let dispersions = unsafe { output_slice(dispersion_output, n_genes, "dispersions")? };
        match &fit.dispersion {
            Some(values) => dispersions.copy_from_slice(values),
            None => dispersions.fill(f64::NAN),
        }
        let disp_fit = unsafe { output_slice(disp_fit_output, n_genes, "dispersion fit")? };
        match &fit.disp_fit {
            Some(values) => disp_fit.copy_from_slice(values),
            None => disp_fit.fill(f64::NAN),
        }
        unsafe { output_slice(base_var_output, n_genes, "base variances")? }
            .copy_from_slice(&fit.base_var);
        let all_zero = unsafe { output_slice(all_zero_output, n_genes, "all-zero flags")? };
        for (target, source) in all_zero.iter_mut().zip(&fit.all_zero) {
            *target = i32::from(*source);
        }
        if n_refit[0] > 0 {
            let plan = replacement.ok_or("replacement refit count has no plan")?;
            pending = Some(Box::new(replacement_counts::ReplacementCounts::new(plan.replacement.replaced_counts)));
        }
        Ok(())
    });
    // Publish only after borrowed R buffers and the fit frame have been released.
    if status == 0 {
        unsafe { *replacement_output = pending.map_or(ptr::null_mut(), Box::into_raw); }
    }
    status
}
