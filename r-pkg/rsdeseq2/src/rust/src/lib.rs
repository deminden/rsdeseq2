use std::cell::RefCell;
use std::ffi::{CString, c_char};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;

use rsdeseq2::prelude::{
    BetaPriorNormalizationFactorWeightInput, BetaPriorRefitOptions, BetaPriorSizeFactorWeightInput,
    BetaPriorVarianceOptions, CooksReplacementOptions, CountMatrix, DeseqBuilder, DeseqResultRow,
    DeseqResults, DesignMatrix, ExecutionMode, FitType, ParametricDispersionTrend, RowMajorMatrix,
    SizeFactorMethod, TestType, WaldAlternative, WaldTestOptions, apply_cooks_cutoff,
    apply_cooks_cutoff_with_low_count_heuristic, base_mean, base_mean_with_weights, base_variance,
    base_variance_with_weights, estimate_beta_prior_variance, estimate_size_factors_with_options,
    fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights,
    fit_glms_with_estimated_beta_prior_variance,
    fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights,
    fit_glms_with_estimated_beta_prior_variance_and_weights, local_vst_inverse_size_factor_mean,
    local_vst_inverse_size_factor_mean_from_normalization_factors,
    normalization_factors_from_size_factors, normalized_counts, normalized_counts_with_factors,
    preprocess_observation_weights, vst_local_from_dispersion_grid, vst_mean, vst_parametric,
    wald_test_coefficient_with_options, wald_test_contrast_with_options,
};

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("unknown Rust backend error").unwrap());
}

fn set_error(error: impl std::fmt::Display) {
    let message = CString::new(error.to_string().replace('\0', " "))
        .unwrap_or_else(|_| CString::new("Rust backend failed").unwrap());
    LAST_ERROR.with(|slot| *slot.borrow_mut() = message);
}

#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_last_error() -> *const c_char {
    LAST_ERROR.with(|slot| slot.borrow().as_ptr())
}

fn run_backend(operation: impl FnOnce() -> Result<(), String>) -> i32 {
    match catch_unwind(AssertUnwindSafe(operation)) {
        Ok(Ok(())) => 0,
        Ok(Err(error)) => {
            set_error(error);
            1
        }
        Err(_) => {
            set_error("Rust backend panicked");
            2
        }
    }
}

unsafe fn input_slice<'a, T>(
    pointer: *const T,
    length: usize,
    name: &str,
) -> Result<&'a [T], String> {
    if length > 0 && pointer.is_null() {
        return Err(format!("{name} pointer is null"));
    }
    Ok(if length == 0 {
        &[]
    } else {
        unsafe { std::slice::from_raw_parts(pointer, length) }
    })
}

unsafe fn output_slice<'a, T>(
    pointer: *mut T,
    length: usize,
    name: &str,
) -> Result<&'a mut [T], String> {
    if length > 0 && pointer.is_null() {
        return Err(format!("{name} output pointer is null"));
    }
    Ok(if length == 0 {
        &mut []
    } else {
        unsafe { std::slice::from_raw_parts_mut(pointer, length) }
    })
}

unsafe fn counts_from_r(
    n_genes: usize,
    n_samples: usize,
    values: *const f64,
) -> Result<CountMatrix, String> {
    let column_major = unsafe {
        input_slice(
            values,
            n_genes
                .checked_mul(n_samples)
                .ok_or("count dimensions overflow")?,
            "counts",
        )?
    };
    let mut row_major = Vec::with_capacity(column_major.len());
    for gene in 0..n_genes {
        for sample in 0..n_samples {
            let value = column_major[gene + sample * n_genes];
            if !value.is_finite()
                || value < 0.0
                || value.fract() != 0.0
                || value > f64::from(u32::MAX)
            {
                return Err(format!(
                    "count at row {}, column {} is not a valid u32 count",
                    gene + 1,
                    sample + 1
                ));
            }
            row_major.push(value as u32);
        }
    }
    CountMatrix::from_row_major_u32(n_genes, n_samples, row_major)
        .map_err(|error| error.to_string())
}

unsafe fn row_major_from_column_major(
    n_rows: usize,
    n_cols: usize,
    values: *const f64,
) -> Result<RowMajorMatrix<f64>, String> {
    let column_major = unsafe {
        input_slice(
            values,
            n_rows
                .checked_mul(n_cols)
                .ok_or("matrix dimensions overflow")?,
            "matrix",
        )?
    };
    let mut row_major = Vec::with_capacity(column_major.len());
    for row in 0..n_rows {
        for col in 0..n_cols {
            row_major.push(column_major[row + col * n_rows]);
        }
    }
    RowMajorMatrix::from_row_major(n_rows, n_cols, row_major).map_err(|error| error.to_string())
}

unsafe fn copy_matrix_to_r(matrix: &RowMajorMatrix<f64>, output: *mut f64) -> Result<(), String> {
    let n_rows = matrix.n_rows();
    let n_cols = matrix.n_cols();
    let out = unsafe {
        output_slice(
            output,
            n_rows
                .checked_mul(n_cols)
                .ok_or("output dimensions overflow")?,
            "matrix",
        )?
    };
    for row in 0..n_rows {
        for col in 0..n_cols {
            out[row + col * n_rows] = *matrix.get(row, col).ok_or("matrix index out of bounds")?;
        }
    }
    Ok(())
}

unsafe fn copy_counts_to_r(counts: &CountMatrix, output: &mut [f64]) -> Result<(), String> {
    if output.len() != counts.n_genes().saturating_mul(counts.n_samples()) {
        return Err("replacement count output has an invalid length".into());
    }
    for gene in 0..counts.n_genes() {
        for sample in 0..counts.n_samples() {
            output[gene + sample * counts.n_genes()] =
                f64::from(counts.as_slice()[gene * counts.n_samples() + sample]);
        }
    }
    Ok(())
}

fn merge_vec<T: Clone>(
    original: &mut [T],
    refit: &[T],
    rows: &[usize],
    field: &str,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    if original.len() != refit.len() || rows.iter().any(|&row| row >= original.len()) {
        return Err(format!("replacement refit {field} dimensions do not match"));
    }
    for &row in rows {
        original[row] = refit[row].clone();
    }
    Ok(())
}

fn merge_optional_vec<T: Clone>(
    original: &mut Option<Vec<T>>,
    refit: &Option<Vec<T>>,
    rows: &[usize],
    field: &str,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    let target = original
        .as_mut()
        .ok_or_else(|| format!("original fit is missing {field}"))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| format!("replacement refit is missing {field}"))?;
    merge_vec(target, source, rows, field)
}

fn merge_option_matrix(
    original: &mut Option<RowMajorMatrix<f64>>,
    refit: &Option<RowMajorMatrix<f64>>,
    rows: &[usize],
    field: &str,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    let target = original
        .as_ref()
        .ok_or_else(|| format!("original fit is missing {field}"))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| format!("replacement refit is missing {field}"))?;
    if target.n_rows() != source.n_rows()
        || target.n_cols() != source.n_cols()
        || rows.iter().any(|&row| row >= target.n_rows())
    {
        return Err(format!("replacement refit {field} dimensions do not match"));
    }
    let mut values = target.as_slice().to_vec();
    for &row in rows {
        let start = row * target.n_cols();
        values[start..start + target.n_cols()]
            .clone_from_slice(&source.as_slice()[start..start + target.n_cols()]);
    }
    *original = Some(
        RowMajorMatrix::from_row_major(target.n_rows(), target.n_cols(), values)
            .map_err(|error| error.to_string())?,
    );
    Ok(())
}

fn merge_refit_fit_rows(
    original: &mut rsdeseq2::prelude::DeseqFit,
    refit: &rsdeseq2::prelude::DeseqFit,
    rows: &[usize],
) -> Result<(), String> {
    macro_rules! merge_required_vec {
        ($field:ident) => {
            merge_vec(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    macro_rules! merge_optional_vec {
        ($field:ident) => {
            merge_optional_vec(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    macro_rules! merge_matrix {
        ($field:ident) => {
            merge_option_matrix(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    merge_required_vec!(base_mean);
    merge_required_vec!(base_var);
    merge_required_vec!(all_zero);
    merge_optional_vec!(disp_gene_est);
    merge_optional_vec!(disp_gene_iter);
    merge_optional_vec!(disp_fit);
    merge_optional_vec!(disp_map);
    merge_optional_vec!(dispersion);
    merge_optional_vec!(disp_iter);
    merge_optional_vec!(disp_outlier);
    merge_optional_vec!(dispersion_converged);
    merge_matrix!(beta);
    merge_matrix!(beta_se);
    merge_matrix!(beta_covariance);
    merge_optional_vec!(beta_converged);
    merge_optional_vec!(beta_iter);
    merge_optional_vec!(log_like);
    merge_optional_vec!(full_deviance);
    merge_matrix!(mu);
    // Keep the original Cook's distances: DESeq2 records those before replacement.
    Ok(())
}

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_estimate_size_factors(
    n_genes: usize,
    n_samples: usize,
    counts: *const f64,
    method: i32,
    geo_means: *const f64,
    control_genes: *const usize,
    n_control: usize,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts)? };
        let geo = if geo_means.is_null() {
            None
        } else {
            let log_means = unsafe { input_slice(geo_means, n_genes, "log geometric means")? };
            Some(
                log_means
                    .iter()
                    .map(|value| value.exp())
                    .collect::<Vec<_>>(),
            )
        };
        let controls = if control_genes.is_null() {
            None
        } else {
            Some(unsafe { input_slice(control_genes, n_control, "control genes")? })
        };
        let method = match method {
            0 => SizeFactorMethod::Ratio,
            1 => SizeFactorMethod::PosCounts,
            _ => return Err("unknown size-factor method".into()),
        };
        let factors = estimate_size_factors_with_options(&counts, method, geo.as_deref(), controls)
            .map_err(|e| e.to_string())?;
        unsafe { output_slice(output, n_samples, "size factors")? }.copy_from_slice(&factors);
        Ok(())
    })
}

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_normalized_counts(
    n_genes: usize,
    n_samples: usize,
    counts: *const f64,
    size_factors: *const f64,
    normalization_factors: *const f64,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts)? };
        let normalized = if normalization_factors.is_null() {
            let sf = unsafe { input_slice(size_factors, n_samples, "size factors")? };
            normalized_counts(&counts, sf).map_err(|e| e.to_string())?
        } else {
            let factors =
                unsafe { row_major_from_column_major(n_genes, n_samples, normalization_factors)? };
            normalized_counts_with_factors(&counts, &factors).map_err(|e| e.to_string())?
        };
        unsafe { copy_matrix_to_r(&normalized, output) }
    })
}

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_base_metadata(
    n_genes: usize,
    n_samples: usize,
    counts: *const f64,
    size_factors: *const f64,
    normalization_factors: *const f64,
    weights: *const f64,
    base_means: *mut f64,
    base_vars: *mut f64,
    all_zero: *mut u8,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts)? };
        let normalized = if normalization_factors.is_null() {
            let sf = unsafe { input_slice(size_factors, n_samples, "size factors")? };
            normalized_counts(&counts, sf).map_err(|e| e.to_string())?
        } else {
            let factors =
                unsafe { row_major_from_column_major(n_genes, n_samples, normalization_factors)? };
            normalized_counts_with_factors(&counts, &factors).map_err(|e| e.to_string())?
        };
        let (means, vars) = if weights.is_null() {
            (
                base_mean(&normalized).map_err(|e| e.to_string())?,
                base_variance(&normalized).map_err(|e| e.to_string())?,
            )
        } else {
            let weights = unsafe { row_major_from_column_major(n_genes, n_samples, weights)? };
            (
                base_mean_with_weights(&normalized, &weights).map_err(|e| e.to_string())?,
                base_variance_with_weights(&normalized, &weights).map_err(|e| e.to_string())?,
            )
        };
        unsafe { output_slice(base_means, n_genes, "base means")? }.copy_from_slice(&means);
        unsafe { output_slice(base_vars, n_genes, "base variances")? }.copy_from_slice(&vars);
        let zero = unsafe { output_slice(all_zero, n_genes, "all-zero flags")? };
        for (gene, flag) in zero.iter_mut().enumerate() {
            *flag = u8::from(counts.is_all_zero_gene(gene).map_err(|e| e.to_string())?);
        }
        Ok(())
    })
}

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
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
    replace_counts_output: *mut f64,
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
) -> i32 {
    run_backend(|| {
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
        let reduced = if selected_test_type == TestType::Lrt {
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
            builder = builder.reduced_design(reduced.clone());
            Some(reduced)
        } else {
            None
        };
        let (fit, results, replacement) = if min_replicates_for_replace == 0 {
            let result = if selected_test_type == TestType::Wald {
                builder.fit_with_results(&counts, &design)
            } else {
                builder.fit_lrt_with_results(&counts, &design, reduced.as_ref().unwrap())
            }
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
            let mut merged_fit = output.original_fit().clone();
            if let Some(refit_fit) = output.refit_fit() {
                merge_refit_fit_rows(&mut merged_fit, refit_fit, &output.refit_plan().refit_rows)?;
            }
            for &row in &output.refit_plan().new_all_zero_rows {
                merged_fit.base_mean[row] = 0.0;
                merged_fit.base_var[row] = 0.0;
                merged_fit.all_zero[row] = true;
            }
            (merged_fit, output.results().clone(), Some(output))
        };
        const COLUMNS: usize = 11;
        let output = unsafe {
            output_slice(
                result_output,
                n_genes
                    .checked_mul(COLUMNS)
                    .ok_or("result dimensions overflow")?,
                "result table",
            )?
        };
        for (gene, row) in results.rows.iter().enumerate() {
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
        let replace_counts = unsafe {
            output_slice(
                replace_counts_output,
                n_genes
                    .checked_mul(n_samples)
                    .ok_or("replacement count dimensions overflow")?,
                "replacement counts",
            )?
        };
        let replace = unsafe { output_slice(replace_output, n_genes, "replacement flags")? };
        let replaceable_samples =
            unsafe { output_slice(replaceable_samples_output, n_samples, "replaceable samples")? };
        let refit_rows = unsafe { output_slice(refit_rows_output, n_genes, "refit rows")? };
        let new_all_zero_rows =
            unsafe { output_slice(new_all_zero_rows_output, n_genes, "new all-zero rows")? };
        let n_refit = unsafe { output_slice(n_refit_output, 1, "replacement refit count")? };
        if let Some(output) = &replacement {
            let plan = output.refit_plan();
            let has_replaceable_samples = plan
                .replacement
                .replaceable_samples
                .iter()
                .any(|&value| value);
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
            if has_replaceable_samples {
                unsafe { copy_counts_to_r(&plan.replacement.replaced_counts, replace_counts)? };
            } else {
                unsafe { copy_counts_to_r(&counts, replace_counts)? };
            }
        } else {
            n_refit[0] = 0;
            replace.fill(2);
            replaceable_samples.fill(0);
            refit_rows.fill(0);
            new_all_zero_rows.fill(0);
            unsafe { copy_counts_to_r(&counts, replace_counts)? };
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
        Ok(())
    })
}

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

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_transform(
    n_genes: usize,
    n_samples: usize,
    counts_values: *const f64,
    size_factors_values: *const f64,
    normalization_factors_values: *const f64,
    observation_weights_values: *const f64,
    n_coefficients: usize,
    design_values: *const f64,
    fit_type: i32,
    transform_type: i32,
    blind: u8,
    nsub: usize,
    threads: usize,
    frozen_intercept_values: *const f64,
    sample_prior_variance: f64,
    frozen_trend_type: i32,
    frozen_trend_parameters: *const f64,
    frozen_disp_fit_values: *const f64,
    frozen_local_curve: *const f64,
    transformed_output: *mut f64,
    size_factors_output: *mut f64,
    intercept_output: *mut f64,
    prior_variance_output: *mut f64,
) -> i32 {
    run_backend(|| {
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts_values)? };
        let design_values =
            unsafe { row_major_from_column_major(n_samples, n_coefficients, design_values)? };
        let design = DesignMatrix::from_row_major(
            n_samples,
            n_coefficients,
            design_values.into_values(),
            None,
        )
        .map_err(|error| error.to_string())?;
        let selected_fit_type = match fit_type {
            0 => FitType::Parametric,
            1 => FitType::Local,
            2 => FitType::Mean,
            _ => return Err("fitType must be parametric, local, or mean".into()),
        };
        let mut builder = DeseqBuilder::new()
            .automatic_dispersion_means()
            .fit_type(selected_fit_type)
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

        let transformed = match transform_type {
            0 => {
                if nsub == 0 {
                    return Err("nsub must be a positive integer".into());
                }
                // DESeq2's `vst()` always fits the fast subset and errors if
                // fewer than `nsub` rows have ordinary normalized mean > 5.
                // The builder's `*_auto` helpers intentionally retain a
                // Rust-specific full-data fallback, so use the strict fast
                // path for the Bioconductor-compatible wrapper.
                let output = builder
                    .fast_vst_glm_mu(&counts, &design, nsub)
                    .map_err(|error| error.to_string())?;
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&output.subset_fit.size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?.fill(f64::NAN);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] = f64::NAN;
                }
                output.transformed
            }
            1 => {
                let (fit, rlog) = if frozen_intercept_values.is_null() {
                    if !frozen_disp_fit_values.is_null() {
                        let mut fit = builder
                            .fit_dispersion_trend_glm_mu(&counts, &design)
                            .map_err(|error| error.to_string())?;
                        fit.disp_fit = Some(unsafe {
                            input_slice(frozen_disp_fit_values, n_genes, "frozen rlog dispFit")?
                                .to_vec()
                        });
                        let rlog = fit
                            .regularized_log_transform(&counts)
                            .map_err(|error| error.to_string())?;
                        (fit, rlog)
                    } else {
                        let output = if blind != 0 {
                            builder.blind_rlog_glm_mu_with_fit(&counts)
                        } else {
                            builder.rlog_glm_mu_with_fit(&counts, &design)
                        }
                        .map_err(|error| error.to_string())?;
                        (output.fit, output.rlog)
                    }
                } else {
                    if !sample_prior_variance.is_finite() || sample_prior_variance <= 0.0 {
                        return Err("sample prior variance must be positive and finite with frozen intercepts".into());
                    }
                    let frozen_intercept = unsafe {
                        input_slice(frozen_intercept_values, n_genes, "frozen rlog intercepts")?
                    };
                    let fit = builder
                        .fit_dispersion_trend_glm_mu(&counts, &design)
                        .map_err(|error| error.to_string())?;
                    let rlog = fit
                        .frozen_rlog(&counts, frozen_intercept, sample_prior_variance)
                        .map_err(|error| error.to_string())?;
                    (fit, rlog)
                };
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&fit.size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?
                        .copy_from_slice(&rlog.intercept);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] =
                        rlog.sample_prior_variance;
                }
                rlog.transformed
            }
            2 => {
                let (transformed, output_size_factors) = if frozen_trend_type >= 0 {
                    let size_factors = if !normalization_factors_values.is_null() {
                        if size_factors_values.is_null() {
                            vec![1.0; n_samples]
                        } else {
                            unsafe {
                                input_slice(size_factors_values, n_samples, "size factors")?
                                    .to_vec()
                            }
                        }
                    } else if !size_factors_values.is_null() {
                        unsafe {
                            input_slice(size_factors_values, n_samples, "size factors")?.to_vec()
                        }
                    } else {
                        estimate_size_factors_with_options(
                            &counts,
                            SizeFactorMethod::Ratio,
                            None,
                            None,
                        )
                        .map_err(|error| error.to_string())?
                    };
                    let normalization_factors = if normalization_factors_values.is_null() {
                        None
                    } else {
                        Some(unsafe {
                            row_major_from_column_major(
                                n_genes,
                                n_samples,
                                normalization_factors_values,
                            )?
                        })
                    };
                    let normalized = match &normalization_factors {
                        Some(factors) => normalized_counts_with_factors(&counts, factors),
                        None => normalized_counts(&counts, &size_factors),
                    }
                    .map_err(|error| error.to_string())?;
                    let transformed =
                        match frozen_trend_type {
                            0 => {
                                let parameters = unsafe {
                                    input_slice(
                                        frozen_trend_parameters,
                                        2,
                                        "frozen parametric trend parameters",
                                    )?
                                };
                                vst_parametric(
                                    &normalized,
                                    ParametricDispersionTrend {
                                        extra_pois: parameters[0],
                                        asympt_disp: parameters[1],
                                    },
                                )
                            }
                            1 => {
                                let curve = unsafe {
                                    input_slice(frozen_local_curve, 1998, "frozen local VST curve")?
                                };
                                let inverse_size_factor_mean = match &normalization_factors {
                                Some(factors) => {
                                    local_vst_inverse_size_factor_mean_from_normalization_factors(
                                        factors,
                                    )
                                }
                                None => local_vst_inverse_size_factor_mean(&size_factors),
                            }
                            .map_err(|error| error.to_string())?;
                                vst_local_from_dispersion_grid(
                                    &normalized,
                                    &curve[..999],
                                    &curve[999..],
                                    inverse_size_factor_mean,
                                )
                            }
                            2 => {
                                let parameters = unsafe {
                                    input_slice(
                                        frozen_trend_parameters,
                                        2,
                                        "frozen mean trend parameters",
                                    )?
                                };
                                vst_mean(&normalized, parameters[0])
                            }
                            _ => return Err(
                                "frozen dispersion trend type must be parametric, local, or mean"
                                    .into(),
                            ),
                        }
                        .map_err(|error| error.to_string())?;
                    (transformed, size_factors)
                } else {
                    let fit = builder
                        .fit_dispersion_trend_glm_mu(&counts, &design)
                        .map_err(|error| error.to_string())?;
                    let transformed = fit.vst(&counts).map_err(|error| error.to_string())?;
                    (transformed, fit.size_factors)
                };
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&output_size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?.fill(f64::NAN);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] = f64::NAN;
                }
                transformed
            }
            _ => return Err("unknown transform type".into()),
        };
        unsafe { copy_matrix_to_r(&transformed, transformed_output) }
    })
}

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_cooks_mask(
    n_genes: usize,
    pvalues: *const f64,
    max_cooks: *const f64,
    cutoff: f64,
    counts: *const f64,
    n_samples: usize,
    cooks: *const f64,
    use_low_count_heuristic: u8,
    output_pvalues: *mut f64,
    output_outlier: *mut u8,
) -> i32 {
    run_backend(|| {
        let pvalues = unsafe { input_slice(pvalues, n_genes, "p-values")? };
        let max_cooks = unsafe { input_slice(max_cooks, n_genes, "maximum Cook's distances")? };
        let mut results = DeseqResults {
            rows: (0..n_genes)
                .map(|gene| DeseqResultRow {
                    gene: None,
                    base_mean: 0.0,
                    log2_fold_change: None,
                    lfc_se: None,
                    stat: None,
                    pvalue: pvalues[gene].is_finite().then_some(pvalues[gene]),
                    padj: None,
                    dispersion: None,
                    converged: None,
                    max_cooks: max_cooks[gene].is_finite().then_some(max_cooks[gene]),
                    cooks_outlier: None,
                    filtered: None,
                })
                .collect(),
            ..DeseqResults::default()
        };
        if use_low_count_heuristic != 0 {
            let count_matrix = unsafe { counts_from_r(n_genes, n_samples, counts)? };
            let cooks_matrix = unsafe { row_major_from_column_major(n_genes, n_samples, cooks)? };
            apply_cooks_cutoff_with_low_count_heuristic(
                &mut results,
                Some(cutoff),
                &count_matrix,
                &cooks_matrix,
            )
            .map_err(|error| error.to_string())?;
        } else {
            apply_cooks_cutoff(&mut results, Some(cutoff)).map_err(|error| error.to_string())?;
        }
        let out_p = unsafe { output_slice(output_pvalues, n_genes, "masked p-values")? };
        let outlier = unsafe { output_slice(output_outlier, n_genes, "Cook's outlier flags")? };
        for (gene, row) in results.rows.iter().enumerate() {
            out_p[gene] = row.pvalue.unwrap_or(f64::NAN);
            outlier[gene] = row.cooks_outlier.map_or(u8::MAX, u8::from);
        }
        Ok(())
    })
}

const DIAGNOSTIC_SCHEMA_NAMES: [&[u8]; 7] = [
    b"betaConv\0",
    b"fullBetaConv\0",
    b"reducedBetaConv\0",
    b"betaIter\0",
    b"reducedBetaIter\0",
    b"deviance\0",
    b"maxCooks\0",
];

#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_placeholder() -> i32 {
    0
}

/// Number of DESeq2-style diagnostic schema names exported by the R bridge.
#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_diagnostic_schema_len() -> usize {
    DIAGNOSTIC_SCHEMA_NAMES.len()
}

/// Return one null-terminated diagnostic schema name by index.
///
/// The returned pointer has static lifetime. A null pointer indicates an
/// out-of-range index.
#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_diagnostic_schema_name(index: usize) -> *const c_char {
    DIAGNOSTIC_SCHEMA_NAMES
        .get(index)
        .map_or(ptr::null(), |name| name.as_ptr().cast())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn diagnostic_schema_contract_is_stable() {
        let names = (0..rsdeseq2_r_diagnostic_schema_len())
            .map(|index| {
                let pointer = rsdeseq2_r_diagnostic_schema_name(index);
                assert!(!pointer.is_null());
                unsafe { CStr::from_ptr(pointer) }
                    .to_str()
                    .unwrap()
                    .to_string()
            })
            .collect::<Vec<_>>();

        assert_eq!(
            names,
            vec![
                "betaConv",
                "fullBetaConv",
                "reducedBetaConv",
                "betaIter",
                "reducedBetaIter",
                "deviance",
                "maxCooks"
            ]
        );
        assert!(rsdeseq2_r_diagnostic_schema_name(names.len()).is_null());
    }
}
