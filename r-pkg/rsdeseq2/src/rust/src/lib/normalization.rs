// Size factors, normalized counts and base metadata exports.

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
