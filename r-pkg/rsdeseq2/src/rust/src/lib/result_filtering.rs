// No R numerical callback: these are the standalone library's entry points.
fn filtering_pvalues(values: &[f64]) -> Vec<Option<f64>> {
    values
        .iter()
        .map(|&value| (!value.is_nan()).then_some(value))
        .collect()
}

/// # Safety
/// C validates both G-element buffers. The missing sentinel is output formatting
/// only; all BH arithmetic is performed by the standalone Rust library.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_bh_adjust(
    genes: usize,
    values: *const f64,
    missing: f64,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        let values = unsafe { input_slice(values, genes, "BH p-values")? };
        let output = unsafe { output_slice(output, genes, "BH output")? };
        let adjusted = rsdeseq2::multiple_testing::bh_adjust(&filtering_pvalues(values));
        for ((target, adjusted), &original) in output.iter_mut().zip(adjusted).zip(values) {
            *target = adjusted.unwrap_or(if original.is_nan() { original } else { missing });
        }
        Ok(())
    })
}

/// # Safety
/// C validates G-element inputs/outputs; metadata buffers have K elements where
/// K is zero when disabled, 50 for default theta, or the supplied theta length.
/// Scalar buffers have lengths 1 and 2. No pointer is retained and no R code runs.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_independent_filtering(
    genes: usize,
    filter: *const f64,
    values: *const f64,
    theta_len: usize,
    theta: *const f64,
    enabled: u8,
    alpha: f64,
    missing: f64,
    padj: *mut f64,
    filtered: *mut u8,
    theta_output: *mut f64,
    rejections: *mut i32,
    smooth: *mut f64,
    selected: *mut i32,
    scalars: *mut f64,
) -> i32 {
    run_backend(|| {
        let filter = unsafe { input_slice(filter, genes, "filter")? };
        let values = unsafe { input_slice(values, genes, "filter p-values")? };
        if enabled > 1 {
            return Err("filter enabled flag must be zero or one".into());
        }
        let theta = if theta.is_null() {
            None
        } else {
            Some(unsafe { input_slice(theta, theta_len, "filter theta")? }.to_vec())
        };
        let options = rsdeseq2::independent_filtering::IndependentFilteringOptions {
            enabled: enabled == 1,
            alpha,
            theta,
        };
        let result = rsdeseq2::independent_filtering::evaluate_independent_filtering(
            filter,
            &filtering_pvalues(values),
            &options,
        )
        .map_err(|e| e.to_string())?;
        let padj = unsafe { output_slice(padj, genes, "filter adjusted p-values")? };
        let filtered = unsafe { output_slice(filtered, genes, "filter flags")? };
        let length = if !options.enabled {
            0
        } else {
            options.theta.as_ref().map_or(50, Vec::len)
        };
        let theta_output = unsafe { output_slice(theta_output, length, "filter theta output")? };
        let rejections = unsafe { output_slice(rejections, length, "filter rejections")? };
        let smooth = unsafe { output_slice(smooth, length, "filter smoothing")? };
        let selected = unsafe { output_slice(selected, 1, "selected filter")? };
        let scalars = unsafe { output_slice(scalars, 2, "filter metadata")? };
        for gene in 0..genes {
            // Preserve NA/NaN payloads only as R transport formatting. Population
            // selection and adjustment have already been decided by the core.
            padj[gene] = rsdeseq2::independent_filtering::format_adjusted_pvalue(
                result.padj[gene], values[gene], filter[gene],
                result.metadata.filter_threshold, missing,
            );
            filtered[gene] = result.filtered[gene].map_or(u8::MAX, u8::from);
        }
        theta_output.copy_from_slice(&result.metadata.theta);
        for (target, &count) in rejections.iter_mut().zip(&result.metadata.num_rejections) {
            *target =
                i32::try_from(count).map_err(|_| "filter rejection count exceeds integer range")?;
        }
        if let Some(curve) = result.metadata.lowess_fit {
            smooth.copy_from_slice(&curve);
        }
        selected[0] = result.metadata.selected_index.map_or(Ok(-1), |index| {
            i32::try_from(index + 1).map_err(|_| "filter index exceeds integer range")
        })?;
        scalars[0] = result.metadata.filter_theta.unwrap_or(missing);
        scalars[1] = result.metadata.filter_threshold.unwrap_or(missing);
        Ok(())
    })
}
