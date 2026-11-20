// Size factors, normalized counts and base metadata exports.
// Optional geometric means use natural-scale units, as in the standalone API.
use super::*;
use rsdeseq2::normalization::{
    estimate_size_factors_from_view, visit_base_metadata, visit_normalized_counts_column_major,
};

/// # Safety
/// The selected input pointer is valid for its complete buffer and lifetime.
unsafe fn read_count_factors<'a>(
    genes: usize,
    samples: usize,
    size_factors: *const f64,
    normalization_factors: *const f64,
) -> Result<(&'a [f64], Option<MatrixView<'a, f64>>), String> {
    if normalization_factors.is_null() {
        // SAFETY: C provides one size factor per sample.
        let size_factors = unsafe { input_slice(size_factors, samples, "size factors")? };
        Ok((size_factors, None))
    } else {
        // SAFETY: C provides a normalization factor for each count cell.
        let factors = unsafe {
            column_major_view(genes, samples, normalization_factors, "normalization factors")?
        };
        Ok((&[], Some(factors)))
    }
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
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        let values = unsafe { column_major_view(n_genes, n_samples, counts, "counts")? };
        let counts = CountView::from_doubles(values).map_err(|error| error.to_string())?;
        let geo = if geo_means.is_null() {
            None
        } else {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            Some(unsafe { input_slice(geo_means, n_genes, "geometric means")? })
        };
        let controls = if control_genes.is_null() {
            None
        } else {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            Some(unsafe { input_slice(control_genes, n_control, "control genes")? })
        };
        let method = match method {
            0 => SizeFactorMethod::Ratio,
            1 => SizeFactorMethod::PosCounts,
            _ => return Err("unknown size-factor method".into()),
        };
        let factors = estimate_size_factors_from_view(counts, method, geo, controls)
            .map_err(|e| e.to_string())?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        unsafe { output_slice(output, n_samples, "size factors")? }.write_copy_of_slice(&factors);
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
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        let values = unsafe { column_major_view(n_genes, n_samples, counts, "counts")? };
        let counts = CountView::from_doubles(values).map_err(|error| error.to_string())?;
        // SAFETY: C validates the supplied factor dimensions and keeps the storage live.
        let (size_factors, factors) = unsafe {
            read_count_factors(n_genes, n_samples, size_factors, normalization_factors)?
        };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let output = unsafe { output_slice(output, values.as_slice().len(), "matrix")? };
        visit_normalized_counts_column_major(counts, size_factors, factors, |index, value| {
            output[index].write(value);
        })
        .map_err(|error| error.to_string())
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
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        let values = unsafe { column_major_view(n_genes, n_samples, counts, "counts")? };
        let counts = CountView::from_doubles(values).map_err(|error| error.to_string())?;
        // SAFETY: C validates the supplied factor dimensions and keeps the storage live.
        let (size_factors, factors) = unsafe {
            read_count_factors(n_genes, n_samples, size_factors, normalization_factors)?
        };
        let weights = if weights.is_null() {
            None
        } else {
            // SAFETY: C supplies one weight per count cell in a readable buffer.
            Some(unsafe { column_major_view(n_genes, n_samples, weights, "weights")? })
        };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let means = unsafe { output_slice(base_means, n_genes, "base means")? };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let vars = unsafe { output_slice(base_vars, n_genes, "base variances")? };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let zero = unsafe { output_slice(all_zero, n_genes, "all-zero flags")? };
        visit_base_metadata(counts, size_factors, factors, weights, |gene, mean, variance, all_zero| {
            means[gene].write(mean);
            vars[gene].write(variance);
            zero[gene].write(u8::from(all_zero));
        })
        .map_err(|error| error.to_string())
    })
}
