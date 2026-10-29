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

// SAFETY: the C entry point must validate the allocated buffer's element count,
// type and lifetime before passing a pointer. Null checks cannot prove length.
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

// Merge fitted rows and diagnostics after Cook's replacement.
include!("lib/refit.rs");

// Size factors, normalized counts and base metadata exports.
include!("lib/normalization.rs");

// Wald/LRT fitting export, replacement and fitted output buffers.
include!("lib/fit.rs");

// Coefficient shrinkage export.
include!("lib/shrinkage.rs");

// Expanded-model contrast shrinkage export.
include!("lib/contrast_shrinkage.rs");

// VST/rlog transformation export and frozen fit inputs.
include!("lib/transform.rs");

// Cook's filtering and stable diagnostic schema exports.
include!("lib/diagnostics.rs");
