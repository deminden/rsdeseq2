#![deny(clippy::undocumented_unsafe_blocks)]

use std::cell::RefCell;
use std::ffi::{CString, c_char};
use std::mem::{MaybeUninit, size_of};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;

use rcompat_numeric::matrix::try_alloc_matrix;
use rsdeseq2::prelude::{
    Analysis, DeseqFit, GlmInput, MeanModel, Offsets, RunOutput, Stage,
    fit_glms_with_beta_prior, fit_glms_with_estimated_beta_prior,
    BetaPriorRefitOptions,
    BetaPriorVarianceOptions, CooksReplacementOptions, CountMatrix, DeseqBuilder, DesignMatrix,
    DispersionTrendFit, FitType, ParametricDispersionTrend, RowMajorMatrix, SizeFactorMethod,
    TestType, WaldAlternative, WaldTestOptions, base_mean, base_mean_with_weights, base_variance,
    base_variance_with_weights, estimate_beta_prior_variance, estimate_size_factors_with_options,
    local_vst_inverse_factor_mean,
    normalized_counts, normalized_counts_with_factors,
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
    // try_borrow: reading the message must not panic even while it is being set.
    guard_value(ptr::null(), || {
        LAST_ERROR.with(|slot| slot.try_borrow().map_or(ptr::null(), |s| s.as_ptr()))
    })
}

/// Run an entry point that returns a value instead of a status code. A panic
/// is caught at the FFI boundary, recorded as the last error and replaced by
/// `fallback`, so unwinding never crosses into C.
fn guard_value<T>(fallback: T, operation: impl FnOnce() -> T) -> T {
    use_r_profile();
    match catch_unwind(AssertUnwindSafe(operation)) {
        Ok(value) => value,
        Err(_) => {
            set_error("Rust backend panicked");
            fallback
        }
    }
}

/// The R package reproduces R's DESeq2 on the host, so it selects the host's R
/// profile once unless the user chose one through `RSDESEQ2_NUMERIC_PROFILE`.
fn use_r_profile() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        if std::env::var_os("RSDESEQ2_NUMERIC_PROFILE").is_none() {
            let _ = rsdeseq2::set_numeric_profile(rsdeseq2::NumericProfile::r_host());
        }
    });
}

fn run_backend(operation: impl FnOnce() -> Result<(), String>) -> i32 {
    use_r_profile();
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

fn check_buffer<T>(pointer: *const T, length: usize, name: &str, kind: &str) -> Result<(), String> {
    if pointer.is_null() {
        return Err(format!("{name} {kind} is null"));
    }
    if !pointer.is_aligned() {
        return Err(format!("{name} {kind} is misaligned"));
    }
    if length
        .checked_mul(size_of::<T>())
        .is_none_or(|bytes| bytes > isize::MAX as usize)
    {
        return Err(format!("{name} {kind} span exceeds isize::MAX bytes"));
    }
    Ok(())
}

/// # Safety
/// For nonempty buffers that pass the pointer and span checks, the caller
/// guarantees `length` readable, initialized elements in one allocation
/// that remains valid and immutable for `'a`.
unsafe fn input_slice<'a, T>(
    pointer: *const T,
    length: usize,
    name: &str,
) -> Result<&'a [T], String> {
    if length == 0 {
        return Ok(&[]);
    }
    check_buffer(pointer, length, name, "pointer")?;
    // SAFETY: the address and span were checked above, and the caller guarantees the storage and lifetime.
    Ok(unsafe { std::slice::from_raw_parts(pointer, length) })
}

/// # Safety
/// For nonempty buffers that pass the pointer and span checks, the caller
/// guarantees storage for `length` elements in one writable allocation,
/// exclusively borrowed for `'a`. The storage may be uninitialized.
unsafe fn output_slice<'a, T>(
    pointer: *mut T,
    length: usize,
    name: &str,
) -> Result<&'a mut [MaybeUninit<T>], String> {
    if length == 0 {
        return Ok(&mut []);
    }
    check_buffer(pointer, length, name, "output pointer")?;
    // SAFETY: the address and span were checked above, and the caller guarantees exclusive storage and lifetime.
    Ok(unsafe { std::slice::from_raw_parts_mut(pointer.cast(), length) })
}

/// # Safety
/// The caller meets the output buffer contract and every element is already
/// initialized as a valid `T`.
unsafe fn inout_slice<'a, T>(
    pointer: *mut T,
    length: usize,
    name: &str,
) -> Result<&'a mut [T], String> {
    // SAFETY: the caller guarantees exclusive storage and lifetime.
    let storage = unsafe { output_slice(pointer, length, name)? };
    // SAFETY: the caller guarantees every element is initialized.
    Ok(unsafe { storage.assume_init_mut() })
}

unsafe fn counts_from_r(
    n_genes: usize,
    n_samples: usize,
    values: *const f64,
) -> Result<CountMatrix, String> {
    // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
    let column_major = unsafe {
        input_slice(
            values,
            n_genes
                .checked_mul(n_samples)
                .ok_or("count dimensions overflow")?,
            "counts",
        )?
    };
    rsdeseq2::matrix::CountView::from_doubles(
        rsdeseq2::matrix::MatrixView::new(
            n_genes,
            n_samples,
            column_major,
            rsdeseq2::matrix::MatrixOrder::ColumnMajor,
        )
        .map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let mut row_major =
        try_alloc_matrix(n_genes, n_samples, "R counts").map_err(|error| error.to_string())?;
    for gene in 0..n_genes {
        for sample in 0..n_samples {
            row_major.push(column_major[gene + sample * n_genes] as u32);
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
    // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
    let column_major = unsafe {
        input_slice(
            values,
            n_rows
                .checked_mul(n_cols)
                .ok_or("matrix dimensions overflow")?,
            "matrix",
        )?
    };
    let mut row_major =
        try_alloc_matrix(n_rows, n_cols, "R matrix").map_err(|error| error.to_string())?;
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
    // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
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
            out[row + col * n_rows]
                .write(*matrix.get(row, col).ok_or("matrix index out of bounds")?);
        }
    }
    Ok(())
}

fn copy_counts_to_r(counts: &CountMatrix, output: &mut [MaybeUninit<f64>]) -> Result<(), String> {
    if output.len() != counts.n_genes().saturating_mul(counts.n_samples()) {
        return Err("replacement count output has an invalid length".into());
    }
    for gene in 0..counts.n_genes() {
        for sample in 0..counts.n_samples() {
            output[gene + sample * counts.n_genes()]
                .write(f64::from(counts.as_slice()[gene * counts.n_samples() + sample]));
        }
    }
    Ok(())
}

// Size factors, normalized counts and base metadata exports.
#[path = "lib/normalization.rs"]
mod normalization_exports;
#[allow(unused_imports)]
pub use normalization_exports::*;

// Wald/LRT fitting export, replacement and fitted output buffers.
#[path = "lib/fit.rs"]
mod fit_exports;
#[allow(unused_imports)]
pub use fit_exports::*;

// Coefficient shrinkage export.
#[path = "lib/shrinkage.rs"]
mod shrinkage_exports;
#[allow(unused_imports)]
pub use shrinkage_exports::*;

// Expanded-model contrast shrinkage export.
#[path = "lib/contrast_shrinkage.rs"]
mod contrast_shrinkage_exports;
#[allow(unused_imports)]
pub use contrast_shrinkage_exports::*;

// VST/rlog transformation export and frozen fit inputs.
#[path = "lib/transform.rs"]
mod transform_exports;
#[allow(unused_imports)]
pub use transform_exports::*;

// Standalone Wald/LRT probabilities and default Cook's cutoff.
#[path = "lib/inference_statistics.rs"]
mod inference_statistics;

// Standalone result adjustment and filtering; no R numerical callbacks.
#[path = "lib/result_filtering.rs"]
mod result_filtering_exports;
#[allow(unused_imports)]
pub use result_filtering_exports::*;

// Cook's filtering and stable diagnostic schema exports; test module is last.
#[path = "lib/diagnostics.rs"]
mod diagnostics_exports;
#[allow(unused_imports)]
pub use diagnostics_exports::*;

#[path = "lib/contrast_inference.rs"]
mod contrast_inference;

#[path = "lib/contrast_zero.rs"]
mod contrast_zero;

#[path = "lib/result_policy.rs"]
mod result_policy;

mod formula_request;
mod model_matrix;
mod row_preparation;

mod replacement_counts;

#[cfg(test)]
mod entry_point_tests {
    use super::check_buffer;
    use std::mem::size_of;
    use std::path::Path;
    use std::ptr;

    #[test]
    fn buffer_checks_reject_invalid_addresses_and_byte_spans() {
        let value = 0.0_f64;
        let pointer = &raw const value;
        assert!(check_buffer(pointer, 1, "input", "pointer").is_ok());
        assert_eq!(
            check_buffer(ptr::null::<f64>(), 1, "input", "pointer").unwrap_err(),
            "input pointer is null"
        );
        assert_eq!(
            check_buffer(pointer.wrapping_byte_add(1), 1, "input", "pointer").unwrap_err(),
            "input pointer is misaligned"
        );
        for length in [isize::MAX as usize / size_of::<f64>() + 1, usize::MAX] {
            assert_eq!(
                check_buffer(pointer, length, "result", "output pointer").unwrap_err(),
                "result output pointer span exceeds isize::MAX bytes"
            );
        }
    }

    /// Every exported `rsdeseq2_r_*` function must catch panics at the FFI
    /// boundary through `run_backend` or `guard_value`.
    #[test]
    fn every_exported_entry_point_is_panic_guarded() {
        fn sources(dir: &Path, out: &mut Vec<(String, String)>) {
            for entry in std::fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    sources(&path, out);
                } else if path.extension().is_some_and(|ext| ext == "rs") {
                    out.push((
                        path.display().to_string(),
                        std::fs::read_to_string(&path).unwrap(),
                    ));
                }
            }
        }
        let mut files = Vec::new();
        sources(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),
            &mut files,
        );
        let marker = "extern \"C\" fn rsdeseq2_r_";
        let (mut total, mut unguarded) = (0, Vec::new());
        for (path, text) in &files {
            let mut rest = text.as_str();
            while let Some(start) = rest.find(marker) {
                let after = &rest[start..];
                let name_end = after.find('(').unwrap();
                let name = &after["extern \"C\" fn ".len()..name_end];
                let open = after.find('{').unwrap();
                let mut depth = 0;
                let mut end = open;
                for (offset, ch) in after[open..].char_indices() {
                    match ch {
                        '{' => depth += 1,
                        '}' => {
                            depth -= 1;
                            if depth == 0 {
                                end = open + offset;
                                break;
                            }
                        }
                        _ => {}
                    }
                }
                let body = &after[open..=end];
                total += 1;
                if !body.contains("run_backend(") && !body.contains("guard_value(") {
                    unguarded.push(format!("{path}: {name}"));
                }
                rest = &after[end..];
            }
        }
        assert!(total >= 34, "found only {total} entry points");
        assert!(
            unguarded.is_empty(),
            "unguarded entry points: {unguarded:?}"
        );
    }
}
