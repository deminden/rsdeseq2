// Borrowed Cook's filtering and stable diagnostic schema exports.
fn finite_value(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

unsafe fn borrowed_cooks_heuristic<'a>(
    genes: usize,
    samples: usize,
    counts: *const f64,
    cooks: *const f64,
) -> Result<rsdeseq2::results::LowCountHeuristic<'a>, String> {
    use rsdeseq2::matrix::{CountView, MatrixOrder, MatrixView};
    let length = genes.checked_mul(samples).ok_or("Cook's shape overflow")?;
    let matrix = |pointer, name| -> Result<_, String> {
        MatrixView::new(
            genes,
            samples,
            unsafe { input_slice(pointer, length, name)? },
            MatrixOrder::ColumnMajor,
        )
        .map_err(|e| e.to_string())
    };
    let counts =
        CountView::from_doubles(matrix(counts, "Cook's counts")?).map_err(|e| e.to_string())?;
    rsdeseq2::results::LowCountHeuristic::new(counts, matrix(cooks, "Cook's distances")?)
        .map_err(|e| e.to_string())
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
        if use_low_count_heuristic > 1 {
            return Err("low-count flag must be zero or one".into());
        }
        let heuristic = if use_low_count_heuristic == 1 {
            Some(unsafe { borrowed_cooks_heuristic(n_genes, n_samples, counts, cooks)? })
        } else {
            None
        };
        let policy = rsdeseq2::results::CooksPolicy::new(Some(cutoff), heuristic)
            .map_err(|e| e.to_string())?;
        let mask = policy
            .evaluate(
                pvalues.iter().copied().map(finite_value),
                max_cooks.iter().copied().map(finite_value),
            )
            .map_err(|e| e.to_string())?;
        let out_p = unsafe { output_slice(output_pvalues, n_genes, "masked p-values")? };
        let outlier = unsafe { output_slice(output_outlier, n_genes, "Cook's outlier flags")? };
        for gene in 0..n_genes {
            out_p[gene] = mask.pvalue[gene].unwrap_or(f64::NAN);
            outlier[gene] = mask.outlier[gene].map_or(u8::MAX, u8::from);
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
