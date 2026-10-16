// Cook's filtering and stable diagnostic schema exports.

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
