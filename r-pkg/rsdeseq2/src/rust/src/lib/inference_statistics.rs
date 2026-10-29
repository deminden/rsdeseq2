//! Borrowed primitive statistics. Every calculation is in the standalone core.
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::glm::{WaldAlternative, WaldTestOptions};

pub(super) unsafe fn borrowed_wald_options(
    genes: usize,
    threshold: f64,
    alternative: i32,
    use_t: u8,
    df_len: usize,
    df: *const f64,
) -> Result<WaldTestOptions, String> {
    let alternative = match alternative {
        0 => WaldAlternative::GreaterAbs,
        1 => WaldAlternative::GreaterAbsUpshot,
        2 => WaldAlternative::GreaterAbs2014,
        3 => WaldAlternative::LessAbs,
        4 => WaldAlternative::Greater,
        5 => WaldAlternative::Less,
        _ => return Err("unknown Wald alternative".into()),
    };
    if use_t > 1 {
        return Err("Wald useT flag must be zero or one".into());
    }
    let options = if use_t == 0 {
        if df_len != 0 || !df.is_null() {
            return Err("missing Wald df buffer".into());
        }
        WaldTestOptions::normal()
    } else {
        let df = unsafe { input_slice(df, df_len, "Wald df")? };
        if df.len() == 1 {
            WaldTestOptions::t_degrees_of_freedom(df[0])
        } else if df.len() == genes {
            WaldTestOptions::t_per_gene_degrees_of_freedom(df.to_vec())
        } else {
            return Err("Wald df must be scalar or one per gene".into());
        }
    }
    .with_lfc_threshold(threshold, alternative);
    Ok(options)
}

/// # Safety
/// C validates G-element beta/SE/output arrays, scalar or G-element df, and
/// scalar options. Buffers are borrowed only for this synchronous call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_wald_statistics(
    genes: usize,
    beta: *const f64,
    se: *const f64,
    threshold: f64,
    alternative: i32,
    use_t: u8,
    df_len: usize,
    df: *const f64,
    missing: f64,
    stat: *mut f64,
    pvalue: *mut f64,
) -> i32 {
    run_backend(|| {
        let options =
            unsafe { borrowed_wald_options(genes, threshold, alternative, use_t, df_len, df)? };
        let beta = unsafe { input_slice(beta, genes, "Wald estimates")? };
        let se = unsafe { input_slice(se, genes, "Wald standard errors")? };
        let result = rsdeseq2::glm::wald_test_from_estimates(beta, se, &options)
            .map_err(|e| e.to_string())?;
        let stat = unsafe { output_slice(stat, genes, "Wald statistics")? };
        let pvalue = unsafe { output_slice(pvalue, genes, "Wald p-values")? };
        for gene in 0..genes {
            stat[gene] = result.stat[gene].unwrap_or(missing);
            pvalue[gene] = result.pvalue[gene].unwrap_or(missing);
        }
        Ok(())
    })
}

/// # Safety
/// C validates G-element statistics/output and scalar or G-element df. No
/// pointer is retained, and the missing sentinel is transport formatting only.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_lrt_pvalues(
    genes: usize,
    statistics: *const f64,
    df_len: usize,
    df: *const f64,
    missing: f64,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        let statistics = unsafe { input_slice(statistics, genes, "LRT statistics")? };
        let df = unsafe { input_slice(df, df_len, "LRT df")? };
        let result = rsdeseq2::glm::lrt_pvalues_from_statistics(statistics, df)
            .map_err(|e| e.to_string())?;
        let output = unsafe { output_slice(output, genes, "LRT p-values")? };
        for (target, value) in output.iter_mut().zip(result) {
            *target = value.unwrap_or(missing);
        }
        Ok(())
    })
}

/// # Safety
/// The caller supplies one writable double and one writable flag. The flag
/// distinguishes an absent saturated-model cutoff from a numerical value.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_default_cooks_cutoff(
    samples: usize,
    coefficients: usize,
    output: *mut f64,
    present: *mut u8,
) -> i32 {
    run_backend(|| {
        let cutoff = rsdeseq2::results::default_cooks_cutoff(samples, coefficients)
            .map_err(|e| e.to_string())?;
        let output = unsafe { output_slice(output, 1, "Cook's cutoff")? };
        let present = unsafe { output_slice(present, 1, "Cook's cutoff flag")? };
        output[0] = cutoff.unwrap_or(0.0);
        present[0] = u8::from(cutoff.is_some());
        Ok(())
    })
}
