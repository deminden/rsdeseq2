//! Synchronous borrowing of R's column-major inference arrays.
use super::contrast_zero::{ContrastZero, apply_zero_context};
use super::inference_statistics::borrowed_wald_options;
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::glm::{CovarianceView, InferenceEstimates, wald_test_contrast_from_estimates};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};

/// # Safety
/// C validates the complete G×P beta/optional SE and G×P×P covariance buffers,
/// P contrast values, scalar options, scalar/G df and G-element outputs.
/// Presence flags are explicit; NULL never selects statistical behavior.
/// An optional zero-context descriptor is validated in C, including every
/// pointed-to count/design/string allocation. No pointer is retained.
/// Input and output buffers must not overlap.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_contrast_inference(
    genes: usize,
    coefficients: usize,
    beta: *const f64,
    has_se: u8,
    se: *const f64,
    has_covariance: u8,
    covariance: *const f64,
    contrast: *const f64,
    use_wald: u8,
    threshold: f64,
    alternative: i32,
    use_t: u8,
    df_len: usize,
    df: *const f64,
    missing: f64,
    lfc: *mut f64,
    lfc_se: *mut f64,
    stat: *mut f64,
    pvalue: *mut f64,
    zero_context: *const ContrastZero,
) -> i32 {
    run_backend(|| {
        if has_se > 1 || has_covariance > 1 || use_wald > 1 {
            return Err("contrast presence flags must be zero or one".into());
        }
        let length = genes
            .checked_mul(coefficients)
            .ok_or("beta shape overflow")?;
        let matrix = |values| -> Result<_, String> {
            MatrixView::new(
                genes,
                coefficients,
                unsafe { input_slice(values, length, "contrast estimates")? },
                MatrixOrder::ColumnMajor,
            )
            .map_err(|e| e.to_string())
        };
        let beta = matrix(beta)?;
        let se = if has_se == 1 { Some(matrix(se)?) } else { None };
        let covariance = if has_covariance == 1 {
            let width = coefficients
                .checked_mul(coefficients)
                .ok_or("covariance width overflow")?;
            let length = genes
                .checked_mul(width)
                .ok_or("covariance shape overflow")?;
            let values = unsafe { input_slice(covariance, length, "contrast covariance")? };
            Some(
                CovarianceView::new(
                    coefficients,
                    MatrixView::new(genes, width, values, MatrixOrder::ColumnMajor)
                        .map_err(|e| e.to_string())?,
                )
                .map_err(|e| e.to_string())?,
            )
        } else {
            None
        };
        let estimates = InferenceEstimates::new(beta, se, covariance).map_err(|e| e.to_string())?;
        let contrast = unsafe { input_slice(contrast, coefficients, "contrast coefficients")? };
        let (mut effect, se, mut inference) = if use_wald == 1 {
            let options =
                unsafe { borrowed_wald_options(genes, threshold, alternative, use_t, df_len, df)? };
            let result = wald_test_contrast_from_estimates(&estimates, contrast, &options)
                .map_err(|e| e.to_string())?;
            (result.log2_fold_change, result.lfc_se, Some(result.wald))
        } else {
            let result = estimates
                .project_contrast(contrast)
                .map_err(|e| e.to_string())?;
            (result.log2_fold_change, result.lfc_se, None)
        };
        unsafe {
            apply_zero_context(
                zero_context,
                genes,
                coefficients,
                contrast,
                &mut effect,
                inference.as_mut(),
            )?;
        }
        let output_effect = unsafe { output_slice(lfc, genes, "contrast effects")? };
        let output_se = unsafe { output_slice(lfc_se, genes, "contrast SEs")? };
        for gene in 0..genes {
            output_effect[gene] = effect[gene].unwrap_or(missing);
            output_se[gene] = se[gene].unwrap_or(missing);
        }
        if let Some(result) = inference {
            let output_stat = unsafe { output_slice(stat, genes, "contrast statistics")? };
            let output_pvalue = unsafe { output_slice(pvalue, genes, "contrast p-values")? };
            for gene in 0..genes {
                output_stat[gene] = result.stat[gene].unwrap_or(missing);
                output_pvalue[gene] = result.pvalue[gene].unwrap_or(missing);
            }
        }
        Ok(())
    })
}
