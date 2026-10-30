//! Bounded synchronous inference; both single and batch calls use this path.
use super::contrast_zero::{ContrastZero, decode_zero, prepare_zero_inputs};
use super::inference_statistics::borrowed_wald_options;
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::glm::{CovarianceView, InferenceEstimates, PreparedContrastView};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};

/// Matches contrast_inference.h. C validates every pointed-to buffer.
#[repr(C)]
pub(super) struct ContrastRequest {
    contrast: *const f64,
    threshold: f64,
    alternative: i32,
    use_wald: u8,
    use_t: u8,
    df_len: usize,
    df: *const f64,
    zero: *const ContrastZero,
    lfc: *mut f64,
    lfc_se: *mut f64,
    stat: *mut f64,
    pvalue: *mut f64,
}

/// # Safety
/// C validates G×P beta/optional SE, G×P×P covariance, 1..=64 request
/// descriptors, each P-element contrast, scalar/G df, all zero contexts and
/// distinct caller-owned G-element output arrays. No pointer is retained or
/// R API called. Outputs must not alias inputs or one another.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_contrast_batch(
    genes: usize,
    coefficients: usize,
    beta: *const f64,
    has_se: u8,
    se: *const f64,
    has_covariance: u8,
    covariance: *const f64,
    count: usize,
    requests: *const ContrastRequest,
    missing: f64,
) -> i32 {
    run_backend(|| {
        if has_se > 1 || has_covariance > 1 || !(1..=64).contains(&count) {
            return Err("invalid contrast batch flags or request count".into());
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
            Some(
                CovarianceView::new(
                    coefficients,
                    MatrixView::new(
                        genes,
                        width,
                        unsafe { input_slice(covariance, length, "contrast covariance")? },
                        MatrixOrder::ColumnMajor,
                    )
                    .map_err(|e| e.to_string())?,
                )
                .map_err(|e| e.to_string())?,
            )
        } else {
            None
        };
        let estimates = InferenceEstimates::new(beta, se, covariance).map_err(|e| e.to_string())?;
        let requests = unsafe { input_slice(requests, count, "contrast requests")? };
        let (counts, design) =
            unsafe { prepare_zero_inputs(requests.iter().map(|r| r.zero), genes, coefficients)? };
        let view =
            PreparedContrastView::new(estimates, counts, design).map_err(|e| e.to_string())?;
        for (index, request) in requests.iter().enumerate() {
            let evaluate = || -> Result<(), String> {
                if request.use_wald > 1 {
                    return Err("Wald mode must be zero or one".into());
                }
                let options = if request.use_wald == 1 {
                    Some(unsafe {
                        borrowed_wald_options(
                            genes,
                            request.threshold,
                            request.alternative,
                            request.use_t,
                            request.df_len,
                            request.df,
                        )?
                    })
                } else {
                    None
                };
                let contrast = unsafe {
                    input_slice(request.contrast, coefficients, "contrast coefficients")?
                };
                let zero = unsafe { decode_zero(request.zero)? };
                let result = view
                    .evaluate(contrast, options.as_ref(), zero.rule())
                    .map_err(|e| e.to_string())?;
                let effect = unsafe { output_slice(request.lfc, genes, "contrast effects")? };
                let se = unsafe { output_slice(request.lfc_se, genes, "contrast SEs")? };
                for gene in 0..genes {
                    effect[gene] = result.projection.log2_fold_change[gene].unwrap_or(missing);
                    se[gene] = result.projection.lfc_se[gene].unwrap_or(missing);
                }
                if let Some(wald) = result.wald {
                    let stat = unsafe { output_slice(request.stat, genes, "contrast statistics")? };
                    let pvalue =
                        unsafe { output_slice(request.pvalue, genes, "contrast p-values")? };
                    for gene in 0..genes {
                        stat[gene] = wald.stat[gene].unwrap_or(missing);
                        pvalue[gene] = wald.pvalue[gene].unwrap_or(missing);
                    }
                }
                Ok(())
            };
            evaluate().map_err(|e| format!("contrast request {}: {e}", index + 1))?;
        }
        Ok(())
    })
}
