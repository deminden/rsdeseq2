//! Decode a validated zero-count context and apply the standalone row policy.
use super::input_slice;
use rsdeseq2::contrasts::{ContrastCounts, ContrastZeroPolicy};
use rsdeseq2::glm::WaldOutput;
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use std::ffi::{CStr, c_char};

/// Matches contrast_zero.h. Every pointed-to allocation is checked in C.
#[repr(C)]
pub(super) struct ContrastZero {
    samples: usize,
    counts: *const f64,
    design: *const f64,
    has_design: u8,
    rule: i32,
    sample_levels: *const *const c_char,
    numerator: *const c_char,
    denominator: *const c_char,
}

// SAFETY: C validates the descriptor, complete matrix lengths and N UTF-8,
// NUL-terminated sample labels. No R call occurs while these borrows are live.
// Output vectors belong to the current Rust projection, not the borrowed fit.
pub(super) unsafe fn apply_zero_context(
    context: *const ContrastZero,
    genes: usize,
    coefficients: usize,
    contrast: &[f64],
    effect: &mut [Option<f64>],
    inference: Option<&mut WaldOutput>,
) -> Result<(), String> {
    let Some(context) = (unsafe { context.as_ref() }) else {
        return Ok(());
    };
    let length = genes
        .checked_mul(context.samples)
        .ok_or("contrast count shape overflow")?;
    let values = unsafe { input_slice(context.counts, length, "contrast counts")? };
    let matrix = MatrixView::new(genes, context.samples, values, MatrixOrder::ColumnMajor)
        .map_err(|e| e.to_string())?;
    let counts = ContrastCounts::from_doubles(matrix).map_err(|e| e.to_string())?;
    let flags = match context.rule {
        1 => {
            if context.has_design > 1 {
                return Err("contrast design flag must be zero or one".into());
            }
            let design = if context.has_design == 1 {
                let length = context
                    .samples
                    .checked_mul(coefficients)
                    .ok_or("contrast design shape overflow")?;
                let values = unsafe { input_slice(context.design, length, "contrast design")? };
                Some(
                    MatrixView::new(
                        context.samples,
                        coefficients,
                        values,
                        MatrixOrder::ColumnMajor,
                    )
                    .map_err(|e| e.to_string())?,
                )
            } else {
                None
            };
            counts.numeric_flags(design, contrast)
        }
        2 => {
            let pointers =
                unsafe { input_slice(context.sample_levels, context.samples, "sample levels")? };
            let text = |pointer: *const c_char| -> Result<&str, String> {
                if pointer.is_null() {
                    return Err("contrast level string is null".into());
                }
                unsafe { CStr::from_ptr(pointer) }
                    .to_str()
                    .map_err(|e| e.to_string())
            };
            let levels = pointers
                .iter()
                .map(|&p| text(p))
                .collect::<Result<Vec<_>, _>>()?;
            counts.factor_flags(
                &levels,
                text(context.numerator)?,
                text(context.denominator)?,
            )
        }
        _ => return Err("unknown contrast zero rule".into()),
    }
    .map_err(|e| e.to_string())?;
    let policy = ContrastZeroPolicy::new(&flags, counts.all_zero()).map_err(|e| e.to_string())?;
    if let Some(inference) = inference {
        policy.apply_wald(effect, &mut inference.stat, &mut inference.pvalue)
    } else {
        policy.apply_effects(effect.iter_mut())
    }
    .map_err(|e| e.to_string())
}
