//! Borrow validated zero metadata and prepare common count/design inputs once.
use super::input_slice;
use rsdeseq2::contrasts::ContrastCounts;
use rsdeseq2::glm::ContrastZeroRule;
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use std::ffi::{CStr, c_char};

/// Matches contrast_zero.h; every pointed-to allocation is checked in C.
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

pub(super) type ZeroInputs<'a> = (Option<ContrastCounts<'a>>, Option<MatrixView<'a, f64>>);

// SAFETY: descriptors/buffers/UTF-8 labels are validated in C. The returned
// borrows are used only during this synchronous call, with no R API callbacks.
pub(super) unsafe fn prepare_zero_inputs<'a>(
    contexts: impl Iterator<Item = *const ContrastZero>,
    genes: usize,
    coefficients: usize,
) -> Result<ZeroInputs<'a>, String> {
    let mut common = None;
    let mut numeric_design = None;
    for pointer in contexts {
        let Some(context) = (unsafe { pointer.as_ref() }) else {
            continue;
        };
        if !(1..=2).contains(&context.rule) || context.has_design > 1 {
            return Err("invalid contrast zero rule or design flag".into());
        }
        let identity = (context.counts, context.samples);
        if common.is_some_and(|old| old != identity) {
            return Err("batch zero contexts must share one count matrix".into());
        }
        common = Some(identity);
        if context.rule == 1 {
            let identity = (context.has_design, context.design);
            if numeric_design.is_some_and(|old| old != identity) {
                return Err("batch numeric zero contexts must share one design".into());
            }
            numeric_design = Some(identity);
        }
    }
    let Some((pointer, samples)) = common else {
        return Ok((None, None));
    };
    let length = genes
        .checked_mul(samples)
        .ok_or("contrast count shape overflow")?;
    let values = unsafe { input_slice(pointer, length, "contrast counts")? };
    let counts = ContrastCounts::from_doubles(
        MatrixView::new(genes, samples, values, MatrixOrder::ColumnMajor)
            .map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let design = if let Some((1, pointer)) = numeric_design {
        let length = samples
            .checked_mul(coefficients)
            .ok_or("contrast design shape overflow")?;
        let values = unsafe { input_slice(pointer, length, "contrast design")? };
        Some(
            MatrixView::new(samples, coefficients, values, MatrixOrder::ColumnMajor)
                .map_err(|e| e.to_string())?,
        )
    } else {
        None
    };
    Ok((Some(counts), design))
}

pub(super) enum DecodedZero<'a> {
    None,
    Numeric,
    Factor {
        levels: Vec<&'a str>,
        numerator: &'a str,
        denominator: &'a str,
    },
}
impl DecodedZero<'_> {
    pub(super) fn rule(&self) -> ContrastZeroRule<'_> {
        match self {
            Self::None => ContrastZeroRule::None,
            Self::Numeric => ContrastZeroRule::Numeric,
            Self::Factor {
                levels,
                numerator,
                denominator,
            } => ContrastZeroRule::Factor {
                levels,
                numerator,
                denominator,
            },
        }
    }
}

pub(super) unsafe fn decode_zero<'a>(
    pointer: *const ContrastZero,
) -> Result<DecodedZero<'a>, String> {
    let Some(context) = (unsafe { pointer.as_ref() }) else {
        return Ok(DecodedZero::None);
    };
    match context.rule {
        1 => Ok(DecodedZero::Numeric),
        2 => {
            let pointers =
                unsafe { input_slice(context.sample_levels, context.samples, "sample levels")? };
            let text = |pointer: *const c_char| -> Result<&'a str, String> {
                if pointer.is_null() {
                    return Err("contrast level string is null".into());
                }
                unsafe { CStr::from_ptr(pointer) }
                    .to_str()
                    .map_err(|e| e.to_string())
            };
            Ok(DecodedZero::Factor {
                levels: pointers
                    .iter()
                    .map(|&p| text(p))
                    .collect::<Result<Vec<_>, _>>()?,
                numerator: text(context.numerator)?,
                denominator: text(context.denominator)?,
            })
        }
        _ => Err("unknown contrast zero rule".into()),
    }
}
