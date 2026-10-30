//! Validated synchronous model-matrix transport, with no R evaluator callbacks.
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::design::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use std::ffi::{CStr, c_char};

/// # Safety
/// Both index pointers reference their stated number of readable i32 values.
/// The distinct output references `columns` writable doubles. C validates the
/// R allocations; Rust validates all one-based indices before shared arithmetic.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_average_contrast(
    columns: usize,
    n_positive: usize,
    positive: *const i32,
    n_negative: usize,
    negative: *const i32,
    output: *mut f64,
) -> i32 {
    run_backend(|| unsafe {
        let indices = |p, n| -> Result<Vec<usize>, String> {
            input_slice(p, n, "average contrast indices")?
                .iter()
                .map(|&i| {
                    usize::try_from(i)
                        .ok()
                        .and_then(|i| i.checked_sub(1))
                        .ok_or_else(|| "average contrast indices must be positive integers".into())
                })
                .collect()
        };
        let values = rsdeseq2::contrasts::averaged_list_contrast(
            columns,
            &indices(positive, n_positive)?,
            &indices(negative, n_negative)?,
        )
        .map_err(|e| e.to_string())?;
        output_slice(output, columns, "average contrast values")?.copy_from_slice(&values);
        Ok(())
    })
}

/// C validates every pointed-to R allocation before invoking Rust.
#[repr(C)]
pub struct InputVariable {
    name: *const c_char,
    kind: i32,
    coding: i32,
    width: usize,
    values: *const f64,
    codes: *const i32,
    levels: *const *const c_char,
    names: *const *const c_char,
    contrast_columns: usize,
    contrast_names: *const *const c_char,
}
unsafe fn string(pointer: *const c_char) -> Result<String, String> {
    if pointer.is_null() {
        return Err("null model variable name/label".into());
    }
    unsafe { CStr::from_ptr(pointer) }
        .to_str()
        .map(str::to_owned)
        .map_err(|e| e.to_string())
}
unsafe fn strings(pointer: *const *const c_char, length: usize) -> Result<Vec<String>, String> {
    unsafe { input_slice(pointer, length, "model variable labels")? }
        .iter()
        .map(|&s| unsafe { string(s) })
        .collect()
}
struct Metadata {
    name: String,
    levels: Vec<String>,
    names: Option<Vec<String>>,
    contrast_names: Option<Vec<String>>,
}
unsafe fn with_plan<T>(
    samples: usize,
    inputs: &[InputVariable],
    n_terms: usize,
    codes: *const i32,
    controls: (i32, usize),
    operation: impl FnOnce(ModelMatrixPlan<'_>) -> Result<T, String>,
) -> Result<T, String> {
    let (intercept, response) = controls;
    if !(0..=1).contains(&intercept) {
        return Err("invalid model intercept flag".into());
    }
    let mut metadata = Vec::with_capacity(inputs.len());
    for input in inputs {
        metadata.push(Metadata {
            name: unsafe { string(input.name)? },
            levels: if input.kind == 1 {
                unsafe { strings(input.levels, input.width)? }
            } else {
                Vec::new()
            },
            names: if input.names.is_null() {
                None
            } else {
                Some(unsafe { strings(input.names, input.width)? })
            },
            contrast_names: if input.contrast_names.is_null() {
                None
            } else {
                Some(unsafe { strings(input.contrast_names, input.contrast_columns)? })
            },
        });
    }
    let view =
        |rows: usize, columns: usize, pointer: *const f64| -> Result<MatrixView<'_, f64>, String> {
            let size = rows
                .checked_mul(columns)
                .ok_or("model variable shape overflow")?;
            MatrixView::new(
                rows,
                columns,
                unsafe { input_slice(pointer, size, "model variable values")? },
                MatrixOrder::ColumnMajor,
            )
            .map_err(|e| e.to_string())
        };
    let mut variables = Vec::with_capacity(inputs.len());
    for (input, meta) in inputs.iter().zip(&metadata) {
        variables.push(match input.kind {
            0 => ModelMatrixVariable::Numeric {
                name: &meta.name,
                values: view(samples, input.width, input.values)?,
                names: meta.names.as_deref(),
            },
            1 => {
                let coding = match input.coding {
                    0 => FactorCoding::Treatment { base: 0 },
                    1 => FactorCoding::Sum,
                    2 => FactorCoding::Helmert,
                    3 => FactorCoding::Polynomial,
                    4 => FactorCoding::Treatment {
                        base: input.width.checked_sub(1).ok_or("empty factor levels")?,
                    },
                    5 => FactorCoding::Matrix {
                        values: view(input.width, input.contrast_columns, input.values)?,
                        names: meta.contrast_names.as_deref(),
                    },
                    _ => return Err("unknown factor coding".into()),
                };
                ModelMatrixVariable::Factor {
                    name: &meta.name,
                    codes: unsafe { input_slice(input.codes, samples, "factor codes")? },
                    levels: &meta.levels,
                    coding,
                }
            }
            _ => return Err("unknown model variable type".into()),
        });
    }
    if inputs.is_empty() && n_terms != 0 {
        return Err("model terms require variables".into());
    }
    let length = inputs
        .len()
        .checked_mul(n_terms)
        .ok_or("model term shape overflow")?;
    let term_codes = unsafe { input_slice(codes, length, "model term codes")? };
    let terms = if inputs.is_empty() {
        Vec::new()
    } else {
        term_codes
            .chunks_exact(inputs.len())
            .map(<[i32]>::to_vec)
            .collect()
    };
    operation(
        ModelMatrixPlan::with_response(
            samples,
            &variables,
            &terms,
            intercept != 0,
            response.checked_sub(1),
        )
        .map_err(|e| e.to_string())?,
    )
}
/// # Safety
/// Input variable pointers reference storage validated by C for the stated
/// dimensions; term codes reference V×T integers. Output pointers reference two
/// distinct writable usize values and one i32. No pointer is retained.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_model_matrix_shape(
    samples: usize,
    n_variables: usize,
    variables: *const InputVariable,
    n_terms: usize,
    codes: *const i32,
    intercept: i32,
    response: usize,
    n_columns: *mut usize,
    name_bytes: *mut usize,
    response_dropped: *mut i32,
) -> i32 {
    run_backend(|| unsafe {
        with_plan(
            samples,
            input_slice(variables, n_variables, "model variables")?,
            n_terms,
            codes,
            (intercept, response),
            |plan| {
                let bytes = plan.column_names().iter().try_fold(0usize, |n, s| {
                    n.checked_add(s.len() + 1).ok_or("model names overflow")
                })?;
                output_slice(n_columns, 1, "model column count")?[0] = plan.n_cols();
                output_slice(name_bytes, 1, "model name size")?[0] = bytes;
                output_slice(response_dropped, 1, "response warning flag")?[0] =
                    i32::from(plan.response_term_dropped());
                Ok(())
            },
        )
    })
}
/// # Safety
/// Inputs follow the shape-call contract. Distinct outputs reference N×K doubles,
/// K integers and name_capacity writable bytes. C allocates those buffers from
/// the validated shape result and performs no R call while Rust borrows them.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_model_matrix_fill(
    samples: usize,
    n_variables: usize,
    variables: *const InputVariable,
    n_terms: usize,
    codes: *const i32,
    intercept: i32,
    response: usize,
    n_columns: usize,
    output: *mut f64,
    assignments: *mut i32,
    names: *mut u8,
    name_capacity: usize,
) -> i32 {
    run_backend(|| unsafe {
        with_plan(
            samples,
            input_slice(variables, n_variables, "model variables")?,
            n_terms,
            codes,
            (intercept, response),
            |plan| {
                if plan.n_cols() != n_columns {
                    return Err("model output shape changed".into());
                }
                let out = output_slice(
                    output,
                    samples
                        .checked_mul(n_columns)
                        .ok_or("model output overflow")?,
                    "model matrix",
                )?;
                let values = plan.values();
                for column in 0..n_columns {
                    for sample in 0..samples {
                        out[sample + column * samples] = values[sample * n_columns + column];
                    }
                }
                let out = output_slice(assignments, n_columns, "model assignments")?;
                for (out, &value) in out.iter_mut().zip(plan.assignments()) {
                    *out = i32::try_from(value).map_err(|e| e.to_string())?;
                }
                let out = output_slice(names, name_capacity, "model names")?;
                let mut used = 0usize;
                for name in plan.column_names() {
                    let end = used
                        .checked_add(name.len() + 1)
                        .ok_or("model names overflow")?;
                    if end > out.len() {
                        return Err("model name buffer too small".into());
                    }
                    out[used..end - 1].copy_from_slice(name.as_bytes());
                    out[end - 1] = 0;
                    used = end;
                }
                Ok(())
            },
        )
    })
}
