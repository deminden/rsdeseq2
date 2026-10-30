//! Primitive request adapter. Numerical implementations remain in the core.
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::design::expression::{
    self, Argument, Data, EvaluationOptions, Expression, Value, frame, wire,
};
use std::ptr;
pub struct FormulaOutput {
    bytes: Vec<u8>,
}
fn fields(v: &Value) -> Result<(&[String], &[Value]), String> {
    let Data::List(values) = &v.data else {
        return Err("formula request must be a named list".into());
    };
    let names = v.names.as_deref().ok_or("formula request requires names")?;
    if names.len() != values.len()
        || names
            .iter()
            .enumerate()
            .any(|(i, n)| names[..i].contains(n))
    {
        return Err("duplicate or invalid formula request fields".into());
    }
    Ok((names, values))
}
fn get<'a>(v: &'a Value, name: &str) -> Result<&'a Value, String> {
    let (n, v) = fields(v)?;
    n.iter()
        .position(|n| n == name)
        .map(|i| &v[i])
        .ok_or_else(|| format!("formula request missing {name}"))
}
fn text(v: &Value) -> Result<&str, String> {
    match &v.data {
        Data::Text(v) if v.len() == 1 => v[0]
            .as_deref()
            .ok_or_else(|| "formula text cannot be missing".into()),
        _ => Err("formula request requires scalar text".into()),
    }
}
fn list(v: &Value) -> Result<&[Value], String> {
    match &v.data {
        Data::List(v) => Ok(v),
        _ => Err("formula request requires a list".into()),
    }
}
fn integer(v: &Value) -> Result<usize, String> {
    match &v.data {
        Data::Number { values, .. }
            if values.len() == 1
                && values[0].is_finite()
                && values[0] >= 0.0
                && values[0] <= i32::MAX as f64
                && values[0].fract() == 0.0 =>
        {
            Ok(values[0] as usize)
        }
        _ => Err("formula request requires a nonnegative integer".into()),
    }
}
fn scalar(n: usize) -> Value {
    Value::new(Data::Number {
        values: vec![n as f64],
        integer: true,
    })
}
fn indices(v: Vec<usize>) -> Value {
    Value::new(Data::Number {
        values: v.into_iter().map(|i| i as f64 + 1.0).collect(),
        integer: true,
    })
}
fn named(names: &[&str], values: Vec<Value>) -> Value {
    let mut v = Value::new(Data::List(values));
    v.names = Some(names.iter().map(|n| (*n).into()).collect());
    v
}
fn ast(v: &Value, depth: usize) -> Result<Expression, String> {
    if depth > 128 {
        return Err("formula AST nesting exceeds 128".into());
    }
    match text(get(v, "type")?)? {
        "literal" => Ok(Expression::Literal(get(v, "value")?.clone())),
        "symbol" => Ok(Expression::Symbol(text(get(v, "name")?)?.into())),
        "call" => Ok(Expression::Call {
            function: text(get(v, "function_name")?)?.into(),
            arguments: list(get(v, "arguments")?)?
                .iter()
                .map(|a| {
                    let n = get(a, "name")?;
                    Ok(Argument {
                        name: if matches!(n.data, Data::Null) {
                            None
                        } else {
                            Some(text(n)?.into())
                        },
                        expression: ast(get(a, "expression")?, depth + 1)?,
                    })
                })
                .collect::<Result<_, String>>()?,
        }),
        _ => Err("invalid formula AST node type".into()),
    }
}
fn options(v: &Value) -> Result<EvaluationOptions, String> {
    Ok(EvaluationOptions {
        collation: text(get(v, "collation")?)?.into(),
    })
}
fn missing_policy(request: &Value) -> Result<frame::MissingPolicy, String> {
    match text(get(request, "na_action")?)? {
        "pass" => Ok(frame::MissingPolicy::Pass),
        "fail" => Ok(frame::MissingPolicy::Fail),
        "omit" => Ok(frame::MissingPolicy::Omit),
        _ => Err("unsupported formula missing-value policy".into()),
    }
}
fn frame_output(result: frame::Frame) -> Value {
    named(
        &["values", "original_rows", "kept", "omitted"],
        vec![
            Value::new(Data::List(result.values)),
            scalar(result.original_rows),
            indices(result.kept_rows),
            indices(result.omitted_rows),
        ],
    )
}
fn dispatch(request: &Value) -> Result<Value, String> {
    let operation = text(get(request, "operation")?)?;
    match operation {
        "evaluate" | "frame_evaluate" => {
            let expressions = list(get(request, "expressions")?)?
                .iter()
                .map(|v| ast(v, 0))
                .collect::<Result<Vec<_>, _>>()?;
            let (names, values) = fields(get(request, "variables")?)?;
            let bindings = names
                .iter()
                .cloned()
                .zip(values.iter().cloned())
                .collect::<Vec<_>>();
            if operation == "frame_evaluate" {
                return frame::evaluate(
                    &expressions,
                    &bindings,
                    &options(request)?,
                    integer(get(request, "data_rows")?)?,
                    missing_policy(request)?,
                )
                .map(frame_output)
                .map_err(|e| e.to_string());
            }
            expression::evaluate_many(&expressions, &bindings, &options(request)?)
                .map(|v| Value::new(Data::List(v)))
                .map_err(|e| e.to_string())
        }
        "frame" => {
            let values = list(get(request, "values")?)?.to_vec();
            let policy = missing_policy(request)?;
            let result = frame::assemble(values, integer(get(request, "data_rows")?)?, policy)
                .map_err(|e| e.to_string())?;
            Ok(frame_output(result))
        }
        "prepare" => {
            let response = integer(get(request, "response")?)?;
            frame::prepare_model_values(
                list(get(request, "values")?)?,
                response.checked_sub(1),
                &options(request)?,
            )
            .map(|v| Value::new(Data::List(v)))
            .map_err(|e| e.to_string())
        }
        "offsets" => {
            let values = list(get(request, "values")?)?;
            let Data::Number {
                values: indices, ..
            } = &get(request, "indices")?.data
            else {
                return Err("formula offset indices must be integer values".into());
            };
            let indices = indices
                .iter()
                .map(|v| {
                    if !v.is_finite() || *v < 1.0 || v.fract() != 0.0 || *v > values.len() as f64 {
                        return Err("invalid formula offset index".into());
                    }
                    Ok(*v as usize - 1)
                })
                .collect::<Result<Vec<_>, String>>()?;
            let base = get(request, "base")?;
            frame::offsets(
                values,
                &indices,
                (!matches!(base.data, Data::Null)).then_some(base),
            )
            .map(|v| v.unwrap_or_else(|| Value::new(Data::Null)))
            .map_err(|e| e.to_string())
        }
        _ => Err("unsupported formula request operation".into()),
    }
}
/// # Safety
/// C validates `input` for `size` bytes and a distinct writable owner slot. No
/// R address is retained. The returned owner must be consumed by formula_drop.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_formula_request(
    input: *const u8,
    size: usize,
    output: *mut *mut FormulaOutput,
) -> i32 {
    run_backend(|| {
        let slot = unsafe { output_slice(output, 1, "formula owner")? };
        if !slot[0].is_null() {
            return Err("formula output owner slot must be empty".into());
        }
        let bytes = unsafe { input_slice(input, size, "formula input")? };
        let request = wire::decode(bytes).map_err(|e| e.to_string())?;
        let result = dispatch(&request)?;
        let bytes = wire::encode(&result).map_err(|e| e.to_string())?;
        slot[0] = Box::into_raw(Box::new(FormulaOutput { bytes }));
        Ok(())
    })
}
/// # Safety
/// Owner is live; outputs are distinct writable C slots. The borrowed byte view
/// remains valid only until the unique owner is consumed, and is immutable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_formula_bytes(
    owner: *const FormulaOutput,
    bytes: *mut *const u8,
    size: *mut usize,
) -> i32 {
    run_backend(|| {
        let owner = unsafe { owner.as_ref() }.ok_or("formula owner is null")?;
        unsafe {
            output_slice(bytes, 1, "formula bytes")?[0] = owner.bytes.as_ptr();
            output_slice(size, 1, "formula size")?[0] = owner.bytes.len();
        }
        Ok(())
    })
}
/// # Safety
/// A non-null slot contains NULL or a unique live owner; clears it before drop.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_formula_drop(slot: *mut *mut FormulaOutput) {
    if slot.is_null() {
        return;
    }
    let owner = unsafe { ptr::replace(&mut *slot, ptr::null_mut()) };
    if !owner.is_null() {
        drop(unsafe { Box::from_raw(owner) });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn primitive_request_owner_is_published_only_on_success_and_consumed_once() {
        let literal = named(
            &["type", "value"],
            vec![
                Value::new(Data::Text(vec![Some("literal".into())])),
                Value::scalar(1.25),
            ],
        );
        let mut bindings = Value::new(Data::List(vec![]));
        bindings.names = Some(vec![]);
        let request = named(
            &["operation", "expressions", "variables", "collation"],
            vec![
                Value::new(Data::Text(vec![Some("evaluate".into())])),
                Value::new(Data::List(vec![literal])),
                bindings,
                Value::new(Data::Text(vec![Some("root".into())])),
            ],
        );
        let bytes = wire::encode(&request).unwrap();
        let mut owner = ptr::null_mut();
        let (mut data, mut size) = (ptr::null(), 0);
        unsafe {
            assert_ne!(rsdeseq2_r_formula_request(ptr::null(), 1, &mut owner), 0);
            assert!(owner.is_null());
            assert_ne!(
                rsdeseq2_r_formula_request(bytes.as_ptr(), bytes.len() - 1, &mut owner),
                0
            );
            assert!(owner.is_null());
            assert_eq!(
                rsdeseq2_r_formula_request(bytes.as_ptr(), bytes.len(), &mut owner),
                0
            );
            let saved = owner;
            assert_ne!(
                rsdeseq2_r_formula_request(bytes.as_ptr(), bytes.len(), &mut owner),
                0
            );
            assert_eq!(owner, saved);
            assert_eq!(rsdeseq2_r_formula_bytes(owner, &mut data, &mut size), 0);
            assert_eq!(
                wire::decode(std::slice::from_raw_parts(data, size)).unwrap(),
                Value::new(Data::List(vec![Value::scalar(1.25)]))
            );
            rsdeseq2_r_formula_drop(&mut owner);
            assert!(owner.is_null());
            rsdeseq2_r_formula_drop(&mut owner);
            rsdeseq2_r_formula_drop(ptr::null_mut());
            assert_ne!(rsdeseq2_r_formula_bytes(owner, &mut data, &mut size), 0);
        }
    }
}
