use super::{Data, Evaluator, Value, arguments, invalid};
use crate::design::numeric_transform::{
    PolynomialCoefficients, ScalePolicy, polynomial_basis, scale_expression,
};
use crate::{
    errors::DeseqError,
    matrix::{MatrixOrder, MatrixView},
};

impl Evaluator<'_> {
    pub(super) fn transform(
        &self,
        name: &str,
        args: &[(Option<String>, Value)],
    ) -> Result<Value, DeseqError> {
        if name == "scale" {
            let a = arguments(args, &["x", "center", "scale"])?;
            let input = a[0].ok_or_else(|| invalid("scale requires x"))?;
            let values = input.numeric()?;
            let [rows, cols] = input.dimensions.unwrap_or([input.len(), 1]);
            let defaults = Value::new(Data::Logical(vec![Some(true)]));
            let center = a[1].unwrap_or(&defaults);
            let scale = a[2].unwrap_or(&defaults);
            let center_values = center.numeric()?;
            let scale_values = scale.numeric()?;
            let policy = |v: &Value| -> Result<Option<bool>, DeseqError> {
                if matches!(&v.data, Data::Logical(_)) {
                    Ok(Some(v.bool_scalar()?))
                } else {
                    Ok(None)
                }
            };
            let center = match policy(center)? {
                Some(true) => ScalePolicy::Automatic,
                Some(false) => ScalePolicy::Disabled,
                None => ScalePolicy::Supplied(&center_values),
            };
            let scale = match policy(scale)? {
                Some(true) => ScalePolicy::Automatic,
                Some(false) => ScalePolicy::Disabled,
                None => ScalePolicy::Supplied(&scale_values),
            };
            let scaled = scale_expression(
                MatrixView::new(rows, cols, &values, MatrixOrder::RowMajor)?,
                center,
                scale,
            )?;
            let mut result = Value::numbers(scaled.values.into_values());
            result.dimensions = Some([rows, cols]);
            result.dimnames = input.dimnames.clone();
            let attribute = |v: Vec<f64>, policy: ScalePolicy<'_>, original: &Value| {
                let mut x = Value::numbers(v);
                x.names = if matches!(policy, ScalePolicy::Automatic) {
                    input.dimnames[1].clone()
                } else {
                    original.names.clone()
                };
                x
            };
            if let Some(v) = scaled.center {
                result.set_attribute(
                    "scaled:center",
                    attribute(v, center, a[1].unwrap_or(&defaults)),
                )
            }
            if let Some(v) = scaled.scale {
                result.set_attribute(
                    "scaled:scale",
                    attribute(v, scale, a[2].unwrap_or(&defaults)),
                )
            }
            return Ok(result);
        }
        let a = arguments(args, &["x", "degree", "coefs", "raw", "simple"])?;
        let input = a[0].ok_or_else(|| invalid("poly requires x"))?;
        let degree = a[1].map(|v| v.scalar_number()).transpose()?.unwrap_or(1.0);
        if !degree.is_finite() || degree < 1.0 || degree > i32::MAX as f64 {
            return Err(invalid("poly degree must be a positive integer"));
        }
        let degree = degree as usize;
        let raw = a[3].map(|v| v.bool_scalar()).transpose()?.unwrap_or(false);
        let simple = a[4].map(|v| v.bool_scalar()).transpose()?.unwrap_or(false);
        let values = input.numeric()?;
        let coefficients = if let Some(v) = a[2] {
            if matches!(v.data, Data::Null) {
                None
            } else {
                let Data::List(fields) = &v.data else {
                    return Err(invalid("poly coefs must be a named list"));
                };
                let names = v
                    .names
                    .as_ref()
                    .ok_or_else(|| invalid("poly coefs require names"))?;
                let get = |name: &str| -> Result<Vec<f64>, DeseqError> {
                    let i = names
                        .iter()
                        .position(|n| n == name)
                        .ok_or_else(|| invalid(format!("poly coefs missing {name}")))?;
                    fields[i].numeric()
                };
                Some(PolynomialCoefficients {
                    alpha: get("alpha")?,
                    norm2: get("norm2")?,
                })
            }
        } else {
            None
        };
        let basis = polynomial_basis(&values, degree, raw, coefficients.as_ref())?;
        let mut result = Value::numbers(basis.values.into_values());
        result.dimensions = Some([values.len(), degree]);
        result.dimnames = [
            input.names.clone(),
            Some((1..=degree).map(|i| i.to_string()).collect()),
        ];
        if !simple {
            if let Some(c) = basis.coefficients {
                let mut coefs = Value::new(Data::List(vec![
                    Value::numbers(c.alpha),
                    Value::numbers(c.norm2),
                ]));
                coefs.names = Some(vec!["alpha".into(), "norm2".into()]);
                result.set_attribute("coefs", coefs)
            }
            result.set_attribute(
                "degree",
                Value::new(Data::Number {
                    values: (1..=degree).map(|i| i as f64).collect(),
                    integer: true,
                }),
            );
            result.set_attribute(
                "class",
                Value::new(Data::Text(vec![Some("poly".into()), Some("matrix".into())])),
            );
        }
        Ok(result)
    }
}
