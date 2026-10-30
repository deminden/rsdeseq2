use super::{Data, Evaluator, NA_REAL, Value, arguments, invalid};
use crate::{design::numeric_transform::scalar_mean, errors::DeseqError};
impl Evaluator<'_> {
    pub(super) fn summary(
        &self,
        name: &str,
        args: &[(Option<String>, Value)],
    ) -> Result<Value, DeseqError> {
        if name == "sum" {
            let controls = args
                .iter()
                .filter(|(n, _)| n.as_deref() == Some("na.rm"))
                .collect::<Vec<_>>();
            if controls.len() > 1 {
                return Err(invalid("duplicate sum na.rm argument"));
            }
            let remove = controls
                .first()
                .map(|(_, v)| v.bool_scalar())
                .transpose()?
                .unwrap_or(false);
            let mut integer = true;
            let mut total = 0.0;
            for (arg, value) in args {
                if arg.as_deref() == Some("na.rm") {
                    continue;
                }
                if value.has_class("factor") || value.has_class("Date") {
                    return Err(invalid("sum requires numeric or logical values"));
                }
                integer &= matches!(
                    value.data,
                    Data::Number { integer: true, .. } | Data::Logical(_) | Data::Null
                );
                for v in value.flattened().numeric()? {
                    if !remove || !v.is_nan() {
                        total += v;
                    }
                }
            }
            if integer && total.is_nan() {
                total = NA_REAL;
            } else if total <= i32::MIN as f64 || total > i32::MAX as f64 {
                integer = false;
            }
            return Ok(Value::new(Data::Number {
                values: vec![total],
                integer,
            }));
        }
        let formals = if name == "mean" {
            vec!["x", "trim", "na.rm"]
        } else {
            vec!["x", "na.rm"]
        };
        let a = arguments(args, &formals)?;
        let input = a[0].ok_or_else(|| invalid("formula summary requires x"))?;
        if name == "mean" && matches!(input.data, Data::Null | Data::Text(_) | Data::List(_)) {
            return Ok(Value::scalar(NA_REAL));
        }
        if input.has_class("factor") {
            if name == "mean" {
                return Ok(Value::scalar(NA_REAL));
            }
            return Err(invalid("sd requires numeric covariates"));
        }
        let mut x = input.flattened().numeric()?;
        let remove = a[if name == "mean" { 2 } else { 1 }]
            .map(Value::bool_scalar)
            .transpose()?
            .unwrap_or(false);
        if remove {
            x.retain(|v| !v.is_nan());
        }
        let value = if name == "sd" {
            if x.len() < 2 || x.iter().any(|v| v.is_nan()) {
                NA_REAL
            } else {
                let mean = scalar_mean(&x);
                (x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (x.len() - 1) as f64)
                    .sqrt()
            }
        } else {
            let trim = if let Some(v) = a[1] {
                if !matches!(v.data, Data::Number { .. }) {
                    return Err(invalid("mean trim must be numeric"));
                }
                v.scalar_number()?
            } else {
                0.0
            };
            if trim.is_nan() {
                return Err(invalid("mean trim must not be missing"));
            }
            if trim > 0.0 && !x.is_empty() {
                if x.iter().any(|v| v.is_nan()) {
                    return Ok(Value::scalar(NA_REAL));
                }
                x.sort_by(f64::total_cmp);
                if trim >= 0.5 {
                    let middle = x.len() / 2;
                    if x.len() % 2 == 0 {
                        scalar_mean(&x[middle - 1..=middle])
                    } else {
                        x[middle]
                    }
                } else {
                    let remove = (x.len() as f64 * trim).floor() as usize;
                    scalar_mean(&x[remove..x.len() - remove])
                }
            } else {
                scalar_mean(&x)
            }
        };
        Ok(Value::scalar(value))
    }
}
