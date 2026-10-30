use super::{Data, Evaluator, Value, arguments, invalid};
use crate::design::numeric_transform::{ScalePolicy, scale_expression};
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
            if input.has_class("factor") {
                return Err(invalid("scale requires numeric covariates"));
            }
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
            let mut result = input.with_numbers(scaled.values.into_values(), false);
            result.names = None;
            result.dimensions = Some([rows, cols]);
            result.dimnames = if input.dimensions.is_some() {
                input.dimnames.clone()
            } else {
                [input.names.clone(), None]
            };
            result.dimname_names = if input.dimensions.is_some() {
                input.dimname_names.clone()
            } else {
                [None, None]
            };
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
        self.polynomial(args)
    }
}
