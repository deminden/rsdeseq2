use super::{Argument, Data, Evaluator, NA_REAL, Value, arguments, invalid};
use crate::errors::DeseqError;

fn kind(value: &Value) -> Result<u8, DeseqError> {
    match &value.data {
        Data::Logical(_) | Data::Null => Ok(0),
        Data::Number { integer: true, .. } => Ok(1),
        Data::Number { integer: false, .. } => Ok(2),
        Data::Text(_) => Ok(3),
        Data::List(_) => Err(invalid("ifelse formula branches must be atomic")),
    }
}

impl Evaluator<'_> {
    pub(super) fn conditional(
        &mut self,
        inputs: &[Argument],
        depth: usize,
    ) -> Result<Value, DeseqError> {
        let inputs = inputs
            .iter()
            .map(|a| (a.name.clone(), &a.expression))
            .collect::<Vec<_>>();
        let bound = arguments(&inputs, &["test", "yes", "no"])?;
        let test = self.expression(
            bound[0].ok_or_else(|| invalid("ifelse requires test"))?,
            depth + 1,
        )?;
        let mut unshaped = test.clone();
        unshaped.dimensions = None;
        unshaped.dimnames = [None, None];
        unshaped.names = None;
        let logical = self.elementary("as.logical", &[(None, unshaped)])?;
        let Data::Logical(conditions) = logical.data else {
            unreachable!()
        };
        let mut output = test.clone();
        output.data = Data::Logical(conditions.clone());
        let any_yes = conditions.contains(&Some(true));
        let any_no = conditions.contains(&Some(false));
        let yes = if any_yes {
            Some(self.expression(
                bound[1].ok_or_else(|| invalid("ifelse requires yes"))?,
                depth + 1,
            )?)
        } else {
            None
        };
        let no = if any_no {
            Some(self.expression(
                bound[2].ok_or_else(|| invalid("ifelse requires no"))?,
                depth + 1,
            )?)
        } else {
            None
        };
        if [yes.as_ref(), no.as_ref()]
            .into_iter()
            .flatten()
            .any(|v| matches!(v.data, Data::Null))
        {
            return Err(invalid("ifelse selected branch has no replacement values"));
        }
        if !any_yes && !any_no {
            return Ok(output);
        }
        if test.len() == 1
            && test.dimensions.is_none()
            && test.names.is_none()
            && test.attributes.is_empty()
        {
            let selected = if any_yes { yes.as_ref() } else { no.as_ref() }.unwrap();
            if selected.len() == 1
                && selected.dimensions.is_none()
                && selected.names.is_none()
                && selected.attributes.is_empty()
            {
                return Ok(selected.clone());
            }
        }
        let output_kind = [yes.as_ref(), no.as_ref()]
            .into_iter()
            .flatten()
            .map(kind)
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .max()
            .unwrap_or(0);
        let prepare = |value: Option<Value>| -> Result<Option<Value>, DeseqError> {
            value
                .map(|v| {
                    if output_kind == 3 {
                        self.characters(&v)
                    } else {
                        Ok(v.flattened())
                    }
                })
                .transpose()
        };
        let yes = prepare(yes)?;
        let no = prepare(no)?;
        // R's branch recycling follows column order even when the primitive
        // result uses row-major storage.
        let positions = if let Some([rows, cols]) = output.dimensions {
            (0..output.len())
                .map(|i| i / cols + (i % cols) * rows)
                .collect::<Vec<_>>()
        } else {
            (0..output.len()).collect()
        };
        let select = |condition| {
            if condition == Some(true) {
                yes.as_ref()
            } else if condition == Some(false) {
                no.as_ref()
            } else {
                None
            }
        };
        output.data = if output_kind == 3 {
            Data::Text(
                conditions
                    .iter()
                    .zip(&positions)
                    .map(|(&condition, &i)| {
                        select(condition).and_then(|v| match &v.data {
                            Data::Text(values) if !values.is_empty() => {
                                values[i % values.len()].clone()
                            }
                            _ => None,
                        })
                    })
                    .collect(),
            )
        } else {
            let yes = yes.as_ref().map(Value::numeric).transpose()?;
            let no = no.as_ref().map(Value::numeric).transpose()?;
            let values = conditions
                .iter()
                .zip(&positions)
                .map(|(&condition, &i)| {
                    let source = match condition {
                        Some(true) => yes.as_ref(),
                        Some(false) => no.as_ref(),
                        None => None,
                    };
                    source
                        .filter(|v| !v.is_empty())
                        .map_or(NA_REAL, |v| v[i % v.len()])
                })
                .collect::<Vec<_>>();
            if output_kind == 0 {
                Data::Logical(
                    values
                        .iter()
                        .map(|v| if v.is_nan() { None } else { Some(*v != 0.0) })
                        .collect(),
                )
            } else {
                Data::Number {
                    values,
                    integer: output_kind == 1,
                }
            }
        };
        Ok(output)
    }
}
