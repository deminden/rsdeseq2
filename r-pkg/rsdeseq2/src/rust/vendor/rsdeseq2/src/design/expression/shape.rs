use super::{Argument, Data, Evaluator, Expression, Value, arguments, invalid};
use crate::errors::DeseqError;

impl Value {
    pub(super) fn flattened(&self) -> Self {
        let indices = if let Some([rows, cols]) = self.dimensions {
            (0..cols)
                .flat_map(|c| (0..rows).map(move |r| r * cols + c))
                .collect::<Vec<_>>()
        } else {
            (0..self.len()).collect()
        };
        self.selected(&indices)
    }
    fn selected(&self, indices: &[usize]) -> Self {
        let mut v = self.clone();
        v.dimensions = None;
        v.dimnames = [None, None];
        v.dimname_names = [None, None];
        v.names = self
            .names
            .as_ref()
            .map(|n| indices.iter().map(|&i| n[i].clone()).collect());
        v.data = match &self.data {
            Data::Number { values, integer } => Data::Number {
                values: indices.iter().map(|&i| values[i]).collect(),
                integer: *integer,
            },
            Data::Logical(x) => Data::Logical(indices.iter().map(|&i| x[i]).collect()),
            Data::Text(x) => Data::Text(indices.iter().map(|&i| x[i].clone()).collect()),
            Data::List(x) => Data::List(indices.iter().map(|&i| x[i].clone()).collect()),
            Data::Null => Data::Null,
        };
        v
    }
}
fn dimension(x: &Value) -> Result<usize, DeseqError> {
    let v = x.scalar_number()?;
    if !v.is_finite() || v < 0.0 || v > i32::MAX as f64 {
        return Err(invalid("formula dimension is outside its supported range"));
    }
    Ok(v as usize)
}
impl Evaluator<'_> {
    pub(super) fn shape(
        &self,
        name: &str,
        args: &[(Option<String>, Value)],
        expressions: &[Argument],
    ) -> Result<Value, DeseqError> {
        match name {
            "list" => {
                let mut v = Value::new(Data::List(args.iter().map(|(_, v)| v.clone()).collect()));
                if args.iter().any(|(n, _)| n.is_some()) {
                    v.names = Some(
                        args.iter()
                            .map(|(n, _)| n.clone().unwrap_or_default())
                            .collect(),
                    )
                }
                Ok(v)
            }
            "c" => {
                let values = args.iter().map(|(_, v)| v.flattened()).collect::<Vec<_>>();
                if values.iter().all(|v| matches!(v.data, Data::Null)) {
                    return Ok(Value::new(Data::Null));
                }
                super::checked_elements(
                    values
                        .iter()
                        .try_fold(0usize, |n, v| n.checked_add(v.len())),
                )?;
                let has_names = args
                    .iter()
                    .zip(&values)
                    .any(|((name, _), v)| !v.is_empty() && (name.is_some() || v.names.is_some()));
                let names = has_names.then(|| {
                    args.iter()
                        .zip(&values)
                        .flat_map(|((name, _), v)| {
                            let prefix = name.as_deref().filter(|s| !s.is_empty());
                            (0..v.len()).map(move |i| {
                                let inner = v
                                    .names
                                    .as_ref()
                                    .map(|n| n[i].as_str())
                                    .filter(|s| !s.is_empty());
                                match (prefix, inner) {
                                    (Some(p), Some(n)) => format!("{p}.{n}"),
                                    (Some(p), None) if v.len() == 1 => p.into(),
                                    (Some(p), None) => format!("{p}{}", i + 1),
                                    (None, Some(n)) => n.into(),
                                    _ => String::new(),
                                }
                            })
                        })
                        .collect()
                });
                let mut result = if values.iter().any(|v| matches!(&v.data, Data::Text(_))) {
                    let mut result = Vec::new();
                    for v in &values {
                        result.extend(self.characters(v)?.text()?)
                    }
                    Value::new(Data::Text(result))
                } else if values
                    .iter()
                    .all(|v| matches!(&v.data, Data::Logical(_) | Data::Null))
                {
                    let mut out = Vec::new();
                    for v in values {
                        if let Data::Logical(v) = v.data {
                            out.extend(v)
                        }
                    }
                    Value::new(Data::Logical(out))
                } else {
                    let integer = values.iter().all(|v| {
                        matches!(
                            &v.data,
                            Data::Number { integer: true, .. } | Data::Logical(_) | Data::Null
                        )
                    });
                    let mut out = Vec::new();
                    for v in values {
                        if !matches!(v.data, Data::Null) {
                            out.extend(v.numeric()?)
                        }
                    }
                    Value::new(Data::Number {
                        values: out,
                        integer,
                    })
                };
                result.names = names;
                Ok(result)
            }
            "length" => {
                let a = arguments(args, &["x"])?;
                Ok(Value::new(Data::Number {
                    values: vec![a[0].ok_or_else(|| invalid("length requires x"))?.len() as f64],
                    integer: true,
                }))
            }
            "seq_len" => {
                let a = arguments(args, &["length.out"])?;
                let n = dimension(a[0].ok_or_else(|| invalid("seq_len requires a length"))?)?;
                super::checked_elements(Some(n))?;
                Ok(Value::new(Data::Number {
                    values: (1..=n).map(|i| i as f64).collect(),
                    integer: true,
                }))
            }
            "rep" => {
                let a = arguments(args, &["x", "times", "length.out", "each"])?;
                let input = a[0].ok_or_else(|| invalid("rep requires x"))?.flattened();
                let each = a[3].map(dimension).transpose()?.unwrap_or(1);
                super::checked_elements(input.len().checked_mul(each))?;
                let mut indices = Vec::new();
                for i in 0..input.len() {
                    indices.extend(std::iter::repeat_n(i, each))
                }
                if let Some(length) = a[2] {
                    let n = dimension(length)?;
                    super::checked_elements(Some(n))?;
                    if indices.is_empty() && n > 0 {
                        return Err(invalid("rep cannot recycle empty values"));
                    }
                    indices = (0..n).map(|i| indices[i % indices.len()]).collect();
                } else if let Some(times) = a[1] {
                    let times = times.numeric()?;
                    let count =
                        |v: f64| -> Result<usize, DeseqError> { dimension(&Value::scalar(v)) };
                    if times.len() == 1 {
                        let n = count(times[0])?;
                        super::checked_elements(indices.len().checked_mul(n))?;
                        let original = indices.clone();
                        for _ in 1..n {
                            indices.extend_from_slice(&original)
                        }
                        if n == 0 {
                            indices.clear()
                        }
                    } else if times.len() == indices.len() {
                        let counts = times
                            .into_iter()
                            .map(count)
                            .collect::<Result<Vec<_>, _>>()?;
                        super::checked_elements(
                            counts.iter().try_fold(0usize, |sum, n| sum.checked_add(*n)),
                        )?;
                        indices = indices
                            .iter()
                            .zip(counts)
                            .flat_map(|(&i, n)| std::iter::repeat_n(i, n))
                            .collect();
                    } else {
                        return Err(invalid("rep times must be scalar or one per value"));
                    }
                }
                Ok(input.selected(&indices))
            }
            "matrix" => {
                let a = arguments(args, &["data", "nrow", "ncol", "byrow", "dimnames"])?;
                let default = Value::scalar(f64::NAN);
                let input = a[0].unwrap_or(&default).flattened();
                let row = a[1].map(dimension).transpose()?;
                let col = a[2].map(dimension).transpose()?;
                let rows = row.unwrap_or_else(|| {
                    if let Some(c) = col {
                        if c == 0 { 0 } else { input.len().div_ceil(c) }
                    } else {
                        1
                    }
                });
                let cols = col.unwrap_or_else(|| {
                    if rows == 0 {
                        0
                    } else {
                        input.len().div_ceil(rows)
                    }
                });
                let size = rows
                    .checked_mul(cols)
                    .ok_or_else(|| invalid("formula matrix shape overflow"))?;
                super::checked_elements(Some(size))?;
                let byrow = a[3].map(|v| v.bool_scalar()).transpose()?.unwrap_or(false);
                if input.is_empty() && size > 0 {
                    return Err(invalid("formula matrix requires values"));
                }
                let indices = (0..size)
                    .map(|i| {
                        let r = i / cols;
                        let c = i % cols;
                        (if byrow { i } else { r + c * rows }) % input.len()
                    })
                    .collect::<Vec<_>>();
                let mut out = input.selected(&indices);
                out.names = None;
                out.attributes.clear();
                out.dimensions = Some([rows, cols]);
                if let Some(v) = a[4]
                    && !matches!(v.data, Data::Null)
                {
                    if let Data::List(parts) = &v.data {
                        if parts.len() != 2 {
                            return Err(invalid("matrix dimnames must have two entries"));
                        }
                        for (i, p) in parts.iter().enumerate() {
                            if !matches!(p.data, Data::Null) {
                                out.dimnames[i] = Some(
                                    p.text()?
                                        .into_iter()
                                        .map(|x| x.unwrap_or_default())
                                        .collect(),
                                );
                                out.dimname_names[i] = p.names.clone();
                            }
                        }
                    } else {
                        return Err(invalid("matrix dimnames must be a list"));
                    }
                }
                Ok(out)
            }
            "cbind" => {
                let include = args
                    .iter()
                    .enumerate()
                    .filter(|(_, (n, _))| n.as_deref() != Some("deparse.level"))
                    .collect::<Vec<_>>();
                let rows = include
                    .iter()
                    .map(|(_, (_, v))| v.dimensions.map_or(v.len(), |d| d[0]))
                    .max()
                    .unwrap_or(0);
                let cols = include
                    .iter()
                    .map(|(_, (_, v))| v.dimensions.map_or(1, |d| d[1]))
                    .sum::<usize>();
                let mut values = vec![0.0; super::checked_elements(rows.checked_mul(cols))?];
                let mut names = Vec::new();
                let mut any_names = false;
                let mut offset = 0;
                let mut integer = true;
                let mut row_names = None;
                let mut row_name_names = None;
                for (index, (name, v)) in include {
                    let width = v.dimensions.map_or(1, |d| d[1]);
                    let n = v.dimensions.map_or(v.len(), |d| d[0]);
                    if row_names.is_none() && n == rows {
                        row_names = if v.dimensions.is_some() {
                            v.dimnames[0].clone()
                        } else {
                            v.names.clone()
                        };
                        row_name_names = if v.dimensions.is_some() {
                            v.dimname_names[0].clone()
                        } else {
                            None
                        };
                    }
                    if n == 0 && rows > 0 {
                        return Err(invalid("cbind cannot recycle empty values"));
                    }
                    if v.dimensions.is_some() && n != rows {
                        return Err(invalid("cbind matrix rows must agree"));
                    }
                    let data = v.numeric()?;
                    integer &= matches!(
                        &v.data,
                        Data::Number { integer: true, .. } | Data::Logical(_)
                    );
                    let automatic = if let Expression::Symbol(s) = &expressions[index].expression {
                        Some(s.clone())
                    } else {
                        None
                    };
                    for c in 0..width {
                        let label = v.dimnames[1].as_ref().map(|v| v[c].clone()).or_else(|| {
                            if width == 1 {
                                name.clone().or_else(|| automatic.clone())
                            } else {
                                None
                            }
                        });
                        any_names |= label.is_some();
                        names.push(label.unwrap_or_default());
                        for r in 0..rows {
                            values[r * cols + offset + c] = data[(r % n) * width + c]
                        }
                    }
                    offset += width;
                }
                let mut out = Value::new(Data::Number { values, integer });
                out.dimensions = Some([rows, cols]);
                out.dimnames[0] = row_names;
                out.dimname_names[0] = row_name_names;
                if any_names {
                    out.dimnames[1] = Some(names)
                }
                Ok(out)
            }
            _ => Err(invalid("unsupported formula shape operation")),
        }
    }
}
