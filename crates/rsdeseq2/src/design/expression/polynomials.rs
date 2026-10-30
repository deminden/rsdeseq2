use super::{Data, Evaluator, Value, invalid};
use crate::design::numeric_transform::{PolynomialCoefficients, polynomial_expression};
use crate::errors::DeseqError;
impl Evaluator<'_> {
    pub(super) fn polynomial(&self, args: &[(Option<String>, Value)]) -> Result<Value, DeseqError> {
        // Only x precedes `...`; controls after it require exact names. A lone
        // scalar in dots is the historical positional degree convention.
        let formals = ["x", "degree", "coefs", "raw", "simple"];
        let mut slots = [None; 5];
        let mut dots = Vec::new();
        for (index, (name, value)) in args.iter().enumerate() {
            if let Some(i) = name
                .as_deref()
                .and_then(|n| formals.iter().position(|f| *f == n))
            {
                if slots[i].replace(value).is_some() {
                    return Err(invalid("duplicate polynomial argument"));
                }
            } else {
                dots.push(index);
            }
        }
        if slots[0].is_none()
            && let Some(position) = dots.iter().position(|i| args[*i].0.is_none())
        {
            slots[0] = Some(&args[dots.remove(position)].1);
        }
        let input = slots[0].ok_or_else(|| invalid("poly requires x"))?;
        let scalar_degree = dots.len() == 1 && args[dots[0]].1.len() == 1;
        let degree = if scalar_degree {
            args[dots[0]].1.scalar_number()?
        } else {
            slots[1]
                .map(Value::scalar_number)
                .transpose()?
                .unwrap_or(1.0)
        };
        if !degree.is_finite() || degree < 1.0 || degree > i32::MAX as f64 {
            return Err(invalid("poly degree must be positive and finite"));
        }
        let degree = degree as usize;
        let raw = slots[3]
            .map(Value::bool_scalar)
            .transpose()?
            .unwrap_or(false);
        let simple = slots[4]
            .map(Value::bool_scalar)
            .transpose()?
            .unwrap_or(false);
        let coefficients = slots[2].filter(|v| !matches!(v.data, Data::Null));
        if input.dimensions.is_none() && (dots.is_empty() || scalar_degree) {
            return self.univariate(input, degree, coefficients, raw, simple);
        }
        let mut columns = Vec::new();
        let mut append = |v: &Value| -> Result<(), DeseqError> {
            if let Some([rows, cols]) = v.dimensions {
                let values = v.numeric()?;
                for c in 0..cols {
                    columns.push(Value::numbers(
                        (0..rows).map(|r| values[r * cols + c]).collect(),
                    ));
                }
            } else {
                columns.push(v.clone());
            }
            Ok(())
        };
        append(input)?;
        if !scalar_degree {
            for i in dots {
                append(&args[i].1)?;
            }
        }
        self.multivariate(&columns, degree, coefficients, raw)
    }

    fn multivariate(
        &self,
        columns: &[Value],
        degree: usize,
        coefficients: Option<&Value>,
        raw: bool,
    ) -> Result<Value, DeseqError> {
        let nd = columns.len();
        if nd == 0 {
            return Err(invalid("poly requires at least one column"));
        }
        let rows = columns[0].len();
        if columns.iter().any(|v| v.len() != rows) {
            return Err(invalid("polynomial arguments must have the same length"));
        }
        let coefficients = coefficients
            .map(|v| match &v.data {
                Data::List(v) if v.len() == nd => Ok(v.as_slice()),
                _ => Err(invalid("wrong number of polynomial coefficient columns")),
            })
            .transpose()?;
        if nd == 1 && coefficients.is_none() {
            return self.univariate(&columns[0], degree, None, raw, false);
        }
        let mut count = 1u128;
        for i in 1..=nd {
            count = count
                .checked_mul((degree + i) as u128)
                .ok_or_else(|| invalid("polynomial size overflow"))?
                / (i as u128);
        }
        let count = usize::try_from(count - 1).map_err(|_| invalid("polynomial size overflow"))?;
        super::checked_elements(rows.checked_mul(count))?;
        super::checked_elements(nd.checked_mul(count))?;
        let mut exponents = Vec::with_capacity(count);
        let mut powers = vec![0; nd];
        fn grid(axis: usize, remaining: usize, powers: &mut [usize], out: &mut Vec<Vec<usize>>) {
            for power in 0..=remaining {
                powers[axis] = power;
                if axis == 0 {
                    if powers.iter().any(|p| *p != 0) {
                        out.push(powers.to_vec());
                    }
                } else {
                    grid(axis - 1, remaining - power, powers, out);
                }
            }
        }
        grid(nd - 1, degree, &mut powers, &mut exponents);
        let mut basis = Vec::with_capacity(nd);
        let mut output_coefs = Vec::new();
        for (i, column) in columns.iter().enumerate() {
            let evaluated = self.univariate(
                column,
                degree,
                coefficients.map(|c| &c[i]),
                if coefficients.is_some() { false } else { raw },
                coefficients.is_some() || raw,
            )?;
            if coefficients.is_none() && !raw {
                output_coefs.push(
                    evaluated
                        .attribute("coefs")
                        .ok_or_else(|| invalid("polynomial recurrence metadata missing"))?
                        .clone(),
                );
            }
            basis.push(evaluated.numeric()?);
        }
        let mut values = vec![0.0; rows * count];
        for (c, powers) in exponents.iter().enumerate() {
            for r in 0..rows {
                let term = |i: usize| {
                    if powers[i] == 0 {
                        1.0
                    } else {
                        basis[i][r * degree + powers[i] - 1]
                    }
                };
                let mut product = term(0);
                for i in 1..nd {
                    product *= term(i);
                }
                values[r * count + c] = product;
            }
        }
        let mut result = Value::numbers(values);
        result.dimensions = Some([rows, count]);
        result.dimnames = [
            if coefficients.is_none() {
                columns[0].names.clone()
            } else {
                None
            },
            Some(
                exponents
                    .iter()
                    .map(|p| p.iter().map(usize::to_string).collect::<Vec<_>>().join("."))
                    .collect(),
            ),
        ];
        let mut weights = Vec::with_capacity(nd);
        let mut weight = 1usize;
        for _ in 0..nd {
            weights.push(weight);
            weight = weight
                .checked_mul(degree + 1)
                .ok_or_else(|| invalid("polynomial grid metadata overflow"))?;
        }
        result.dimname_names[1] = Some(
            exponents
                .iter()
                .map(|p| {
                    let index = p
                        .iter()
                        .zip(&weights)
                        .try_fold(1usize, |sum, (p, w)| {
                            p.checked_mul(*w).and_then(|v| sum.checked_add(v))
                        })
                        .ok_or_else(|| invalid("polynomial grid metadata overflow"))?;
                    Ok(index.to_string())
                })
                .collect::<Result<_, DeseqError>>()?,
        );
        if coefficients.is_none() {
            result.set_attribute(
                "degree",
                Value::numbers(
                    exponents
                        .iter()
                        .map(|p| p.iter().sum::<usize>() as f64)
                        .collect(),
                ),
            );
            if !raw {
                result.set_attribute("coefs", Value::new(Data::List(output_coefs)));
            }
            result.set_attribute(
                "class",
                Value::new(Data::Text(vec![Some("poly".into()), Some("matrix".into())])),
            );
        }
        Ok(result)
    }
    fn univariate(
        &self,
        input: &Value,
        degree: usize,
        coefficients: Option<&Value>,
        raw: bool,
        simple: bool,
    ) -> Result<Value, DeseqError> {
        if input.has_class("factor") && !input.has_class("ordered") {
            return Err(invalid("poly requires numeric or ordered covariates"));
        }
        let values = input.numeric()?;
        super::checked_elements(values.len().checked_mul(degree.saturating_add(1)))?;
        let coefficients = if let Some(v) = coefficients.filter(|_| !raw) {
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
        let basis = polynomial_expression(&values, degree, raw, coefficients.as_ref())?;
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
