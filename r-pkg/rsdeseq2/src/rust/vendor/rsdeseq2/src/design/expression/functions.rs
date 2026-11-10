use super::{Argument, Data, Evaluator, Value, arguments, invalid};
use crate::errors::DeseqError;

pub(super) fn supported(name: &str) -> bool {
    matches!(
        name,
        "(" | "I"
            | "offset"
            | "+"
            | "-"
            | "*"
            | "/"
            | "^"
            | "%%"
            | "%/%"
            | "!"
            | "&"
            | "|"
            | "<"
            | "<="
            | ">"
            | ">="
            | "=="
            | "!="
            | "log"
            | "log2"
            | "log10"
            | "log1p"
            | "exp"
            | "expm1"
            | "sqrt"
            | "abs"
            | "sign"
            | "floor"
            | "ceiling"
            | "trunc"
            | "as.numeric"
            | "as.double"
            | "as.integer"
            | "as.logical"
            | "as.character"
            | "is.na"
            | "is.nan"
            | "is.finite"
            | "is.infinite"
            | "c"
            | "list"
            | "cbind"
            | "matrix"
            | "rep"
            | "seq_len"
            | "length"
            | "scale"
            | "poly"
            | "factor"
            | "ordered"
            | "as.factor"
            | "as.ordered"
            | "droplevels"
            | "relevel"
            | "mean"
            | "sum"
            | "sd"
            | "ifelse"
            | "pmax"
            | "pmin"
    )
}
impl Evaluator<'_> {
    pub(super) fn call(
        &mut self,
        name: &str,
        args: &[(Option<String>, Value)],
        expressions: &[Argument],
    ) -> Result<Value, DeseqError> {
        match name {
            "(" | "I" | "offset" => {
                let a = arguments(args, &["x"])?;
                let mut x = a[0]
                    .ok_or_else(|| invalid("formula function requires x"))?
                    .clone();
                if name == "I" {
                    let mut classes = vec![Some("AsIs".into())];
                    if let Some(v) = x.attribute("class") {
                        classes.extend(
                            v.text()?
                                .into_iter()
                                .filter(|s| s.as_deref() != Some("AsIs")),
                        )
                    }
                    x.set_attribute("class", Value::new(Data::Text(classes)))
                }
                Ok(x)
            }
            "+" | "-" | "!" if args.len() == 1 => super::arithmetic::unary(name, &args[0].1),
            "+" | "-" | "*" | "/" | "^" | "%%" | "%/%" | "&" | "|" | "<" | "<=" | ">" | ">="
            | "==" | "!=" => {
                if args.len() != 2 {
                    return Err(invalid("binary formula operation requires two arguments"));
                }
                if matches!(name, "<" | "<=" | ">" | ">=" | "==" | "!=") {
                    self.comparison(name, &args[0].1, &args[1].1)
                } else {
                    super::arithmetic::binary(name, &args[0].1, &args[1].1)
                }
            }
            "scale" | "poly" => self.transform(name, args),
            "factor" | "ordered" | "as.factor" | "as.ordered" | "droplevels" | "relevel" => {
                self.factor(name, args)
            }
            "c" | "list" | "cbind" | "matrix" | "rep" | "seq_len" | "length" => {
                self.shape(name, args, expressions)
            }
            "mean" | "sum" | "sd" => self.summary(name, args),
            "pmax" | "pmin" => super::extrema::parallel(name, args),
            _ => self.elementary(name, args),
        }
    }
    pub(super) fn elementary(
        &self,
        name: &str,
        args: &[(Option<String>, Value)],
    ) -> Result<Value, DeseqError> {
        let formals = if name == "log" {
            vec!["x", "base"]
        } else {
            vec!["x"]
        };
        let a = arguments(args, &formals)?;
        let input = a[0].ok_or_else(|| invalid("formula transform requires x"))?;
        if name == "as.character" {
            return self.characters(input);
        }
        if let Data::Text(text) = &input.data
            && matches!(
                name,
                "is.na" | "is.nan" | "is.finite" | "is.infinite" | "as.logical"
            )
        {
            let mut out = if name == "as.logical" {
                input.flattened()
            } else {
                input.clone()
            };
            out.data = Data::Logical(
                text.iter()
                    .map(|v| {
                        if name == "as.logical" {
                            match v.as_deref() {
                                Some("TRUE" | "True" | "true" | "T") => Some(true),
                                Some("FALSE" | "False" | "false" | "F") => Some(false),
                                _ => None,
                            }
                        } else {
                            Some(name == "is.na" && v.is_none())
                        }
                    })
                    .collect(),
            );
            out.attributes.clear();
            if name == "as.logical" {
                out.names = None;
            }
            return Ok(out);
        }
        if name == "as.numeric" || name == "as.double" || name == "as.integer" {
            let flattened = input.flattened();
            let input = &flattened;
            let values = if matches!(&input.data, Data::Text(_)) {
                input
                    .text()?
                    .iter()
                    .map(|x| {
                        x.as_ref()
                            .and_then(|s| s.trim().parse::<f64>().ok())
                            .unwrap_or(super::NA_REAL)
                    })
                    .collect()
            } else {
                input.numeric()?
            };
            let integer = name == "as.integer";
            let values = if integer {
                values
                    .iter()
                    .map(|v| {
                        super::super::numeric_transform::integer_covariate(*v)
                            .unwrap_or(super::NA_REAL)
                    })
                    .collect()
            } else {
                values
            };
            return Ok(Value::new(Data::Number { values, integer }));
        }
        let flattened;
        let input = if name == "as.logical" {
            flattened = input.flattened();
            &flattened
        } else {
            input
        };
        let values = input.numeric()?;
        if matches!(
            name,
            "is.na" | "is.nan" | "is.finite" | "is.infinite" | "as.logical"
        ) {
            let mut result = input.clone();
            result.data = Data::Logical(
                values
                    .iter()
                    .map(|v| {
                        if name == "as.logical" && v.is_nan() {
                            None
                        } else {
                            Some(match name {
                                "is.na" => v.is_nan(),
                                "is.nan" => v.is_nan() && !super::is_na(*v),
                                "is.finite" => v.is_finite(),
                                "is.infinite" => v.is_infinite(),
                                _ => *v != 0.0,
                            })
                        }
                    })
                    .collect(),
            );
            result.attributes.clear();
            if name == "as.logical" {
                result.names = None;
            }
            return Ok(result);
        }
        if input.has_class("factor") {
            return Err(invalid("numeric formula transform received a factor"));
        }
        let base = if name == "log" {
            a[1].map(|v| v.scalar_number())
                .transpose()?
                .unwrap_or(std::f64::consts::E)
        } else {
            std::f64::consts::E
        };
        let out = values
            .iter()
            .map(|&v| match name {
                "log" => {
                    if base == 2.0 {
                        v.log2()
                    } else if base == 10.0 {
                        v.log10()
                    } else if base == std::f64::consts::E {
                        v.ln()
                    } else {
                        v.ln() / base.ln()
                    }
                }
                "log2" => v.log2(),
                "log10" => v.log10(),
                "log1p" => crate::math::elementary::language_log1p(v),
                "exp" => crate::math::elementary::language_exp(v),
                "expm1" => crate::math::elementary::language_expm1(v),
                "sqrt" => v.sqrt(),
                "abs" => v.abs(),
                "sign" => {
                    if v.is_nan() {
                        v
                    } else if v == 0.0 {
                        0.0
                    } else {
                        v.signum()
                    }
                }
                "floor" => v.floor(),
                "ceiling" => v.ceil(),
                "trunc" => v.trunc(),
                _ => f64::NAN,
            })
            .collect();
        Ok(input.with_numbers(out, false))
    }
}
