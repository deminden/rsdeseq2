use super::{Data, Evaluator, Value, arguments, invalid};
use crate::errors::DeseqError;

fn numeric_label(x: f64, integer: bool) -> Option<String> {
    if x.is_nan() {
        return if super::is_na(x) {
            None
        } else {
            Some("NaN".into())
        };
    }
    if x == 0.0 {
        return Some("0".into());
    }
    if x.is_infinite() {
        return Some(if x > 0.0 { "Inf" } else { "-Inf" }.into());
    }
    if integer {
        return Some(format!("{x:.0}"));
    }
    let scientific = format!("{x:.14e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent = exponent.parse::<i32>().unwrap();
    let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
    let sci = format!(
        "{mantissa}e{}{:02}",
        if exponent >= 0 { "+" } else { "-" },
        exponent.abs()
    );
    let decimals = (14 - exponent).max(0) as usize;
    let mut fixed = format!("{x:.decimals$}");
    if fixed.contains('.') {
        fixed = fixed
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
    Some(if fixed.len() <= sci.len() { fixed } else { sci })
}
impl Evaluator<'_> {
    pub(super) fn characters(&self, input: &Value) -> Result<Value, DeseqError> {
        let input = input.flattened();
        let values = if input.has_class("factor") {
            let levels = input
                .attribute("levels")
                .ok_or_else(|| invalid("factor levels missing"))?
                .text()?;
            input
                .numeric()?
                .iter()
                .map(|v| {
                    if v.is_nan() {
                        Ok(None)
                    } else {
                        let i = *v as usize;
                        if *v != i as f64 || i == 0 || i > levels.len() {
                            Err(invalid("factor code outside levels"))
                        } else {
                            Ok(levels[i - 1].clone())
                        }
                    }
                })
                .collect::<Result<Vec<_>, DeseqError>>()?
        } else {
            match &input.data {
                Data::Text(v) => v.clone(),
                Data::Logical(v) => v
                    .iter()
                    .map(|v| v.map(|b| if b { "TRUE" } else { "FALSE" }.into()))
                    .collect(),
                Data::Number { values, integer } => {
                    values.iter().map(|&v| numeric_label(v, *integer)).collect()
                }
                Data::Null => Vec::new(),
                _ => return Err(invalid("cannot coerce formula list to text")),
            }
        };
        Ok(Value::new(Data::Text(values)))
    }
    pub(super) fn factor(
        &self,
        name: &str,
        args: &[(Option<String>, Value)],
    ) -> Result<Value, DeseqError> {
        let formals = if name == "relevel" {
            vec!["x", "ref"]
        } else {
            vec!["x", "levels", "labels", "exclude", "ordered", "nmax"]
        };
        let a = if matches!(name, "as.factor" | "as.ordered") {
            let mut a = arguments(args, &["x"])?;
            a.resize(6, None);
            a
        } else if name == "droplevels" {
            let a = arguments(args, &["x", "exclude"])?;
            vec![a[0], None, None, a[1], None, None]
        } else {
            arguments(args, &formals)?
        };
        let input = a[0].ok_or_else(|| invalid("factor transform requires x"))?;
        if name == "as.factor" && input.has_class("factor")
            || name == "as.ordered" && input.has_class("ordered")
        {
            return Ok(input.clone());
        }
        let values = self.characters(input)?.text()?;
        let original = input.attribute("levels").map(|v| v.text()).transpose()?;
        let ordered = if matches!(name, "ordered" | "as.ordered") {
            true
        } else if name == "relevel" {
            if input.has_class("ordered") {
                return Err(invalid("relevel requires an unordered factor"));
            }
            false
        } else {
            a.get(4)
                .copied()
                .flatten()
                .map(|v| v.bool_scalar())
                .transpose()?
                .unwrap_or(input.has_class("ordered"))
        };
        let mut levels = if name == "relevel" {
            let mut levels = original.ok_or_else(|| invalid("relevel requires a factor"))?;
            let reference = a[1].ok_or_else(|| invalid("relevel requires ref"))?;
            let index = if matches!(reference.data, Data::Text(_)) {
                let r = reference.text()?;
                if r.len() != 1 {
                    return Err(invalid("relevel ref must have length one"));
                }
                levels
                    .iter()
                    .position(|v| *v == r[0])
                    .ok_or_else(|| invalid("relevel ref is not a level"))?
            } else {
                let x = reference.scalar_number()?;
                if !x.is_finite() || x < 1.0 || x > levels.len() as f64 {
                    return Err(invalid("relevel ref is outside levels"));
                }
                x as usize - 1
            };
            let reference = levels.remove(index);
            levels.insert(0, reference);
            levels
        } else if let Some(explicit) =
            a[1].filter(|_| !matches!(name, "as.factor" | "as.ordered" | "droplevels"))
        {
            self.characters(explicit)?.text()?
        } else if let Some(levels) = original {
            levels.into_iter().filter(|v| values.contains(v)).collect()
        } else if matches!(&input.data, Data::Number { .. }) {
            let mut numbers = input.numeric()?;
            numbers.retain(|v| !super::is_na(*v));
            numbers.sort_by(f64::total_cmp);
            numbers.dedup_by(|a, b| *a == *b || a.is_nan() && b.is_nan());
            let integer = matches!(&input.data, Data::Number { integer: true, .. });
            numbers
                .into_iter()
                .map(|v| numeric_label(v, integer))
                .collect()
        } else {
            let mut levels = values.iter().filter_map(Clone::clone).collect::<Vec<_>>();
            levels.sort_by(|a, b| {
                self.collator
                    .as_ref()
                    .map_or_else(|| a.cmp(b), |c| c.compare(a, b))
            });
            levels.dedup();
            let mut levels = levels.into_iter().map(Some).collect::<Vec<_>>();
            if values.contains(&None) {
                levels.push(None);
            }
            levels
        };
        let exclude = if name == "relevel" {
            Vec::new()
        } else {
            match a[3] {
                Some(v) if matches!(v.data, Data::Null) => Vec::new(),
                Some(v) => self.characters(v)?.text()?,
                None => vec![None],
            }
        };
        levels.retain(|v| !exclude.contains(v));
        for (i, v) in levels.iter().enumerate() {
            if levels[..i].contains(v) {
                return Err(invalid("factor levels must be unique"));
            }
        }
        let labels = if name != "relevel" {
            a[2].map(|v| self.characters(v).and_then(|v| v.text()))
                .transpose()?
        } else {
            None
        };
        let labels = if let Some(labels) = labels {
            if labels.len() == levels.len() {
                labels
            } else if labels.len() == 1 {
                let prefix = labels[0].as_deref().unwrap_or("NA");
                (1..=levels.len())
                    .map(|i| Some(format!("{prefix}{i}")))
                    .collect()
            } else {
                return Err(invalid(
                    "factor labels must match the level count or be a prefix",
                ));
            }
        } else {
            levels.clone()
        };
        let mut output_levels = Vec::new();
        for v in &labels {
            if !output_levels.contains(v) {
                output_levels.push(v.clone())
            }
        }
        let codes = values
            .iter()
            .map(|v| {
                levels
                    .iter()
                    .position(|l| l == v)
                    .and_then(|i| output_levels.iter().position(|v| *v == labels[i]))
                    .map_or(super::NA_REAL, |i| i as f64 + 1.0)
            })
            .collect();
        let mut result = Value::new(Data::Number {
            values: codes,
            integer: true,
        });
        result.names = input.names.clone();
        result.set_attribute("levels", Value::new(Data::Text(output_levels)));
        result.set_attribute(
            "class",
            Value::new(Data::Text(if ordered {
                vec![Some("ordered".into()), Some("factor".into())]
            } else {
                vec![Some("factor".into())]
            })),
        );
        Ok(result)
    }
}
