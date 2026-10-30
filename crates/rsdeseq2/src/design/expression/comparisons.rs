use super::{Data, Evaluator, Value, invalid};
use crate::errors::DeseqError;
use std::cmp::Ordering;
impl Evaluator<'_> {
    pub(super) fn comparison(
        &self,
        name: &str,
        left: &Value,
        right: &Value,
    ) -> Result<Value, DeseqError> {
        let factor = left.has_class("factor") || right.has_class("factor");
        if !factor && !matches!(left.data, Data::Text(_)) && !matches!(right.data, Data::Text(_)) {
            let mut result = super::arithmetic::binary(name, left, right)?;
            if left.has_class("Date") || right.has_class("Date") {
                result.attributes.clear();
            }
            return Ok(result);
        }
        if left.dimensions.is_some()
            && right.dimensions.is_some()
            && left.dimensions != right.dimensions
        {
            return Err(invalid("non-conformable formula arrays"));
        }
        let characters = |v: &Value| {
            let mut x = v.clone();
            x.dimensions = None;
            x.dimnames = [None, None];
            x.dimname_names = [None, None];
            self.characters(&x)?.text()
        };
        let a = characters(left)?;
        let b = characters(right)?;
        let n = if a.is_empty() || b.is_empty() {
            0
        } else {
            a.len().max(b.len())
        };
        let template = if left.dimensions.is_some() {
            left
        } else if right.dimensions.is_some() || b.len() > a.len() {
            right
        } else {
            left
        };
        if template.dimensions.is_some() && template.len() != n {
            return Err(invalid("formula vector is longer than its matrix"));
        }
        let equal = matches!(name, "==" | "!=");
        let levels = |v: &Value| {
            v.attribute("levels")
                .ok_or_else(|| invalid("factor levels missing"))?
                .text()
        };
        let ordered = left.has_class("ordered") || right.has_class("ordered");
        if factor && !equal && !ordered {
            return Err(invalid("relational comparison requires ordered factors"));
        }
        if left.has_class("factor") && right.has_class("factor") {
            let x = levels(left)?;
            let y = levels(right)?;
            if x.len() != y.len() || x.iter().any(|v| !y.contains(v)) || ordered && x != y {
                return Err(invalid("comparison factors have different level sets"));
            }
        }
        let ordered_levels = if !equal && ordered {
            Some(levels(if left.has_class("ordered") {
                left
            } else {
                right
            })?)
        } else {
            None
        };
        let values = (0..n)
            .map(|i| {
                let (Some(a), Some(b)) = (&a[i % a.len()], &b[i % b.len()]) else {
                    return None;
                };
                let order = if let Some(levels) = &ordered_levels {
                    let x = levels.iter().position(|v| v.as_deref() == Some(a));
                    let y = levels.iter().position(|v| v.as_deref() == Some(b));
                    let (Some(x), Some(y)) = (x, y) else {
                        return None;
                    };
                    x.cmp(&y)
                } else if equal {
                    a.cmp(b)
                } else {
                    self.collator
                        .as_ref()
                        .map_or_else(|| a.cmp(b), |c| c.compare(a, b))
                };
                Some(match name {
                    "==" => order == Ordering::Equal,
                    "!=" => order != Ordering::Equal,
                    "<" => order == Ordering::Less,
                    "<=" => order != Ordering::Greater,
                    ">" => order == Ordering::Greater,
                    ">=" => order != Ordering::Less,
                    _ => unreachable!(),
                })
            })
            .collect();
        let mut result = template.clone();
        result.data = Data::Logical(values);
        result.attributes.clear();
        if a.len() == b.len() {
            result.names = left.names.clone().or_else(|| right.names.clone());
        }
        if n == 0 {
            result.dimensions = None;
            result.names = None;
            result.dimnames = [None, None];
            result.dimname_names = [None, None];
        }
        Ok(result)
    }
}
