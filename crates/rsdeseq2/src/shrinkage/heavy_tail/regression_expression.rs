//! Restricted regression functions are syntax/data, never foreign callbacks.
use super::super::expression_likelihood::{scalar, validate};
use super::{density, prior::invalid};
use crate::{
    DeseqError,
    design::expression::{self, Argument, Data, EvaluationOptions, Expression, Value},
    matrix::MatrixView,
};
#[derive(Clone, Debug)]
pub struct RestrictedRegressionLikelihood {
    arguments: Vec<String>,
    expression: Expression,
    constants: Vec<(String, Value)>,
}
fn special(function: &str) -> Option<&str> {
    let (owner, name) = function
        .split_once("::")
        .map_or((None, function), |(a, b)| (Some(a), b));
    let source = match name {
        "dnorm" | "dnbinom" | "dbinom" | "dt" => "stats",
        "%*%" | "lgamma" | "as.numeric" | "as.double" => "base",
        _ => return None,
    };
    (owner.is_none() || owner == Some(source)).then_some(name)
}
fn formals(name: &str) -> &'static [&'static str] {
    match name {
        "dnorm" => &["x", "mean", "sd", "log"],
        "dnbinom" => &["x", "size", "prob", "mu", "log"],
        "dbinom" => &["x", "size", "prob", "log"],
        "dt" => &["x", "df", "ncp", "log"],
        _ => &["x"],
    }
}
fn validate_node(
    e: &Expression,
    bindings: &[&str],
    depth: usize,
    nodes: &mut usize,
) -> Result<(), DeseqError> {
    *nodes += 1;
    if depth > 128 || *nodes > 4096 {
        return Err(invalid("regression expression exceeds nesting/node limit"));
    }
    if let Expression::Call {
        function,
        arguments,
    } = e
    {
        if let Some(name) = special(function) {
            if name == "%*%" {
                if arguments.len() != 2 || arguments.iter().any(|a| a.name.is_some()) {
                    return Err(invalid("matrix product requires two positional operands"));
                }
            } else {
                let args = arguments
                    .iter()
                    .map(|a| (a.name.clone(), ()))
                    .collect::<Vec<_>>();
                let bound = expression::arguments(&args, formals(name))?;
                if bound[0].is_none()
                    || (matches!(name, "dnbinom" | "dbinom" | "dt") && bound[1].is_none())
                    || (name == "dbinom" && bound[2].is_none())
                {
                    return Err(invalid("missing likelihood density argument"));
                }
                if name == "dnbinom" && (bound[2].is_some() == bound[3].is_some()) {
                    return Err(invalid(
                        "negative-binomial density requires exactly one of prob or mu",
                    ));
                }
            }
        } else {
            let shell = Expression::Call {
                function: function.clone(),
                arguments: arguments
                    .iter()
                    .map(|a| Argument {
                        name: a.name.clone(),
                        expression: Expression::Literal(Value::numbers(vec![1.0])),
                    })
                    .collect(),
            };
            validate(&shell, bindings, 0, &mut 0)?;
        }
        for a in arguments {
            validate_node(&a.expression, bindings, depth + 1, nodes)?;
        }
        Ok(())
    } else {
        validate(e, bindings, 0, &mut 0)
    }
}
fn numeric(v: &Value) -> Result<Vec<f64>, DeseqError> {
    match &v.data {
        Data::Number { values, .. } => Ok(values.clone()),
        Data::Logical(v) => Ok(v
            .iter()
            .map(|v| {
                v.map(f64::from)
                    .unwrap_or(super::super::location_missing::NA)
            })
            .collect()),
        Data::Null => Ok(vec![]),
        _ => Err(invalid("numeric likelihood values are required")),
    }
}
fn product(a: &Value, b: &Value) -> Result<Value, DeseqError> {
    let av = numeric(a)?;
    let bv = numeric(b)?;
    let (ar, ac, br, bc) = match (a.dimensions, b.dimensions) {
        (Some([ar, ac]), Some([br, bc])) => (ar, ac, br, bc),
        (Some([ar, ac]), None) if ac == bv.len() => (ar, ac, bv.len(), 1),
        (Some([ar, 1]), None) => (ar, 1, 1, bv.len()),
        (None, Some([br, bc])) if av.len() == br => (1, av.len(), br, bc),
        (None, Some([1, bc])) => (av.len(), 1, 1, bc),
        (None, None) if av.len() == bv.len() => (1, av.len(), bv.len(), 1),
        _ => return Err(invalid("nonconformable likelihood matrix product")),
    };
    if ac != br {
        return Err(invalid("nonconformable likelihood matrix product"));
    }
    let len = ar
        .checked_mul(bc)
        .ok_or_else(|| invalid("likelihood matrix product shape overflows"))?;
    let values = (0..len)
        .map(|k| {
            let i = k / bc;
            let j = k % bc;
            (0..ac).map(|q| av[i * ac + q] * bv[q * bc + j]).sum()
        })
        .collect();
    let mut out = Value::numbers(values);
    out.dimensions = Some([ar, bc]);
    Ok(out)
}
fn distribution(name: &str, args: &[(Option<String>, Value)]) -> Result<Value, DeseqError> {
    let f = formals(name);
    let bound = expression::arguments(args, f)?;
    let required =
        |index: usize| bound[index].ok_or_else(|| invalid("missing likelihood density argument"));
    required(0)?;
    if matches!(name, "dnbinom" | "dbinom" | "dt") {
        required(1)?;
    }
    if name == "dbinom" {
        required(2)?;
    }
    if name == "dnbinom" && (bound[2].is_some() == bound[3].is_some()) {
        return Err(invalid(
            "negative-binomial density requires exactly one of prob or mu",
        ));
    }
    let log = bound[f.len() - 1].map(scalar).transpose()?.unwrap_or(0.0) != 0.0;
    let defaults = match name {
        "dnorm" => vec![f64::NAN, 0.0, 1.0, 0.0],
        "dt" => vec![f64::NAN, f64::NAN, 0.0, 0.0],
        _ => vec![f64::NAN; f.len()],
    };
    let values = bound
        .iter()
        .enumerate()
        .map(|(i, v)| v.map(numeric).unwrap_or_else(|| Ok(vec![defaults[i]])))
        .collect::<Result<Vec<_>, _>>()?;
    let used = values
        .iter()
        .take(f.len() - 1)
        .enumerate()
        .filter(|(i, _)| name != "dnbinom" || *i != if bound[3].is_some() { 2 } else { 3 })
        .map(|(_, v)| v)
        .collect::<Vec<_>>();
    if used.iter().any(|v| v.is_empty()) {
        return Ok(Value::numbers(vec![]));
    }
    let n = used.iter().map(|v| v.len()).max().unwrap();
    let mut out = Value::numbers(
        (0..n)
            .map(|i| {
                let v = values.iter().map(|v| v[i % v.len()]).collect::<Vec<_>>();
                if used
                    .iter()
                    .any(|a| super::super::location_missing::is_na(a[i % a.len()]))
                {
                    return super::super::location_missing::NA;
                }
                let value = match name {
                    "dnorm" => density::normal(v[0], v[1], v[2]),
                    "dbinom" => density::binomial(v[0], v[1], v[2]),
                    "dt" => density::student(v[0], v[1], v[2]),
                    "dnbinom" => density::negative_binomial(
                        v[0],
                        v[1],
                        if bound[3].is_some() {
                            v[3]
                        } else {
                            v[1] * (1.0 - v[2]) / v[2]
                        },
                    ),
                    _ => unreachable!(),
                };
                if log { value } else { value.exp() }
            })
            .collect(),
    );
    if let Some(template) = bound.iter().flatten().find(|v| v.len() == n) {
        out.dimensions = template.dimensions;
    }
    Ok(out)
}
impl RestrictedRegressionLikelihood {
    /// Import a five-argument (y,x,beta,param,offset) function without executing it.
    pub fn new(
        arguments: Vec<String>,
        expression: Expression,
        constants: Vec<(String, Value)>,
    ) -> Result<Self, DeseqError> {
        if arguments.len() != 5 {
            return Err(invalid("regression likelihood requires five arguments"));
        }
        let mut names = vec![];
        for name in arguments.iter().chain(constants.iter().map(|(n, _)| n)) {
            if name.is_empty() || names.contains(&name.as_str()) {
                return Err(invalid(
                    "regression function bindings must be unique and nonempty",
                ));
            }
            names.push(name.as_str());
        }
        for (_, v) in &constants {
            scalar(v)?;
        }
        validate_node(&expression, &names, 0, &mut 0)?;
        Ok(Self {
            arguments,
            expression,
            constants,
        })
    }
    fn node(&self, e: &Expression, bindings: &[(String, Value)]) -> Result<Value, DeseqError> {
        match e {
            Expression::Literal(v) => Ok(v.clone()),
            Expression::Symbol(s) => bindings
                .iter()
                .find(|(n, _)| n == s)
                .map(|(_, v)| v.clone())
                .ok_or_else(|| invalid("uncaptured likelihood variable")),
            Expression::Call {
                function,
                arguments,
            } => {
                let is_ifelse = function == "ifelse" || function == "base::ifelse";
                let args = if is_ifelse {
                    let inputs = arguments
                        .iter()
                        .map(|a| (a.name.clone(), &a.expression))
                        .collect::<Vec<_>>();
                    let b = expression::arguments(&inputs, &["test", "yes", "no"])?;
                    let test = self.node(
                        b[0].ok_or_else(|| invalid("ifelse requires test"))?,
                        bindings,
                    )?;
                    let conditions = numeric(&test)?;
                    let yes = conditions.iter().any(|v| !v.is_nan() && *v != 0.0);
                    let no = conditions.contains(&0.0);
                    vec![
                        (None, test),
                        (
                            None,
                            if yes {
                                self.node(
                                    b[1].ok_or_else(|| invalid("ifelse requires yes"))?,
                                    bindings,
                                )?
                            } else {
                                Value::numbers(vec![0.0])
                            },
                        ),
                        (
                            None,
                            if no {
                                self.node(
                                    b[2].ok_or_else(|| invalid("ifelse requires no"))?,
                                    bindings,
                                )?
                            } else {
                                Value::numbers(vec![0.0])
                            },
                        ),
                    ]
                } else {
                    arguments
                        .iter()
                        .map(|a| {
                            self.node(&a.expression, bindings)
                                .map(|v| (a.name.clone(), v))
                        })
                        .collect::<Result<Vec<_>, _>>()?
                };
                match special(function) {
                    Some("%*%") => product(&args[0].1, &args[1].1),
                    Some("lgamma") => {
                        let b = expression::arguments(&args, &["x"])?;
                        let input = b[0].ok_or_else(|| invalid("lgamma requires x"))?;
                        let mut out = input.clone();
                        out.data = Data::Number {
                            values: numeric(input)?.into_iter().map(libm::lgamma).collect(),
                            integer: false,
                        };
                        Ok(out)
                    }
                    Some(name @ ("dnorm" | "dnbinom" | "dbinom" | "dt")) => {
                        distribution(name, &args)
                    }
                    _ => expression::evaluate(
                        &Expression::Call {
                            function: function.clone(),
                            arguments: args
                                .into_iter()
                                .map(|(name, v)| Argument {
                                    name,
                                    expression: Expression::Literal(v),
                                })
                                .collect(),
                        },
                        &[],
                        &EvaluationOptions {
                            collation: "C".into(),
                        },
                    ),
                }
            }
        }
    }
    pub fn evaluate(
        &self,
        y: &[f64],
        design: MatrixView<'_, f64>,
        beta: &[f64],
        parameters: Option<&[f64]>,
        offset: &[f64],
    ) -> Result<Vec<f64>, DeseqError> {
        if design.n_rows() != y.len() || design.n_cols() != beta.len() || offset.len() != y.len() {
            return Err(invalid("likelihood row shapes do not match"));
        }
        let mut x = Value::numbers(
            (0..design.n_rows())
                .flat_map(|i| (0..design.n_cols()).map(move |j| *design.at(i, j)))
                .collect(),
        );
        x.dimensions = Some([design.n_rows(), design.n_cols()]);
        let values = [
            Value::numbers(y.to_vec()),
            x,
            Value::numbers(beta.to_vec()),
            parameters
                .map(|v| Value::numbers(v.to_vec()))
                .unwrap_or_else(|| Value::new(Data::Null)),
            Value::numbers(offset.to_vec()),
        ];
        let bindings = self
            .arguments
            .iter()
            .cloned()
            .zip(values)
            .chain(self.constants.clone())
            .collect::<Vec<_>>();
        let result = numeric(&self.node(&self.expression, &bindings)?)?;
        if result.len() != y.len() {
            return Err(invalid("likelihood must return one value per observation"));
        }
        Ok(result)
    }
}
