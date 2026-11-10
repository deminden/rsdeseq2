//! Restricted numeric functions imported as data, never foreign callbacks.
use super::{expression_quadrature::integrate, location_missing::NA};
use crate::{
    DeseqError,
    design::expression::{Data, EvaluationOptions, Expression, Value, evaluate},
    math::elementary::language_exp,
};

/// An elementwise arithmetic expression with scalar captured constants.
///
/// Only pure arithmetic, comparisons, elementary math, `ifelse`, `pmin` and
/// `pmax` are accepted. Arrays, reductions, assignments and unknown calls are
/// rejected during construction, including in unused branches.
#[derive(Clone, Debug)]
pub struct RestrictedFunction {
    arguments: Vec<String>,
    expression: Expression,
    constants: Vec<(String, Value)>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("restricted likelihood expression: {reason}"),
    }
}
pub(super) fn scalar(value: &Value) -> Result<f64, DeseqError> {
    if value.dimensions.is_some() || !value.attributes.is_empty() {
        return Err(invalid("only unclassed scalar values are supported"));
    }
    match &value.data {
        Data::Number { values, .. } => match values.as_slice() {
            [v] => Ok(*v),
            _ => Err(invalid("numeric value must have length one")),
        },
        Data::Logical(v) => match v.as_slice() {
            [Some(v)] => Ok(f64::from(*v)),
            [None] => Ok(NA),
            _ => Err(invalid("logical value must have length one")),
        },
        _ => Err(invalid("only numeric and logical values are supported")),
    }
}
pub(super) fn validate(
    expression: &Expression,
    bindings: &[&str],
    depth: usize,
    nodes: &mut usize,
) -> Result<(), DeseqError> {
    *nodes += 1;
    if depth > 128 || *nodes > 4096 {
        return Err(invalid("expression exceeds nesting or node limit"));
    }
    match expression {
        Expression::Literal(v) => {
            scalar(v)?;
        }
        Expression::Symbol(s) => {
            if !bindings.contains(&s.as_str()) {
                return Err(invalid("expression refers to an uncaptured variable"));
            }
        }
        Expression::Call {
            function,
            arguments,
        } => {
            let name = if let Some(name) = function.strip_prefix("base::") {
                name
            } else {
                function.as_str()
            };
            let (minimum, maximum, formals): (_, _, &[&str]) = match name {
                "(" | "I" | "abs" | "sqrt" | "exp" | "expm1" | "log1p" | "log2" | "log10"
                | "sign" | "floor" | "ceiling" | "trunc" | "!" => (1, 1, &["x"]),
                "log" => (1, 2, &["x", "base"]),
                "+" | "-" => (1, 2, &[]),
                "*" | "/" | "^" | "%%" | "%/%" | "<" | "<=" | ">" | ">=" | "==" | "!=" | "&"
                | "|" => (2, 2, &[]),
                "ifelse" => (3, 3, &["test", "yes", "no"]),
                "pmin" | "pmax" => (1, usize::MAX, &[]),
                _ => {
                    return Err(invalid(
                        "unsupported call; only pure elementwise arithmetic is allowed",
                    ));
                }
            };
            if arguments.len() < minimum || arguments.len() > maximum {
                return Err(invalid("unsupported function arity"));
            }
            let mut seen = vec![];
            for a in arguments {
                if let Some(name) = &a.name {
                    if !formals.contains(&name.as_str()) || seen.contains(&name.as_str()) {
                        return Err(invalid("unsupported or duplicate argument name"));
                    }
                    seen.push(name.as_str());
                }
                validate(&a.expression, bindings, depth + 1, nodes)?;
            }
        }
    }
    Ok(())
}
impl RestrictedFunction {
    /// Validate syntax, bindings and scalar types without executing foreign code.
    pub fn new(
        arguments: Vec<String>,
        expression: Expression,
        constants: Vec<(String, Value)>,
    ) -> Result<Self, DeseqError> {
        if arguments.is_empty() || arguments.len() > 2 {
            return Err(invalid("one or two scalar arguments are required"));
        }
        let mut bindings = vec![];
        for name in &arguments {
            if name.is_empty() || bindings.contains(&name.as_str()) {
                return Err(invalid("argument names must be nonempty and unique"));
            }
            bindings.push(name.as_str());
        }
        for (name, value) in &constants {
            if name.is_empty() || bindings.contains(&name.as_str()) {
                return Err(invalid("captured names must be nonempty and unique"));
            }
            scalar(value)?;
            bindings.push(name.as_str());
        }
        validate(&expression, &bindings, 0, &mut 0)?;
        let function = Self {
            arguments,
            expression,
            constants,
        };
        // All accepted primitives are total on scalar real inputs (domain
        // failures produce missing values). This also checks named binding.
        function.evaluate(&vec![1.0; function.arguments.len()])?;
        Ok(function)
    }
    /// Number of declared scalar arguments.
    pub fn arity(&self) -> usize {
        self.arguments.len()
    }
    /// Evaluate with the same standalone expression engine used by formulas.
    pub fn evaluate(&self, arguments: &[f64]) -> Result<f64, DeseqError> {
        if arguments.len() != self.arguments.len() {
            return Err(invalid("argument count differs from the imported function"));
        }
        let mut bindings = self.constants.clone();
        bindings.extend(
            self.arguments
                .iter()
                .zip(arguments)
                .map(|(name, &x)| (name.clone(), Value::numbers(vec![x]))),
        );
        scalar(&evaluate(
            &self.expression,
            &bindings,
            &EvaluationOptions::default(),
        )?)
    }
    fn value(&self, arguments: &[f64]) -> f64 {
        // Construction proves scalar types and call structure. Domain failures
        // remain numeric missing values; no runtime function lookup is possible.
        self.evaluate(arguments).unwrap_or(f64::NAN)
    }
}
/// Standalone custom likelihood leaves. Missing second-moment expressions retain
/// the original helper's missing SD behavior. Missing first moments use native
/// finite/infinite-interval quadrature of the supplied density expression.
#[derive(Clone, Debug)]
pub struct RestrictedLikelihood {
    pub(super) name: String,
    log_cdf: RestrictedFunction,
    log_density: RestrictedFunction,
    first: Option<RestrictedFunction>,
    second: Option<RestrictedFunction>,
}
impl RestrictedLikelihood {
    /// Build a likelihood from one-argument log-CDF/log-density functions and
    /// optional two-argument truncated first/second moment functions.
    pub fn new(
        name: String,
        log_cdf: RestrictedFunction,
        log_density: RestrictedFunction,
        first: Option<RestrictedFunction>,
        second: Option<RestrictedFunction>,
    ) -> Result<Self, DeseqError> {
        if name.is_empty()
            || log_cdf.arity() != 1
            || log_density.arity() != 1
            || first.as_ref().is_some_and(|f| f.arity() != 2)
            || second.as_ref().is_some_and(|f| f.arity() != 2)
        {
            return Err(invalid("likelihood name and function arities are invalid"));
        }
        Ok(Self {
            name,
            log_cdf,
            log_density,
            first,
            second,
        })
    }
    /// Declared likelihood name, retained as metadata for method validation.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Evaluate its supplied log-CDF expression in Rust.
    pub fn log_cdf(&self, x: f64) -> f64 {
        self.log_cdf.value(&[x])
    }
    /// Evaluate its supplied log-density expression in Rust.
    pub fn log_density(&self, x: f64) -> f64 {
        self.log_density.value(&[x])
    }
    /// First or second truncated moment; unsupported orders are errors.
    pub fn truncated_moment(&self, a: f64, b: f64, order: u8) -> Result<f64, DeseqError> {
        let expression = match order {
            1 => self.first.as_ref(),
            2 => self.second.as_ref(),
            _ => return Err(invalid("only first and second moments are supported")),
        };
        if let Some(f) = expression {
            return f.evaluate(&[a, b]);
        }
        if order == 2 {
            return Ok(NA);
        }
        if a == b {
            return Ok(a);
        }
        let denominator = language_exp(self.log_cdf(b)) - language_exp(self.log_cdf(a));
        if denominator.is_nan() {
            return Err(invalid("truncated CDF difference is missing"));
        }
        if denominator == 0.0 {
            return Ok(NA);
        }
        Ok(integrate(|x| x * language_exp(self.log_density(x)), a, b)
            .map(|v| v / denominator)
            .unwrap_or(NA))
    }
}

fn numeric_values(value: &Value) -> Vec<f64> {
    match &value.data {
        Data::Number { values, .. } => values.clone(),
        Data::Logical(v) => v.iter().map(|v| v.map(f64::from).unwrap_or(NA)).collect(),
        _ => vec![],
    }
}
fn expression_warnings(
    expression: &Expression,
    bindings: &[(String, Value)],
    warnings: &mut Vec<String>,
) {
    let Expression::Call {
        function,
        arguments,
    } = expression
    else {
        return;
    };
    let name = function.strip_prefix("base::").unwrap_or(function);
    if name == "ifelse" {
        let mut ordered = [None, None, None];
        for argument in arguments {
            if let Some(name) = argument.name.as_deref() {
                let i = match name {
                    "test" => 0,
                    "yes" => 1,
                    _ => 2,
                };
                ordered[i] = Some(&argument.expression);
            }
        }
        for argument in arguments {
            if argument.name.is_none()
                && let Some(slot) = ordered.iter_mut().find(|v| v.is_none())
            {
                *slot = Some(&argument.expression);
            }
        }
        if let [Some(test), Some(yes), Some(no)] = ordered {
            expression_warnings(test, bindings, warnings);
            if let Ok(value) = evaluate(test, bindings, &EvaluationOptions::default()) {
                let values = numeric_values(&value);
                if values.iter().any(|v| !v.is_nan() && *v != 0.0) {
                    expression_warnings(yes, bindings, warnings);
                }
                if values.contains(&0.0) {
                    expression_warnings(no, bindings, warnings);
                }
            }
        }
        return;
    }
    for argument in arguments {
        expression_warnings(&argument.expression, bindings, warnings);
    }
    if !matches!(name, "log" | "log1p" | "log2" | "log10" | "sqrt") {
        return;
    }
    let Ok(result) = evaluate(expression, bindings, &EvaluationOptions::default()) else {
        return;
    };
    let output = numeric_values(&result);
    let inputs = arguments
        .iter()
        .filter_map(|a| {
            evaluate(&a.expression, bindings, &EvaluationOptions::default())
                .ok()
                .map(|v| numeric_values(&v))
        })
        .collect::<Vec<_>>();
    if output.iter().enumerate().any(|(i, v)| {
        v.is_nan()
            && inputs
                .iter()
                .all(|a| !a.is_empty() && !a[i % a.len()].is_nan())
    }) {
        warnings.push("NaNs produced".into());
    }
}
impl RestrictedFunction {
    fn warnings(&self, arguments: &[&[f64]]) -> Vec<String> {
        let mut bindings = self.constants.clone();
        bindings.extend(
            self.arguments
                .iter()
                .zip(arguments)
                .map(|(n, v)| (n.clone(), Value::numbers(v.to_vec()))),
        );
        let mut warnings = vec![];
        expression_warnings(&self.expression, &bindings, &mut warnings);
        warnings
    }
}
impl RestrictedLikelihood {
    pub(super) fn cdf_warnings(&self, x: &[f64]) -> Vec<String> {
        self.log_cdf.warnings(&[x])
    }
    pub(super) fn density_warnings(&self, x: &[f64]) -> Vec<String> {
        self.log_density.warnings(&[x])
    }
    pub(super) fn moment_warnings(&self, a: &[f64], b: &[f64], order: u8) -> Vec<String> {
        let f = if order == 1 {
            self.first.as_ref()
        } else {
            self.second.as_ref()
        };
        f.map(|f| f.warnings(&[a, b])).unwrap_or_default()
    }
}
