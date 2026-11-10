//! Standalone evaluation of built-in formula expressions. No R execution or
//! user-defined function callback is part of this numerical interface.
use crate::errors::DeseqError;
use icu_collator::{Collator, CollatorBorrowed};
use icu_locale::Locale;

mod arithmetic;
mod comparisons;
mod conditional;
mod extrema;
mod factors;
pub mod frame;
mod functions;
mod polynomials;
mod shape;
mod summaries;
mod transforms;
mod validation;
pub mod wire;
use validation::{function_name, validate_expression, validate_value};

/// Atomic/list storage. Numeric matrices use row-major order; missing logical
/// and text values remain distinct from empty values.
#[derive(Clone, Debug, PartialEq)]
pub enum Data {
    /// Double values, or integer values represented exactly as doubles.
    Number { values: Vec<f64>, integer: bool },
    /// Logical values with explicit missingness.
    Logical(Vec<Option<bool>>),
    /// UTF-8 text with explicit missingness.
    Text(Vec<Option<String>>),
    /// Nested primitive lists, such as polynomial recurrence coefficients.
    List(Vec<Value>),
    /// R's absent argument value.
    Null,
}
/// A formula value with interface metadata. Attributes contain only primitive
/// values; they never retain foreign runtime pointers or executable closures.
#[derive(Clone, Debug, PartialEq)]
pub struct Value {
    /// Stored values.
    pub data: Data,
    /// Optional matrix shape.
    pub dimensions: Option<[usize; 2]>,
    /// Optional element/list names in storage order (row-major for matrices).
    pub names: Option<Vec<String>>,
    /// Optional row/column names.
    pub dimnames: [Option<Vec<String>>; 2],
    /// Optional names attached to the row/column label vectors themselves.
    /// Multivariate polynomial column labels use these grid-row identifiers.
    pub dimname_names: [Option<Vec<String>>; 2],
    /// Class, factor levels, recurrence/scaling and other primitive metadata.
    pub attributes: Vec<(String, Value)>,
}
impl Value {
    /// Construct an ordinary double vector.
    pub fn numbers(values: Vec<f64>) -> Self {
        Self::new(Data::Number {
            values,
            integer: false,
        })
    }
    /// Construct a scalar.
    pub fn scalar(value: f64) -> Self {
        Self::numbers(vec![value])
    }
    /// Construct a value without metadata.
    pub fn new(mut data: Data) -> Self {
        if let Data::Number {
            values,
            integer: true,
        } = &mut data
        {
            for v in values {
                *v = integer_storage(*v);
            }
        }
        Self {
            data,
            dimensions: None,
            names: None,
            dimnames: [None, None],
            dimname_names: [None, None],
            attributes: Vec::new(),
        }
    }
    /// Number of stored elements.
    pub fn len(&self) -> usize {
        match &self.data {
            Data::Number { values, .. } => values.len(),
            Data::Logical(v) => v.len(),
            Data::Text(v) => v.len(),
            Data::List(v) => v.len(),
            Data::Null => 0,
        }
    }
    /// Whether storage is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Return a primitive attribute.
    pub fn attribute(&self, name: &str) -> Option<&Value> {
        self.attributes
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v)
    }
    /// Replace an attribute, preserving the order of the remaining metadata.
    pub fn set_attribute(&mut self, name: &str, value: Value) {
        self.attributes.retain(|(n, _)| n != name);
        self.attributes.push((name.to_owned(), value));
    }
    pub(super) fn numeric(&self) -> Result<Vec<f64>, DeseqError> {
        match &self.data {
            Data::Number { values, .. } => Ok(values.clone()),
            Data::Logical(v) => Ok(v
                .iter()
                .map(|v| v.map_or(NA_REAL, |b| f64::from(u8::from(b))))
                .collect()),
            Data::Null => Ok(Vec::new()),
            _ => Err(invalid("formula operation requires numeric values")),
        }
    }
    pub(super) fn scalar_number(&self) -> Result<f64, DeseqError> {
        let x = self.numeric()?;
        if x.len() != 1 {
            return Err(invalid("formula argument must be a numeric scalar"));
        };
        Ok(x[0])
    }
    pub(super) fn bool_scalar(&self) -> Result<bool, DeseqError> {
        let v = self.scalar_number()?;
        if !v.is_finite() {
            return Err(invalid("formula flag must be nonmissing"));
        };
        Ok(v != 0.0)
    }
    pub(super) fn text(&self) -> Result<Vec<Option<String>>, DeseqError> {
        match &self.data {
            Data::Text(v) => Ok(v.clone()),
            _ => Err(invalid("formula argument requires text")),
        }
    }
    pub(super) fn has_class(&self, class: &str) -> bool {
        matches!(self.attribute("class").map(|v|&v.data),Some(Data::Text(v)) if v.iter().any(|s|s.as_deref()==Some(class)))
    }
    pub(super) fn with_numbers(&self, mut values: Vec<f64>, integer: bool) -> Self {
        if integer {
            for v in &mut values {
                *v = integer_storage(*v);
            }
        }
        let mut v = self.clone();
        v.data = Data::Number { values, integer };
        v
    }
}
/// A built-in expression or a raw data/environment binding.
#[derive(Clone, Debug, PartialEq)]
pub enum Expression {
    /// Already supplied primitive value.
    Literal(Value),
    /// Exact variable lookup, before any name cleanup.
    Symbol(String),
    /// Built-in function/operator with positional or named arguments.
    Call {
        function: String,
        arguments: Vec<Argument>,
    },
}
/// One call argument.
#[derive(Clone, Debug, PartialEq)]
pub struct Argument {
    /// An optional formal argument name.
    pub name: Option<String>,
    /// Value-producing expression.
    pub expression: Expression,
}
/// Locale metadata for categorical ordering; numeric operations do not consult
/// the calling frontend or host R runtime.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaluationOptions {
    /// ICU locale. "root" and "und" select Unicode root collation; "C" uses
    /// byte order. Other locale identifiers are parsed as Unicode locales.
    pub collation: String,
}
impl Default for EvaluationOptions {
    fn default() -> Self {
        Self {
            collation: "root".into(),
        }
    }
}
pub(super) fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}
// Prevent short symbolic requests (e.g. seq_len(INT_MAX)) from requesting
// unbounded intermediate storage before the bounded transport can reject it.
pub(super) fn checked_elements(count: Option<usize>) -> Result<usize, DeseqError> {
    count
        .filter(|n| *n <= wire::MAX_BYTES / std::mem::size_of::<f64>())
        .ok_or_else(|| invalid("formula expression output exceeds size limit"))
}
pub(super) const NA_REAL: f64 = f64::from_bits(0x7ff0_0000_0000_07a2);
// Integer storage has one missing sentinel and no signed zero. Canonicalize
// its double carrier without changing any ordinary floating-point value.
pub(super) fn integer_storage(value: f64) -> f64 {
    if value.is_nan() {
        NA_REAL
    } else if value == 0.0 {
        0.0
    } else {
        value
    }
}
pub(super) fn is_na(value: f64) -> bool {
    value.is_nan() && value.to_bits() as u32 == 1954
}
pub(super) struct Evaluator<'a> {
    variables: &'a [(String, Value)],
    collator: Option<CollatorBorrowed<'static>>,
}
/// Evaluate one expression using the standalone built-in engine. Unknown/custom
/// functions are rejected and never dispatched to another language runtime.
pub fn evaluate(
    expression: &Expression,
    variables: &[(String, Value)],
    options: &EvaluationOptions,
) -> Result<Value, DeseqError> {
    evaluate_many(std::slice::from_ref(expression), variables, options).map(|mut v| v.remove(0))
}

/// Evaluate a bounded set with shared bindings and collation initialization.
pub fn evaluate_many(
    expressions: &[Expression],
    variables: &[(String, Value)],
    options: &EvaluationOptions,
) -> Result<Vec<Value>, DeseqError> {
    for (index, (name, value)) in variables.iter().enumerate() {
        if name.is_empty() || variables[..index].iter().any(|(n, _)| n == name) {
            return Err(invalid("formula bindings require unique nonempty names"));
        }
        validate_value(value, 0)?;
    }
    for expression in expressions {
        validate_expression(expression, 0)?;
    }
    let collation = if matches!(options.collation.as_str(), "root" | "") {
        "und"
    } else {
        &options.collation
    };
    let collator = if collation == "C" {
        None
    } else {
        Some(
            Collator::try_new(
                collation
                    .parse::<Locale>()
                    .map_err(|e| invalid(e.to_string()))?
                    .into(),
                Default::default(),
            )
            .map_err(|e| invalid(e.to_string()))?,
        )
    };
    let mut evaluator = Evaluator {
        variables,
        collator,
    };
    expressions
        .iter()
        .map(|e| evaluator.expression(e, 0))
        .collect()
}
impl Evaluator<'_> {
    fn expression(&mut self, expression: &Expression, depth: usize) -> Result<Value, DeseqError> {
        if depth > 128 {
            return Err(invalid("formula expression nesting exceeds 128"));
        }
        match expression {
            Expression::Literal(v) => {
                validate_value(v, depth)?;
                Ok(v.clone())
            }
            Expression::Symbol(name) => self
                .variables
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, v)| v.clone())
                .ok_or_else(|| invalid(format!("unknown formula variable '{name}'"))),
            Expression::Call {
                function,
                arguments,
            } => {
                let name = function_name(function)?;
                let value = if name == "ifelse" {
                    self.conditional(arguments, depth)?
                } else {
                    let mut args = Vec::with_capacity(arguments.len());
                    for a in arguments {
                        args.push((a.name.clone(), self.expression(&a.expression, depth + 1)?));
                    }
                    self.call(name, &args, arguments)?
                };
                // A built-in can construct shape/name metadata too. Validate
                // it before a parent expression performs indexed operations.
                validate_value(&value, depth)?;
                Ok(value)
            }
        }
    }
}

pub(crate) fn arguments<'a, T>(
    args: &'a [(Option<String>, T)],
    formals: &[&str],
) -> Result<Vec<Option<&'a T>>, DeseqError> {
    let mut bound = vec![None; formals.len()];
    for (name, value) in args.iter().filter(|(name, _)| name.is_some()) {
        let name = name.as_deref().unwrap();
        let matches = formals
            .iter()
            .enumerate()
            .filter(|(_, f)| f.starts_with(name))
            .collect::<Vec<_>>();
        let i = formals
            .iter()
            .position(|f| *f == name)
            .or_else(|| (matches.len() == 1).then(|| matches[0].0))
            .ok_or_else(|| invalid(format!("unknown or ambiguous formula argument '{name}'")))?;
        if bound[i].replace(value).is_some() {
            return Err(invalid(format!("duplicate formula argument '{name}'")));
        }
    }
    for (_, value) in args.iter().filter(|(name, _)| name.is_none()) {
        let i = bound
            .iter()
            .position(Option::is_none)
            .ok_or_else(|| invalid("too many formula arguments"))?;
        bound[i] = Some(value);
    }
    Ok(bound)
}
