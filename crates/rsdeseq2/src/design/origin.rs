//! Formula provenance affects result policy without changing a numeric design.
use crate::errors::DeseqError;
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq, Eq)]
enum Source {
    Matrix,
    Formula(Arc<[String]>),
    Legacy,
}

/// Origin of a fitted model matrix, separate from its numeric values.
/// Formula variables describe the original expression, including removed terms
/// and offsets. Primitive matrices do not implicitly acquire formula semantics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DesignOrigin(Source);
impl Default for DesignOrigin {
    fn default() -> Self {
        Self::matrix()
    }
}
impl DesignOrigin {
    /// A directly supplied matrix, with no formula-only result exceptions.
    pub fn matrix() -> Self {
        Self(Source::Matrix)
    }
    /// Explicit compatibility context for older stored fits without provenance.
    /// This preserves their former partition-based inference and does not prove
    /// that they originated from a particular formula.
    pub fn legacy() -> Self {
        Self(Source::Legacy)
    }
    /// Store identifiers extracted from a formula, deduplicated in occurrence
    /// order as in R's `all.vars()`. Empty variables describe an intercept-only
    /// formula; empty identifiers are invalid.
    pub fn formula(variables: Vec<String>) -> Result<Self, DeseqError> {
        let mut unique = Vec::with_capacity(variables.len());
        for variable in variables {
            if variable.is_empty() || variable.contains('\0') {
                return Err(DeseqError::InvalidOptions {
                    reason: "formula origin variables must be nonempty identifiers without NUL"
                        .into(),
                });
            }
            if !unique.contains(&variable) {
                unique.push(variable);
            }
        }
        Ok(Self(Source::Formula(unique.into())))
    }
    /// Original formula identifiers, or None for matrix/legacy origins.
    pub fn variables(&self) -> Option<&[String]> {
        match &self.0 {
            Source::Formula(v) => Some(v),
            _ => None,
        }
    }
    /// Whether this is an explicitly unverified legacy context.
    pub fn is_legacy(&self) -> bool {
        matches!(self.0, Source::Legacy)
    }
}
