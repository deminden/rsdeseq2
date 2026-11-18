use super::*;
// Private state and derived-covariate types shared by the formula implementation.

pub(super) struct ExpandedFormulaDesignState<'a> {
    pub(super) frame_order: Vec<String>,
    pub(super) term_order: Option<Vec<Vec<String>>>,
    pub(super) has_intercept: bool,
    pub(super) selected_factors: Vec<ExpandedFactorSpec<'a>>,
    pub(super) selected_numeric_covariates: Vec<ExpandedNumericSpec<'a>>,
    pub(super) factor_interactions: Vec<ExpandedFactorInteractionSpec<'a>>,
    pub(super) factor_numeric_interactions: Vec<ExpandedFactorNumericInteractionSpec<'a>>,
    pub(super) numeric_interactions: Vec<ExpandedNumericInteractionSpec<'a>>,
    pub(super) higher_order_interactions: Vec<FormulaHigherOrderInteractionSpec<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct FormulaHigherOrderInteractionSpec<'a> {
    pub(super) variables: Vec<FormulaVariableRef<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct FormulaDerivedNumericCovariate {
    pub(super) name: String,
    pub(super) values: Vec<f64>,
    pub(super) group: Option<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct FormulaDerivedFactorCovariate {
    pub(super) name: String,
    pub(super) sample_levels: Vec<String>,
    pub(super) reference: String,
    pub(super) levels: Option<Vec<String>>,
    pub(super) primitive: super::super::expression::Value,
}

pub(super) type FormulaDerivedNumericTerm = (String, Vec<f64>);
pub(super) type FormulaNumericTransformExpansion = (String, Vec<FormulaDerivedNumericTerm>);

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum FormulaScaleOption {
    Auto,
    Disabled,
    Explicit(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum FormulaVariableRef<'a> {
    Factor(&'a str),
    Numeric(&'a str),
}

impl FormulaVariableRef<'_> {
    pub(super) fn name(&self) -> &str {
        match self {
            Self::Factor(name) | Self::Numeric(name) => name,
        }
    }
}

impl Default for ExpandedFormulaDesignState<'_> {
    fn default() -> Self {
        Self {
            frame_order: Vec::new(),
            term_order: None,
            has_intercept: true,
            selected_factors: Vec::new(),
            selected_numeric_covariates: Vec::new(),
            factor_interactions: Vec::new(),
            factor_numeric_interactions: Vec::new(),
            numeric_interactions: Vec::new(),
            higher_order_interactions: Vec::new(),
        }
    }
}
