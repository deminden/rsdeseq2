// Private state and derived-covariate types shared by the formula implementation.

struct ExpandedFormulaDesignState<'a> {
    frame_order: Vec<String>,
    has_intercept: bool,
    selected_factors: Vec<ExpandedFactorSpec<'a>>,
    selected_numeric_covariates: Vec<ExpandedNumericSpec<'a>>,
    factor_interactions: Vec<ExpandedFactorInteractionSpec<'a>>,
    factor_numeric_interactions: Vec<ExpandedFactorNumericInteractionSpec<'a>>,
    numeric_interactions: Vec<ExpandedNumericInteractionSpec<'a>>,
    higher_order_interactions: Vec<FormulaHigherOrderInteractionSpec<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
struct FormulaHigherOrderInteractionSpec<'a> {
    variables: Vec<FormulaVariableRef<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
struct FormulaDerivedNumericCovariate {
    name: String,
    values: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
struct FormulaDerivedFactorCovariate {
    name: String,
    sample_levels: Vec<String>,
    reference: String,
    levels: Option<Vec<String>>,
    primitive: super::expression::Value,
}

type FormulaDerivedNumericTerm = (String, Vec<f64>);
type FormulaNumericTransformExpansion = (String, Vec<FormulaDerivedNumericTerm>);

#[derive(Clone, Copy, Debug, PartialEq)]
enum FormulaScaleOption {
    Auto,
    Disabled,
    Explicit(f64),
}

#[derive(Clone, Debug, PartialEq)]
enum FormulaVariableRef<'a> {
    Factor(&'a str),
    Numeric(&'a str),
}

impl FormulaVariableRef<'_> {
    fn name(&self) -> &str {
        match self {
            Self::Factor(name) | Self::Numeric(name) => name,
        }
    }
}

impl Default for ExpandedFormulaDesignState<'_> {
    fn default() -> Self {
        Self {
            frame_order: Vec::new(),
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
