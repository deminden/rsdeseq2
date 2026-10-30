// A factor in an interaction is contrasted only when the complementary term
// is present. This metadata rule is shared by both compatibility layouts.
fn formula_term_variables<'a>(state: &'a ExpandedFormulaDesignState<'_>) -> Vec<Vec<&'a str>> {
    let mut terms = Vec::new();
    terms.extend(state.selected_factors.iter().map(|f| vec![f.factor]));
    terms.extend(
        state
            .selected_numeric_covariates
            .iter()
            .map(|v| vec![v.name]),
    );
    terms.extend(
        state
            .factor_interactions
            .iter()
            .map(|t| vec![t.left_factor, t.right_factor]),
    );
    terms.extend(
        state
            .factor_numeric_interactions
            .iter()
            .map(|t| vec![t.factor, t.numeric]),
    );
    terms.extend(
        state
            .numeric_interactions
            .iter()
            .map(|t| vec![t.left_numeric, t.right_numeric]),
    );
    terms.extend(
        state
            .higher_order_interactions
            .iter()
            .map(|t| t.variables.iter().map(FormulaVariableRef::name).collect()),
    );
    terms
}

fn formula_factor_has_complement(
    state: &ExpandedFormulaDesignState<'_>,
    factor: &str,
    variables: &[&str],
) -> bool {
    formula_term_variables(state).iter().any(|term| {
        term.len() + 1 == variables.len()
            && !term.contains(&factor)
            && variables
                .iter()
                .filter(|&&v| v != factor)
                .all(|v| term.contains(v))
    })
}
