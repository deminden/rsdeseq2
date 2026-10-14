/// Build an expanded additive design from a primitive DESeq2-style formula subset.
///
/// Supported right-hand-side terms are intercept-preserving main effects
/// (`condition`, `dose`), intercept-only `1`, pairwise interactions
/// (`condition:dose`), nested shorthand (`condition/batch`), R nesting
/// operator terms (`batch %in% condition`), and `*` shorthand for main effects
/// plus interactions (`condition*dose`).
/// Interaction variables can appear without corresponding main effects. The
/// reported standard-design interaction columns then follow R model-matrix
/// treatment-coding shape for the supported primitive terms. Intercept removal
/// with `0` or `-1` and R-style intercept restoration with `+1` or `-0` are
/// supported for these primitive terms. Formula interactions can contain two
/// or more variables. Primitive `- term` subtraction, signed additive
/// parenthesized groups such as `(condition + batch - batch)`, and R-like
/// unary `+` / `-` sign sequences are supported for the same term subset.
/// Signed additive groups are reduced inside supported `*`, `:`, `/`, `%in%`,
/// and power expressions before operator products are formed.
/// Plain `I(numeric)` and signed `I(-numeric)` terms, simple
/// `I(numeric op scalar)` arithmetic, and integer numeric power transforms,
/// raw polynomial transforms, and default orthogonal polynomial transforms are
/// materialized as derived numeric covariates. Supported
/// parenthesized formula powers such as `(condition + batch + dose)^2` expand
/// into main effects plus interactions up to the requested order, with signed
/// additive bases normalized before expansion. Dot formula terms such as `.`,
/// `condition:.`, `condition*.`, and `.^k` use the supplied factors and
/// numeric covariates as the base variable set.
/// Supported
/// `offset(numeric)` and single-vector `offset(transform(numeric))` terms are
/// parsed by [`expanded_formula_design_with_offsets`]; this compatibility
/// helper returns only the design matrix.
pub fn expanded_formula_design<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    Ok(expanded_formula_design_with_offsets(formula, factors, numeric_covariates)?.design)
}

/// Build an expanded formula design and return supported per-sample offsets.
///
/// Supported offsets are `offset(numeric)` terms where `numeric` names a
/// supplied numeric covariate, plus single-vector offsets from the supported
/// numeric transform subset such as `offset(log2(numeric))` or
/// `offset(I(numeric + other_numeric))`. Multiple offset terms are summed
/// sample-wise.
pub fn expanded_formula_design_with_offsets<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    let rhs = formula_rhs(formula)?;
    let (rhs, offsets) = extract_formula_offsets(rhs, factors, numeric_covariates)?;
    let (rhs, derived_factor_covariates) = expand_formula_factor_transform_terms(&rhs, factors)?;
    let mut all_factors =
        Vec::with_capacity(factors.len() + derived_factor_covariates.len());
    all_factors.extend(factors.iter().cloned());
    all_factors.extend(derived_factor_covariates.iter().map(|factor| ExpandedFactorSpec {
        factor: factor.name.as_str(),
        sample_levels: factor.sample_levels.as_slice(),
        reference: factor.reference.as_str(),
        levels: factor.levels.as_deref(),
    }));
    validate_formula_factor_name_collisions(&all_factors)?;
    let (rhs, derived_numeric_covariates) =
        expand_formula_numeric_transform_terms(&rhs, &all_factors, numeric_covariates)?;
    let mut all_numeric_covariates =
        Vec::with_capacity(numeric_covariates.len() + derived_numeric_covariates.len());
    all_numeric_covariates.extend(numeric_covariates.iter().cloned());
    all_numeric_covariates.extend(derived_numeric_covariates.iter().map(|covariate| {
        ExpandedNumericSpec {
            name: covariate.name.as_str(),
            values: covariate.values.as_slice(),
        }
    }));
    validate_formula_numeric_covariate_name_collisions(&all_numeric_covariates)?;
    validate_formula_factor_numeric_name_collisions(&all_factors, &all_numeric_covariates)?;
    let mut state = ExpandedFormulaDesignState::default();

    let rhs = expand_parenthesized_formula_terms(&rhs)?;
    for (sign, raw_term) in split_formula_signed_terms(&rhs)? {
        let term = raw_term.trim();
        if term.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula '{formula}' contains an empty term"),
            });
        }
        if sign < 0 {
            remove_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if term == "1" {
            state.has_intercept = true;
            continue;
        }
        if term == "0" || term == "-1" {
            state.has_intercept = false;
            continue;
        }
        if formula_contains_top_level(term, '^')? {
            add_power_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if formula_contains_top_level(term, '(')? || formula_contains_top_level(term, ')')? {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula term '{term}' is not supported by the primitive parser"),
            });
        }
        if formula_contains_top_level(term, '/')? {
            add_nested_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if formula_contains_top_level(term, '*')? {
            add_star_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if formula_contains_top_level(term, ':')? {
            add_interaction_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        add_main_formula_term(
            term,
            &all_factors,
            &all_numeric_covariates,
            &mut state.selected_factors,
            &mut state.selected_numeric_covariates,
        )?;
    }

    let model_frame = formula_model_frame_from_specs(&all_factors, &all_numeric_covariates);
    let design =
        expanded_formula_design_from_state(&state, &all_factors, &all_numeric_covariates)?;
    Ok(ExpandedFormulaDesignWithOffsets {
        design,
        offsets,
        model_frame,
    })
}

fn validate_formula_factor_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(), DeseqError> {
    for (idx, factor) in factors.iter().enumerate() {
        for previous in &factors[..idx] {
            if previous.factor == factor.factor {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", factor.factor),
                });
            }
            let previous_aliases = r_like_name_candidates(previous.factor);
            let current_aliases = r_like_name_candidates(factor.factor);
            if previous_aliases
                .iter()
                .any(|previous_alias| current_aliases.contains(previous_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        factor.factor
                    ),
                });
            }
        }
    }
    Ok(())
}

fn validate_formula_numeric_covariate_name_collisions(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), DeseqError> {
    for (idx, covariate) in numeric_covariates.iter().enumerate() {
        for previous in &numeric_covariates[..idx] {
            if previous.name == covariate.name {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", covariate.name),
                });
            }
            let previous_aliases = r_like_name_candidates(previous.name);
            let current_aliases = r_like_name_candidates(covariate.name);
            if previous_aliases
                .iter()
                .any(|previous_alias| current_aliases.contains(previous_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        covariate.name
                    ),
                });
            }
        }
    }
    Ok(())
}

fn validate_formula_factor_numeric_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), DeseqError> {
    for factor in factors {
        for covariate in numeric_covariates {
            if factor.factor == covariate.name {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", factor.factor),
                });
            }
            let factor_aliases = r_like_name_candidates(factor.factor);
            let covariate_aliases = r_like_name_candidates(covariate.name);
            if factor_aliases
                .iter()
                .any(|factor_alias| covariate_aliases.contains(factor_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        covariate.name
                    ),
                });
            }
        }
    }
    Ok(())
}

fn formula_model_frame_from_specs(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> FormulaModelFrame {
    FormulaModelFrame {
        factors: factors
            .iter()
            .map(|factor| FormulaFactorColumn {
                name: factor.factor.to_string(),
                sample_levels: factor.sample_levels.to_vec(),
                levels: factor.levels.map(<[String]>::to_vec),
                reference: Some(factor.reference.to_string()),
            })
            .collect(),
        numeric_covariates: numeric_covariates
            .iter()
            .map(|covariate| FormulaNumericColumn {
                name: covariate.name.to_string(),
                values: covariate.values.to_vec(),
            })
            .collect(),
    }
}

/// Return whether a primitive formula contains at least one `offset(...)` term.
///
/// This syntax-level helper validates the formula shape enough to distinguish
/// real offset terms from ordinary variable names before any model-frame lookup
/// is required.
pub fn formula_has_offset_terms(formula: &str) -> Result<bool, DeseqError> {
    let rhs = formula_rhs(formula)?;
    for (_, term) in split_formula_signed_terms(rhs)? {
        if formula_term_has_offset(&term)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Build an expanded formula design from owned model-frame metadata.
///
/// This is the wrapper/object-facing companion to [`expanded_formula_design`].
/// It derives borrowed factor/numeric specs from an owned
/// [`FormulaModelFrame`], inferring each factor reference from the first
/// declared factor level when available, otherwise the first observed sample
/// level, when the caller did not supply one.
pub fn expanded_formula_design_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    Ok(expanded_formula_design_with_offsets_from_model_frame(formula, model_frame)?.design)
}

/// Build an expanded formula design plus offsets from owned model-frame metadata.
pub fn expanded_formula_design_with_offsets_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    let resolved = ResolvedFormulaModelFrame::new(model_frame)?;
    expanded_formula_design_with_offsets(
        formula,
        &resolved.factor_specs,
        &resolved.numeric_specs,
    )
}

struct ResolvedFormulaModelFrame<'a> {
    factor_specs: Vec<ExpandedFactorSpec<'a>>,
    numeric_specs: Vec<ExpandedNumericSpec<'a>>,
}

impl<'a> ResolvedFormulaModelFrame<'a> {
    fn new(model_frame: &'a FormulaModelFrame) -> Result<Self, DeseqError> {
        validate_formula_model_frame(model_frame)?;
        let factor_specs = model_frame
            .factors
            .iter()
            .map(|factor| {
                let reference = formula_model_frame_factor_reference(factor)
                    .expect("model-frame factor reference validated before borrowing specs");
                ExpandedFactorSpec {
                    factor: factor.name.as_str(),
                    sample_levels: factor.sample_levels.as_slice(),
                    reference,
                    levels: factor.levels.as_deref(),
                }
            })
            .collect::<Vec<_>>();
        let numeric_specs = model_frame
            .numeric_covariates
            .iter()
            .map(|covariate| ExpandedNumericSpec {
                name: covariate.name.as_str(),
                values: covariate.values.as_slice(),
            })
            .collect::<Vec<_>>();
        Ok(Self {
            factor_specs,
            numeric_specs,
        })
    }
}

fn validate_formula_model_frame(model_frame: &FormulaModelFrame) -> Result<(), DeseqError> {
    let n_samples = formula_model_frame_sample_count(model_frame)?;
    for (idx, factor) in model_frame.factors.iter().enumerate() {
        validate_formula_model_frame_column_name(&factor.name)?;
        if model_frame.factors[..idx]
            .iter()
            .any(|previous| previous.name == factor.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula factor '{}' appears more than once", factor.name),
            });
        }
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "formula factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        let reference = formula_model_frame_factor_reference(factor)?;
        validate_factor_design_inputs_with_levels(
            &factor.name,
            &factor.sample_levels,
            reference,
            factor.levels.as_deref(),
        )?;
    }
    for (idx, covariate) in model_frame.numeric_covariates.iter().enumerate() {
        validate_formula_model_frame_column_name(&covariate.name)?;
        if model_frame.numeric_covariates[..idx]
            .iter()
            .any(|previous| previous.name == covariate.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' appears more than once",
                    covariate.name
                ),
            });
        }
        if model_frame
            .factors
            .iter()
            .any(|factor| factor.name == covariate.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula column '{}' cannot be both a factor and numeric covariate",
                    covariate.name
                ),
            });
        }
        if model_frame.factors.iter().any(|factor| {
            let factor_aliases = r_like_name_candidates(&factor.name);
            let covariate_aliases = r_like_name_candidates(&covariate.name);
            factor_aliases
                .iter()
                .any(|factor_alias| covariate_aliases.contains(factor_alias))
        }) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula column '{}' resolves ambiguously after R-style cleanup",
                    covariate.name
                ),
            });
        }
        if covariate.values.len() != n_samples {
            return Err(invalid_dimensions(
                "formula numeric covariate values",
                n_samples,
                covariate.values.len(),
            ));
        }
        if !covariate.values.iter().all(|value| value.is_finite()) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' contains non-finite values",
                    covariate.name
                ),
            });
        }
    }
    Ok(())
}

fn formula_model_frame_factor_reference(
    factor: &FormulaFactorColumn,
) -> Result<&str, DeseqError> {
    if let Some(reference) = factor.reference.as_deref() {
        return resolve_formula_model_frame_reference_alias(factor, reference);
    }
    factor
        .levels
        .as_deref()
        .and_then(|levels| levels.first().map(String::as_str))
        .or_else(|| factor.sample_levels.first().map(String::as_str))
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{}' requires at least one sample level",
                factor.name
            ),
        })
}

fn resolve_formula_model_frame_reference_alias<'a>(
    factor: &'a FormulaFactorColumn,
    requested: &str,
) -> Result<&'a str, DeseqError> {
    let levels = formula_model_frame_factor_level_candidates(factor);
    resolve_formula_reference_alias_from_levels(&factor.name, &levels, requested)
}

fn formula_model_frame_factor_level_candidates(factor: &FormulaFactorColumn) -> Vec<&str> {
    let mut levels = Vec::new();
    if let Some(declared) = &factor.levels {
        for level in declared {
            if !levels.iter().any(|candidate| candidate == &level.as_str()) {
                levels.push(level.as_str());
            }
        }
    }
    for level in &factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == &level.as_str()) {
            levels.push(level.as_str());
        }
    }
    levels
}

fn resolve_formula_reference_alias_from_levels<'a>(
    factor_name: &str,
    levels: &[&'a str],
    requested: &str,
) -> Result<&'a str, DeseqError> {
    let exact = levels
        .iter()
        .copied()
        .filter(|level| *level == requested)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [level] => return Ok(*level),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula factor '{factor_name}' reference level '{requested}' appears more than once"
                ),
            });
        }
    }

    let aliases = levels
        .iter()
        .copied()
        .filter(|level| {
            r_like_name_candidates(level)
                .into_iter()
                .any(|candidate| candidate == requested)
        })
        .collect::<Vec<_>>();
    match aliases.as_slice() {
        [level] => Ok(*level),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' does not contain reference level '{requested}'"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' reference level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}

fn formula_model_frame_sample_count(model_frame: &FormulaModelFrame) -> Result<usize, DeseqError> {
    if let Some(factor) = model_frame.factors.first() {
        if factor.sample_levels.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula factor '{}' requires at least one sample level",
                    factor.name
                ),
            });
        }
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = model_frame.numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' requires at least one sample value",
                    covariate.name
                ),
            });
        }
        return Ok(covariate.values.len());
    }
    Err(DeseqError::InvalidOptions {
        reason: "formula model frame requires at least one column".to_string(),
    })
}

struct ExpandedFormulaDesignState<'a> {
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

fn formula_rhs(formula: &str) -> Result<&str, DeseqError> {
    let trimmed = formula.trim();
    let Some(rhs) = trimmed.strip_prefix('~') else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula '{formula}' must start with '~'"),
        });
    };
    let rhs = rhs.trim();
    if rhs.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula right-hand side must be non-empty".to_string(),
        });
    }
    Ok(rhs)
}

fn extract_formula_offsets(
    rhs: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let signed_terms = split_formula_signed_terms(rhs)?;
    let n_samples = numeric_covariates
        .first()
        .map(|covariate| covariate.values.len())
        .unwrap_or(0);
    let mut offsets = vec![0.0; n_samples];
    let mut saw_offset = false;
    let mut remaining_terms = Vec::new();
    for (sign, term) in signed_terms {
        extract_formula_offset_term(
            sign,
            &term,
            factors,
            numeric_covariates,
            &mut offsets,
            &mut saw_offset,
            &mut remaining_terms,
        )?;
    }
    if remaining_terms.is_empty() {
        remaining_terms.push("1".to_string());
    }
    if !saw_offset {
        offsets.clear();
    }
    Ok((join_formula_terms(&remaining_terms), offsets))
}

fn formula_term_has_offset(term: &str) -> Result<bool, DeseqError> {
    let term = term.trim();
    if formula_offset_expression(term)?.is_some() {
        return Ok(true);
    }
    let Some(stripped) = formula_outer_group_for_offset(term)? else {
        return Ok(false);
    };
    for (_, expanded_term) in split_formula_signed_terms(stripped)? {
        if formula_term_has_offset(&expanded_term)? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn extract_formula_offset_term(
    sign: i8,
    term: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    offsets: &mut Vec<f64>,
    saw_offset: &mut bool,
    remaining_terms: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let term = term.trim();
    if let Some(stripped) = formula_outer_group_for_offset(term)? {
        for (inner_sign, expanded_term) in split_formula_signed_terms(stripped)? {
            extract_formula_offset_term(
                sign * inner_sign,
                &expanded_term,
                factors,
                numeric_covariates,
                offsets,
                saw_offset,
                remaining_terms,
            )?;
        }
        return Ok(());
    }

    let Some(offset_values) = formula_offset_values(term, factors, numeric_covariates)? else {
        if sign < 0 {
            remaining_terms.push(format!("- {term}"));
        } else {
            remaining_terms.push(term.to_string());
        }
        return Ok(());
    };
    if sign < 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula offset terms cannot be subtracted".to_string(),
        });
    }
    *saw_offset = true;
    if offsets.is_empty() {
        offsets.resize(offset_values.len(), 0.0);
    }
    if offset_values.len() != offsets.len() {
        return Err(invalid_dimensions(
            "formula offset values",
            offsets.len(),
            offset_values.len(),
        ));
    }
    for (idx, value) in offset_values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offset '{term}' is non-finite at sample {idx}"),
            });
        }
        offsets[idx] += value;
        if !offsets[idx].is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offsets sum to a non-finite value at sample {idx}"),
            });
        }
    }
    Ok(())
}

fn formula_outer_group_for_offset(term: &str) -> Result<Option<&str>, DeseqError> {
    if formula_offset_expression(term)?.is_some() {
        return Ok(None);
    }
    if !term.contains("offset") {
        return Ok(None);
    }
    let stripped = strip_formula_outer_parentheses(term)?;
    Ok((stripped != term).then_some(stripped))
}

fn formula_offset_values(
    term: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<Option<Vec<f64>>, DeseqError> {
    let Some(expression) = formula_offset_expression(term)? else {
        return Ok(None);
    };
    if let Ok(offset_name) = formula_variable_name(expression) {
        let covariate = numeric_covariates
            .iter()
            .find(|candidate| candidate.name == offset_name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula offset numeric covariate '{offset_name}' is not present in supplied design metadata"
                ),
            })?;
        return Ok(Some(covariate.values.to_vec()));
    }
    for transform in FORMULA_NUMERIC_TRANSFORMS {
        let Some(inner) = formula_numeric_transform_inner_expression(expression, transform) else {
            continue;
        };
        let (_replacement, transformed_covariates) =
            formula_numeric_transform_term(transform, inner.trim(), factors, numeric_covariates)?;
        let [(_name, values)] = transformed_covariates.as_slice() else {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offset term '{term}' must produce one numeric vector"),
            });
        };
        return Ok(Some(values.clone()));
    }
    Err(DeseqError::InvalidOptions {
        reason: format!(
            "formula offset expression '{expression}' is not a supported numeric covariate or transform"
        ),
    })
}

fn formula_offset_expression(term: &str) -> Result<Option<&str>, DeseqError> {
    let term = term.trim();
    let Some(after_open) = formula_call_after_open(term, 0, "offset") else {
        if term == "offset" {
            return Err(DeseqError::InvalidOptions {
                reason: "formula offset term 'offset' must be offset(numeric)".to_string(),
            });
        }
        return Ok(None);
    };
    let Some(inner) = term[after_open..].strip_suffix(')') else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula offset term '{term}' must be offset(numeric)"),
        });
    };
    let expression = inner.trim();
    if expression.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula offset term '{term}' must be offset(numeric)"),
        });
    }
    Ok(Some(expression))
}

fn expand_formula_numeric_transform_terms(
    rhs: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<FormulaDerivedNumericCovariate>), DeseqError> {
    let mut expanded = String::with_capacity(rhs.len());
    let mut derived = Vec::new();
    let mut remainder = rhs;
    while let Some((start, after_open_idx, transform)) = next_formula_numeric_transform(remainder) {
        expanded.push_str(&remainder[..start]);
        let after_open = &remainder[after_open_idx..];
        let close = formula_transform_close_index(after_open)?;
        let expression = after_open[..close].trim();
        let (replacement, transformed_covariates) =
            formula_numeric_transform_term(transform, expression, factors, numeric_covariates)?;
        for (name, values) in transformed_covariates {
            if !derived
                .iter()
                .any(|candidate: &FormulaDerivedNumericCovariate| candidate.name == name)
            {
                derived.push(FormulaDerivedNumericCovariate { name, values });
            }
        }
        expanded.push_str(&replacement);
        remainder = &after_open[close + 1..];
    }
    expanded.push_str(remainder);
    Ok((expanded, derived))
}

fn formula_transform_close_index(after_open: &str) -> Result<usize, DeseqError> {
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for (idx, character) in after_open.char_indices() {
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        match character {
            '(' => depth += 1,
            ')' if depth == 0 => return Ok(idx),
            ')' => depth -= 1,
            _ => {}
        }
    }
    if in_backticks {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform has unbalanced backtick-quoted variable name".to_string(),
        });
    }
    if quote.is_some() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform has unbalanced quoted string".to_string(),
        });
    }
    Err(DeseqError::InvalidOptions {
        reason: "formula transform has unbalanced parentheses".to_string(),
    })
}

#[derive(Clone, Copy, Debug)]
struct FormulaNumericTransform {
    name: &'static str,
    label: &'static str,
    apply: fn(f64) -> f64,
}

const FORMULA_NUMERIC_TRANSFORMS: [FormulaNumericTransform; 11] = [
    FormulaNumericTransform {
        name: "poly",
        label: "poly",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        name: "scale",
        label: "scale",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        name: "log10",
        label: "log10",
        apply: f64::log10,
    },
    FormulaNumericTransform {
        name: "log1p",
        label: "log1p",
        apply: f64::ln_1p,
    },
    FormulaNumericTransform {
        name: "log2",
        label: "log2",
        apply: f64::log2,
    },
    FormulaNumericTransform {
        name: "sqrt",
        label: "sqrt",
        apply: f64::sqrt,
    },
    FormulaNumericTransform {
        name: "as.numeric",
        label: "as.numeric",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        name: "as.double",
        label: "as.double",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        name: "as.integer",
        label: "as.integer",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        name: "log",
        label: "log",
        apply: f64::ln,
    },
    FormulaNumericTransform {
        name: "I",
        label: "I",
        apply: std::convert::identity,
    },
];

fn next_formula_numeric_transform(rhs: &str) -> Option<(usize, usize, FormulaNumericTransform)> {
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    let mut idx = 0_usize;
    while idx < rhs.len() {
        let character = rhs[idx..]
            .chars()
            .next()
            .expect("idx is inside rhs char boundary");
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            idx += character.len_utf8();
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            idx += character.len_utf8();
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            idx += character.len_utf8();
            continue;
        }
        if !in_backticks {
            for transform in FORMULA_NUMERIC_TRANSFORMS {
                if let Some(after_open_idx) = formula_call_after_open(rhs, idx, transform.name) {
                    return Some((idx, after_open_idx, transform));
                }
            }
        }
        idx += character.len_utf8();
    }
    None
}

fn formula_numeric_transform_inner_expression(
    expression: &str,
    transform: FormulaNumericTransform,
) -> Option<&str> {
    let trimmed = expression.trim();
    let after_open = formula_call_after_open(trimmed, 0, transform.name)?;
    trimmed[after_open..].strip_suffix(')')
}

fn formula_call_after_open(text: &str, idx: usize, name: &str) -> Option<usize> {
    if !formula_call_has_left_boundary(text, idx) || !text[idx..].starts_with(name) {
        return None;
    }
    let mut after_name = idx + name.len();
    while after_name < text.len() {
        let character = text[after_name..]
            .chars()
            .next()
            .expect("after_name is inside text char boundary");
        if !character.is_whitespace() {
            break;
        }
        after_name += character.len_utf8();
    }
    if text[after_name..].starts_with('(') {
        Some(after_name + 1)
    } else {
        None
    }
}

fn formula_call_has_left_boundary(text: &str, idx: usize) -> bool {
    if idx == 0 {
        return true;
    }
    let Some(previous) = text[..idx].chars().next_back() else {
        return true;
    };
    !(previous.is_alphanumeric() || previous == '_' || previous == '.')
}

fn formula_numeric_transform_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    match transform.label {
        "I" => {
            let (name, values) =
                formula_numeric_identity_or_power_term(expression, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
        "poly" => formula_numeric_poly_term(expression, numeric_covariates),
        "scale" => {
            let (name, values) = formula_numeric_scale_term(expression, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
        _ => {
            let (name, values) =
                formula_numeric_function_term(transform, expression, factors, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
    }
}

fn formula_numeric_function_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let (numeric_text, base) = formula_numeric_function_arguments(transform, expression)?;
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    if matches!(transform.label, "as.numeric" | "as.double" | "as.integer") {
        if let Some(covariate) = maybe_formula_numeric_covariate(Some(numeric_name), numeric_covariates)? {
            return formula_numeric_covariate_function_values(transform, expression, covariate, base);
        }
        return formula_factor_numeric_coercion_term(transform, expression, numeric_name, factors);
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    formula_numeric_covariate_function_values(transform, expression, covariate, base)
}

fn formula_numeric_covariate_function_values(
    transform: FormulaNumericTransform,
    expression: &str,
    covariate: &ExpandedNumericSpec<'_>,
    base: Option<f64>,
) -> Result<(String, Vec<f64>), DeseqError> {
    let numeric_name = covariate.name;
    let mut values = Vec::with_capacity(covariate.values.len());
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        let transformed = if let Some(base) = base {
            value.log(base)
        } else {
            (transform.apply)(value)
        };
        if !transformed.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform '{}({expression})' produced non-finite value at sample {idx}",
                    transform.label
                ),
            });
        }
        values.push(transformed);
    }
    let name = if let Some(base) = base {
        format!("{numeric_name}_log_base_{}", formula_numeric_label(base))
    } else {
        format!("{}_{}", numeric_name, transform.label)
    };
    Ok((name, values))
}

fn formula_factor_numeric_coercion_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factor_name: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let factor = find_formula_factor(factor_name, factors)?;
    let levels = formula_factor_code_levels(factor);
    let values = factor
        .sample_levels
        .iter()
        .enumerate()
        .map(|(sample, level)| {
            levels
                .iter()
                .position(|candidate| candidate == level)
                .map(|idx| idx as f64 + 1.0)
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform '{}({expression})' has sample level '{}' outside factor levels at sample {sample}",
                        transform.label, level
                    ),
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((format!("{}_{}", factor.factor, transform.label), values))
}

fn formula_factor_code_levels(factor: &ExpandedFactorSpec<'_>) -> Vec<String> {
    if let Some(levels) = factor.levels {
        return levels.to_vec();
    }
    let mut levels = Vec::new();
    for level in factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.clone());
        }
    }
    levels
}

fn formula_numeric_function_arguments(
    transform: FormulaNumericTransform,
    expression: &str,
) -> Result<(String, Option<f64>), DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    let mut numeric_text = None;
    let mut base = None;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        numeric_text.is_some(),
                        transform.label,
                        "x",
                        expression,
                    )?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "base" if transform.label == "log" => {
                    reject_duplicate_formula_transform_argument(
                        base.is_some(),
                        transform.label,
                        "base",
                        expression,
                    )?;
                    base = Some(parse_formula_log_base(value, expression)?);
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform '{}({expression})' has unsupported argument '{argument}'",
                            transform.label
                        ),
                    });
                }
            }
        }
        match positional_idx {
            0 if numeric_text.is_none() => {
                numeric_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if transform.label == "log" && base.is_none() => {
                base = Some(parse_formula_log_base(trimmed, expression)?);
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' has unsupported argument '{argument}'",
                transform.label
            ),
        });
    }
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' must provide a numeric covariate",
                transform.label
            ),
        });
    };
    Ok((numeric_text, base))
}

fn parse_formula_log_base(value: &str, expression: &str) -> Result<f64, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    let parsed = value.parse::<f64>().map_err(|_| DeseqError::InvalidOptions {
        reason: format!("formula transform 'log({expression})' base must be finite and positive"),
    })?;
    if !parsed.is_finite() || parsed <= 0.0 || parsed == 1.0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'log({expression})' base must be finite and positive"),
        });
    }
    Ok(parsed)
}

fn formula_numeric_label(value: f64) -> String {
    value.to_string().replace('.', "_")
}

fn find_formula_numeric_covariate<'a>(
    numeric_name: &str,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<&'a ExpandedNumericSpec<'a>, DeseqError> {
    let exact = numeric_covariates
        .iter()
        .filter(|candidate| candidate.name == numeric_name)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [covariate] => return Ok(*covariate),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula numeric covariate '{numeric_name}' appears more than once"),
            });
        }
    }

    let aliases = numeric_covariates
        .iter()
        .filter(|candidate| {
            r_like_name_candidates(candidate.name)
                .into_iter()
                .any(|candidate| candidate == numeric_name)
        })
        .collect::<Vec<_>>();
    match aliases.as_slice() {
        [covariate] => Ok(*covariate),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' is not present in supplied design metadata"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}

fn maybe_formula_numeric_covariate<'a>(
    numeric_name: Option<&str>,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Option<&'a ExpandedNumericSpec<'a>>, DeseqError> {
    let Some(numeric_name) = numeric_name else {
        return Ok(None);
    };
    let matches_any_covariate = numeric_covariates
        .iter()
        .any(|candidate| candidate.name == numeric_name)
        || numeric_covariates.iter().any(|candidate| {
            r_like_name_candidates(candidate.name)
                .into_iter()
                .any(|candidate| candidate == numeric_name)
        });
    if !matches_any_covariate {
        return Ok(None);
    }
    find_formula_numeric_covariate(numeric_name, numeric_covariates).map(Some)
}

fn formula_numeric_scale_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    let mut numeric_text = None;
    let mut center = FormulaScaleOption::Auto;
    let mut scale = FormulaScaleOption::Auto;
    let mut center_seen = false;
    let mut scale_seen = false;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if let Some(value) = normalized.strip_prefix("center=") {
            reject_duplicate_formula_transform_argument(
                center_seen,
                "scale",
                "center",
                expression,
            )?;
            center_seen = true;
            center = parse_formula_scale_option(value, "center", expression)?;
            continue;
        }
        if let Some(value) = normalized.strip_prefix("scale=") {
            reject_duplicate_formula_transform_argument(
                scale_seen,
                "scale",
                "scale",
                expression,
            )?;
            scale_seen = true;
            scale = parse_formula_scale_option(value, "scale", expression)?;
            continue;
        }
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        numeric_text.is_some(),
                        "scale",
                        "x",
                        expression,
                    )?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "center" => {
                    reject_duplicate_formula_transform_argument(
                        center_seen,
                        "scale",
                        "center",
                        expression,
                    )?;
                    center_seen = true;
                    center = parse_formula_scale_option(value, "center", expression)?;
                    continue;
                }
                "scale" => {
                    reject_duplicate_formula_transform_argument(
                        scale_seen,
                        "scale",
                        "scale",
                        expression,
                    )?;
                    scale_seen = true;
                    scale = parse_formula_scale_option(value, "scale", expression)?;
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform 'scale({expression})' has unsupported argument '{argument}'"
                        ),
                    });
                }
            }
        }
        match positional_idx {
            0 if numeric_text.is_none() => {
                numeric_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if !center_seen => {
                center_seen = true;
                center = parse_formula_scale_option(&normalized, "center", expression)?;
                positional_idx += 1;
                continue;
            }
            2 if !scale_seen => {
                scale_seen = true;
                scale = parse_formula_scale_option(&normalized, "scale", expression)?;
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' must provide a numeric covariate"
            ),
        });
    };
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    if scale == FormulaScaleOption::Auto && covariate.values.len() < 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' requires at least two samples"
            ),
        });
    }
    let mut sum = 0.0_f64;
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula transform 'scale({expression})' received non-finite value at sample {idx}"),
            });
        }
        sum += value;
        if !sum.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'scale({expression})' produced non-finite center"
                ),
            });
        }
    }
    let center_value = match center {
        FormulaScaleOption::Auto => sum / covariate.values.len() as f64,
        FormulaScaleOption::Disabled => 0.0,
        FormulaScaleOption::Explicit(value) => value,
    };
    let divisor = match scale {
        FormulaScaleOption::Auto => {
            formula_scale_auto_divisor(expression, covariate.values, center_value)?
        }
        FormulaScaleOption::Disabled => 1.0,
        FormulaScaleOption::Explicit(value) => value,
    };
    let values = covariate
        .values
        .iter()
        .map(|value| (value - center_value) / divisor)
        .collect::<Vec<_>>();
    let label = match (center, scale) {
        (FormulaScaleOption::Auto, FormulaScaleOption::Auto) => "scale",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Auto) => "scale_uncentered",
        (FormulaScaleOption::Auto, FormulaScaleOption::Disabled) => "center",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Disabled) => "identity",
        (FormulaScaleOption::Explicit(_), FormulaScaleOption::Disabled) => "centered",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Explicit(_)) => "scaled",
        _ => "centered_scaled",
    };
    Ok((format!("{numeric_name}_{label}"), values))
}

fn formula_scale_auto_divisor(
    expression: &str,
    values: &[f64],
    center: f64,
) -> Result<f64, DeseqError> {
    let mut sum_squares = 0.0_f64;
    for value in values.iter().copied() {
        let adjusted = value - center;
        let square = adjusted * adjusted;
        if !square.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'scale({expression})' produced non-finite scale"
                ),
            });
        }
        sum_squares += square;
        if !sum_squares.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'scale({expression})' produced non-finite scale"
                ),
            });
        }
    }
    let divisor = (sum_squares / (values.len() as f64 - 1.0)).sqrt();
    if !divisor.is_finite() || divisor <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'scale({expression})' requires non-zero scale"),
        });
    }
    Ok(divisor)
}

fn parse_formula_scale_option(
    value: &str,
    argument: &str,
    expression: &str,
) -> Result<FormulaScaleOption, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    match value.to_ascii_lowercase().as_str() {
        "true" | "t" => Ok(FormulaScaleOption::Auto),
        "false" | "f" => Ok(FormulaScaleOption::Disabled),
        _ => {
            let parsed = value
                .parse::<f64>()
                .map_err(|_| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'scale({expression})' argument '{argument}' must be TRUE, FALSE, or a finite number"
                    ),
                })?;
            if !parsed.is_finite() || (argument == "scale" && parsed <= 0.0) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'scale({expression})' argument '{argument}' must be TRUE, FALSE, or a finite positive scale"
                    ),
                });
            }
            Ok(FormulaScaleOption::Explicit(parsed))
        }
    }
}

fn formula_numeric_poly_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    if arguments.len() < 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' must provide numeric and degree"
            ),
        });
    }
    let mut numeric_text = None;
    let mut degree_text = None;
    let mut raw = false;
    let mut raw_seen = false;
    let mut simple_seen = false;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if normalized == "raw=true" || normalized == "raw=t" {
            reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
            raw_seen = true;
            raw = true;
            continue;
        }
        if normalized == "raw=false" || normalized == "raw=f" {
            reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
            raw_seen = true;
            raw = false;
            continue;
        }
        if normalized == "simple=true"
            || normalized == "simple=t"
            || normalized == "simple=false"
            || normalized == "simple=f"
        {
            reject_duplicate_formula_poly_argument(simple_seen, "simple", expression)?;
            simple_seen = true;
            continue;
        }
        if let Some(named_degree) = normalized.strip_prefix("degree=") {
            reject_duplicate_formula_poly_argument(degree_text.is_some(), "degree", expression)?;
            degree_text = Some(named_degree.to_string());
            continue;
        }
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_poly_argument(numeric_text.is_some(), "x", expression)?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "degree" => {
                    reject_duplicate_formula_poly_argument(
                        degree_text.is_some(),
                        "degree",
                        expression,
                    )?;
                    degree_text = Some(value.to_string());
                    continue;
                }
                "raw" => {
                    reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
                    raw_seen = true;
                    raw = parse_formula_bool_option(value, "raw", expression)?;
                    continue;
                }
                "simple" => {
                    reject_duplicate_formula_poly_argument(simple_seen, "simple", expression)?;
                    simple_seen = true;
                    parse_formula_bool_option(value, "simple", expression)?;
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform 'poly({expression})' has unsupported argument '{argument}'"
                        ),
                    });
                }
            }
        }
        match positional_idx {
            0 if numeric_text.is_none() => {
                numeric_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if degree_text.is_none() => {
                degree_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        if positional_idx < 2 {
            positional_idx += 1;
            continue;
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let Some(degree_text) = degree_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'poly({expression})' must provide a degree"),
        });
    };
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' must provide a numeric covariate"
            ),
        });
    };
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    let degree = degree_text
        .parse::<i32>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula polynomial degree '{degree_text}' must be an integer"),
        })?;
    if !(1..=16).contains(&degree) {
        return Err(DeseqError::InvalidOptions {
            reason: "formula polynomial degrees must be integers from 1 through 16".to_string(),
        });
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    if raw {
        return formula_numeric_raw_poly_columns(expression, numeric_name, covariate.values, degree);
    }
    formula_numeric_orthogonal_poly_columns(expression, numeric_name, covariate.values, degree)
}

fn reject_duplicate_formula_poly_argument(
    seen: bool,
    argument: &str,
    expression: &str,
) -> Result<(), DeseqError> {
    reject_duplicate_formula_transform_argument(seen, "poly", argument, expression)
}

fn reject_duplicate_formula_transform_argument(
    seen: bool,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<(), DeseqError> {
    if seen {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{transform}({expression})' argument '{argument}' matched multiple times"
            ),
        });
    }
    Ok(())
}

fn split_formula_named_argument(argument: &str) -> Option<(String, &str)> {
    let equals = formula_top_level_equals_index(argument)?;
    let (key, value_with_equals) = argument.split_at(equals);
    let value = &value_with_equals[1..];
    let key = key
        .split_whitespace()
        .collect::<String>()
        .to_ascii_lowercase();
    Some((key, value.trim()))
}

fn formula_top_level_equals_index(argument: &str) -> Option<usize> {
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for (idx, character) in argument.char_indices() {
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        match character {
            '(' => depth += 1,
            ')' => depth -= 1,
            '=' if depth == 0 => return Some(idx),
            _ => {}
        }
    }
    None
}

fn parse_formula_bool_option(
    value: &str,
    argument: &str,
    expression: &str,
) -> Result<bool, DeseqError> {
    match value
        .split_whitespace()
        .collect::<String>()
        .to_ascii_lowercase()
        .as_str()
    {
        "true" | "t" => Ok(true),
        "false" | "f" => Ok(false),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' argument '{argument}' must be TRUE or FALSE"
            ),
        }),
    }
}

fn parse_formula_string_argument(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<String, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?.trim();
    let parsed = if value.starts_with('"') || value.starts_with('\'') {
        parse_formula_quoted_string(value, transform, argument, expression)?
    } else {
        value.to_string()
    };
    if parsed.is_empty() || parsed.contains('`') {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{transform}({expression})' argument '{argument}' must be a simple level string"
            ),
        });
    }
    Ok(parsed)
}

fn parse_formula_quoted_string(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<String, DeseqError> {
    let quote = value
        .chars()
        .next()
        .expect("quoted formula string is not empty");
    let Some(inner) = value
        .strip_prefix(quote)
        .and_then(|candidate| candidate.strip_suffix(quote))
    else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{transform}({expression})' argument '{argument}' must be a balanced quoted string"
            ),
        });
    };
    let mut parsed = String::with_capacity(inner.len());
    let mut chars = inner.chars();
    while let Some(character) = chars.next() {
        if character == '\\' {
            let Some(escaped) = chars.next() else {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform '{transform}({expression})' argument '{argument}' has a trailing string escape"
                    ),
                });
            };
            parsed.push(match escaped {
                'n' => '\n',
                'r' => '\r',
                't' => '\t',
                other => other,
            });
        } else {
            parsed.push(character);
        }
    }
    Ok(parsed)
}

fn parse_formula_string_vector_argument(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<Vec<String>, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?.trim();
    let Some(inner) = value
        .strip_prefix("c(")
        .and_then(|candidate| candidate.strip_suffix(')'))
    else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{transform}({expression})' argument '{argument}' must be c(\"level\", ...)"
            ),
        });
    };
    let mut levels = Vec::new();
    for level in split_formula_transform_arguments(inner)? {
        let parsed = parse_formula_string_argument(&level, transform, argument, expression)?;
        if levels.iter().any(|candidate| candidate == &parsed) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform '{transform}({expression})' argument '{argument}' contains duplicate level '{parsed}'"
                ),
            });
        }
        levels.push(parsed);
    }
    if levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{transform}({expression})' argument '{argument}' must contain at least one level"
            ),
        });
    }
    Ok(levels)
}

fn formula_string_vector_label(values: &[String]) -> String {
    let quoted = values
        .iter()
        .map(|value| format!("\"{}\"", formula_quote_string(value)))
        .collect::<Vec<_>>();
    format!("c({})", quoted.join(", "))
}

fn formula_quote_string(value: &str) -> String {
    let mut quoted = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '\\' => quoted.push_str("\\\\"),
            '"' => quoted.push_str("\\\""),
            '\n' => quoted.push_str("\\n"),
            '\r' => quoted.push_str("\\r"),
            '\t' => quoted.push_str("\\t"),
            other => quoted.push(other),
        }
    }
    quoted
}

fn formula_numeric_raw_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let mut replacement_terms = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for exponent in 1..=degree {
        let name = format!("{numeric_name}_poly_{exponent}");
        let mut column_values = Vec::with_capacity(values.len());
        for (idx, value) in values.iter().copied().enumerate() {
            let transformed = value.powi(exponent);
            if !transformed.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'poly({expression})' produced non-finite value at sample {idx}"
                    ),
                });
            }
            column_values.push(transformed);
        }
        replacement_terms.push(formula_variable_term(&name));
        derived.push((name, column_values));
    }
    Ok((format!("({})", replacement_terms.join(" + ")), derived))
}

fn formula_numeric_orthogonal_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    validate_formula_poly_values(expression, values, degree)?;
    let mut basis = Vec::with_capacity(degree as usize + 1);
    let constant_norm = (values.len() as f64).sqrt();
    basis.push(vec![1.0 / constant_norm; values.len()]);
    let mut replacement_terms = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for exponent in 1..=degree {
        let mut column_values = Vec::with_capacity(values.len());
        for (idx, value) in values.iter().copied().enumerate() {
            let transformed = value.powi(exponent);
            if !transformed.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'poly({expression})' produced non-finite value at sample {idx}"
                    ),
                });
            }
            column_values.push(transformed);
        }
        for previous in &basis {
            let projection = formula_dot_product(&column_values, previous);
            for (value, previous_value) in column_values.iter_mut().zip(previous) {
                *value -= projection * previous_value;
            }
        }
        let norm = formula_dot_product(&column_values, &column_values).sqrt();
        if !norm.is_finite() || norm <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'poly({expression})' produced a dependent polynomial column"
                ),
            });
        }
        for value in &mut column_values {
            *value /= norm;
        }
        let name = format!("poly({numeric_name}, {degree}){exponent}");
        replacement_terms.push(formula_variable_term(&name));
        derived.push((name, column_values.clone()));
        basis.push(column_values);
    }
    Ok((format!("({})", replacement_terms.join(" + ")), derived))
}

fn validate_formula_poly_values(
    expression: &str,
    values: &[f64],
    degree: i32,
) -> Result<(), DeseqError> {
    let mut unique = Vec::new();
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'poly({expression})' received non-finite value at sample {idx}"
                ),
            });
        }
        if !unique.contains(&value) {
            unique.push(value);
        }
    }
    if unique.len() <= degree as usize {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' degree must be less than the number of unique values"
            ),
        });
    }
    Ok(())
}

fn formula_dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}

fn split_formula_transform_arguments(expression: &str) -> Result<Vec<String>, DeseqError> {
    split_formula_top_level(expression, ',')
}

fn formula_numeric_identity_or_power_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    if let Some((key, value)) = split_formula_named_argument(expression.trim()) {
        if key == "x" {
            return formula_numeric_identity_or_power_term(value, numeric_covariates);
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' has unsupported argument '{key}'"),
        });
    }
    let Some((numeric_name, exponent_text)) = expression.split_once('^') else {
        if let Some((left, op, right)) = split_formula_numeric_scalar_expression(expression) {
            return formula_numeric_scalar_expression_term(
                expression,
                left,
                op,
                right,
                numeric_covariates,
            );
        }
        return formula_numeric_identity_term(expression, numeric_covariates);
    };
    formula_numeric_power_term(expression, numeric_name, exponent_text, numeric_covariates)
}

fn split_formula_numeric_scalar_expression(expression: &str) -> Option<(&str, char, &str)> {
    let mut in_backticks = false;
    for (idx, op) in expression.char_indices() {
        if op == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        if idx == 0 || !matches!(op, '+' | '-' | '*' | '/') {
            continue;
        }
        let left = expression[..idx].trim();
        let right = expression[idx + op.len_utf8()..].trim();
        if left.is_empty() || right.is_empty() {
            continue;
        }
        return Some((left, op, right));
    }
    None
}

fn formula_numeric_scalar_expression_term(
    expression: &str,
    left: &str,
    op: char,
    right: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let left_name = formula_variable_name(left).ok();
    let right_name = formula_variable_name(right).ok();
    let left_covariate = maybe_formula_numeric_covariate(left_name, numeric_covariates)?;
    let right_covariate = maybe_formula_numeric_covariate(right_name, numeric_covariates)?;
    let (numeric_name, scalar_text, scalar_on_left, covariate) =
        match (left_covariate, right_covariate) {
            (Some(covariate), None) => (covariate.name, right, false, covariate),
            (None, Some(covariate)) => (covariate.name, left, true, covariate),
            (Some(left_covariate), Some(right_covariate)) => {
                return formula_numeric_binary_expression_term(
                    expression,
                    left_covariate.name,
                    op,
                    right_covariate.name,
                    left_covariate,
                    right_covariate,
                );
            }
            (None, None) => {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'I({expression})' must contain one supplied numeric covariate"
                    ),
                });
            }
        };
    let scalar = parse_formula_scalar(scalar_text, expression)?;
    if op == '/' && !scalar_on_left && scalar == 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'I({expression})' scalar '{scalar_text}' must be non-zero for division"
            ),
        });
    }
    let mut values = Vec::with_capacity(covariate.values.len());
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' received non-finite value at sample {idx}"
                ),
            });
        }
        let transformed = match (scalar_on_left, op) {
            (_, '+') => value + scalar,
            (false, '-') => value - scalar,
            (true, '-') => scalar - value,
            (_, '*') => value * scalar,
            (false, '/') => value / scalar,
            (true, '/') => scalar / value,
            _ => unreachable!("validated formula scalar operator"),
        };
        if !transformed.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' produced non-finite value at sample {idx}"
                ),
            });
        }
        values.push(transformed);
    }
    Ok((
        format!(
            "{}_{}_{}",
            numeric_name,
            formula_scalar_operator_label(op, scalar_on_left),
            formula_scalar_value_label(scalar_text)
        ),
        values,
    ))
}

fn formula_numeric_binary_expression_term(
    expression: &str,
    left_name: &str,
    op: char,
    right_name: &str,
    left: &ExpandedNumericSpec<'_>,
    right: &ExpandedNumericSpec<'_>,
) -> Result<(String, Vec<f64>), DeseqError> {
    if left.values.len() != right.values.len() {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'I({expression})' received numeric covariates with different sample counts"
            ),
        });
    }
    let mut values = Vec::with_capacity(left.values.len());
    for (idx, (&left_value, &right_value)) in left.values.iter().zip(right.values.iter()).enumerate()
    {
        if !left_value.is_finite() || !right_value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' received non-finite value at sample {idx}"
                ),
            });
        }
        let transformed = match op {
            '+' => left_value + right_value,
            '-' => left_value - right_value,
            '*' => left_value * right_value,
            '/' => left_value / right_value,
            _ => unreachable!("validated formula binary operator"),
        };
        if !transformed.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' produced non-finite value at sample {idx}"
                ),
            });
        }
        values.push(transformed);
    }
    Ok((
        format!(
            "{}_{}_{}",
            left_name,
            formula_binary_operator_label(op),
            right_name
        ),
        values,
    ))
}

fn parse_formula_scalar(value: &str, expression: &str) -> Result<f64, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    let scalar = value
        .parse::<f64>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' scalar '{value}' must be finite"),
        })?;
    if !scalar.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' scalar '{value}' must be finite"),
        });
    }
    Ok(scalar)
}

fn formula_binary_operator_label(op: char) -> &'static str {
    match op {
        '+' => "plus",
        '-' => "minus",
        '*' => "times",
        '/' => "div",
        _ => unreachable!("validated formula binary operator"),
    }
}

fn formula_scalar_operator_label(op: char, scalar_on_left: bool) -> &'static str {
    match (scalar_on_left, op) {
        (_, '+') => "plus",
        (false, '-') => "minus",
        (true, '-') => "rminus",
        (_, '*') => "times",
        (false, '/') => "div",
        (true, '/') => "rdiv",
        _ => unreachable!("validated formula scalar operator"),
    }
}

fn formula_scalar_value_label(value: &str) -> String {
    let mut label = String::new();
    for ch in value.trim().chars() {
        match ch {
            '-' => label.push_str("neg"),
            '+' => label.push_str("pos"),
            '.' => label.push('p'),
            ch if ch.is_ascii_alphanumeric() => label.push(ch),
            _ => label.push('_'),
        }
    }
    label
}

fn formula_numeric_identity_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let (numeric_term, sign, label) = match expression.trim().as_bytes().first().copied() {
        Some(b'+') => (expression.trim()[1..].trim(), 1.0, "identity"),
        Some(b'-') => (expression.trim()[1..].trim(), -1.0, "neg"),
        _ => (expression.trim(), 1.0, "identity"),
    };
    let numeric_name = formula_variable_name(numeric_term)?;
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    let mut values = Vec::with_capacity(covariate.values.len());
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' received non-finite value at sample {idx}"
                ),
            });
        }
        values.push(sign * value);
    }
    Ok((format!("{numeric_name}_{label}"), values))
}

fn formula_numeric_power_term(
    expression: &str,
    numeric_name: &str,
    exponent_text: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let numeric_name = formula_variable_name(numeric_name.trim())?;
    let exponent = exponent_text
        .trim()
        .parse::<i32>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula transform exponent '{exponent_text}' must be an integer"),
        })?;
    if !(2..=16).contains(&exponent) {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform powers must be integers from 2 through 16".to_string(),
        });
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    let mut values = Vec::with_capacity(covariate.values.len());
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        let transformed = value.powi(exponent);
        if !transformed.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' produced non-finite value at sample {idx}"
                ),
            });
        }
        values.push(transformed);
    }
    Ok((format!("{numeric_name}_pow_{exponent}"), values))
}

fn formula_variable_term(name: &str) -> String {
    if validate_formula_variable(name).is_ok() {
        name.to_string()
    } else {
        format!("`{name}`")
    }
}
