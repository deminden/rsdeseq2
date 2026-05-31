/// Build an expanded additive design from a primitive DESeq2-style formula subset.
///
/// Supported right-hand-side terms are intercept-preserving main effects
/// (`condition`, `dose`), intercept-only `1`, pairwise interactions
/// (`condition:dose`), nested shorthand (`condition/batch`), and `*`
/// shorthand for main effects plus interactions (`condition*dose`).
/// Interaction variables can appear without corresponding main effects. The
/// reported standard-design interaction columns then follow R model-matrix
/// treatment-coding shape for the supported primitive terms. Intercept removal
/// with `0` or `-1` is supported for these primitive terms. Formula
/// interactions can contain two or more variables. Primitive `- term`
/// subtraction is supported for the same term subset. Integer numeric power
/// transforms are materialized as derived numeric covariates. Supported
/// `offset(numeric)` terms are parsed by [`expanded_formula_design_with_offsets`];
/// this compatibility helper returns only the design surface.
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
/// supplied numeric covariate. Multiple offset terms are summed sample-wise.
pub fn expanded_formula_design_with_offsets<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    let rhs = formula_rhs(formula)?;
    let (rhs, offsets) = extract_formula_offsets(rhs, numeric_covariates)?;
    let (rhs, derived_numeric_covariates) =
        expand_formula_numeric_transform_terms(&rhs, numeric_covariates)?;
    let mut all_numeric_covariates =
        Vec::with_capacity(numeric_covariates.len() + derived_numeric_covariates.len());
    all_numeric_covariates.extend(numeric_covariates.iter().cloned());
    all_numeric_covariates.extend(derived_numeric_covariates.iter().map(|covariate| {
        ExpandedNumericSpec {
            name: covariate.name.as_str(),
            values: covariate.values.as_slice(),
        }
    }));
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
            remove_formula_term(term, factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if term == "1" {
            continue;
        }
        if term == "0" || term == "-1" {
            state.has_intercept = false;
            continue;
        }
        if term.contains('^') || term.contains('(') || term.contains(')') {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula term '{term}' is not supported by the primitive parser"),
            });
        }
        if term.contains('/') {
            add_nested_formula_term(term, factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if term.contains('*') {
            add_star_formula_term(term, factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        if term.contains(':') {
            add_interaction_formula_term(term, factors, &all_numeric_covariates, &mut state)?;
            continue;
        }
        add_main_formula_term(
            term,
            factors,
            &all_numeric_covariates,
            &mut state.selected_factors,
            &mut state.selected_numeric_covariates,
        )?;
    }

    let design = expanded_formula_design_from_state(&state, factors, &all_numeric_covariates)?;
    Ok(ExpandedFormulaDesignWithOffsets { design, offsets })
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
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let signed_terms = split_formula_signed_terms(rhs)?;
    let n_samples = numeric_covariates
        .first()
        .map(|covariate| covariate.values.len())
        .unwrap_or(0);
    let mut offsets = vec![0.0; n_samples];
    let mut remaining_terms = Vec::new();
    for (sign, term) in signed_terms {
        let Some(offset_name) = formula_offset_name(&term)? else {
            if sign < 0 {
                remaining_terms.push(format!("- {term}"));
            } else {
                remaining_terms.push(term);
            }
            continue;
        };
        if sign < 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "formula offset terms cannot be subtracted".to_string(),
            });
        }
        let covariate = numeric_covariates
            .iter()
            .find(|candidate| candidate.name == offset_name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula offset numeric covariate '{offset_name}' is not present in supplied design metadata"
                ),
            })?;
        if offsets.is_empty() {
            offsets.resize(covariate.values.len(), 0.0);
        }
        if covariate.values.len() != offsets.len() {
            return Err(invalid_dimensions(
                "formula offset values",
                offsets.len(),
                covariate.values.len(),
            ));
        }
        for (idx, value) in covariate.values.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula offset '{offset_name}' is non-finite at sample {idx}"),
                });
            }
            offsets[idx] += value;
            if !offsets[idx].is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula offsets sum to a non-finite value at sample {idx}"),
                });
            }
        }
    }
    if remaining_terms.is_empty() {
        remaining_terms.push("1".to_string());
    }
    Ok((join_formula_terms(&remaining_terms), offsets))
}

fn formula_offset_name(term: &str) -> Result<Option<&str>, DeseqError> {
    let term = term.trim();
    if !term.starts_with("offset") {
        return Ok(None);
    }
    let Some(inner) = term
        .strip_prefix("offset(")
        .and_then(|value| value.strip_suffix(')'))
    else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula offset term '{term}' must be offset(numeric)"),
        });
    };
    let inner = inner.trim();
    validate_formula_variable(inner)?;
    Ok(Some(inner))
}

fn expand_formula_numeric_transform_terms(
    rhs: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<FormulaDerivedNumericCovariate>), DeseqError> {
    let mut expanded = String::with_capacity(rhs.len());
    let mut derived = Vec::new();
    let mut remainder = rhs;
    while let Some((start, transform)) = next_formula_numeric_transform(remainder) {
        expanded.push_str(&remainder[..start]);
        let after_open = &remainder[start + transform.prefix.len()..];
        let Some(close) = after_open.find(')') else {
            return Err(DeseqError::InvalidOptions {
                reason: "formula transform has unbalanced parentheses".to_string(),
            });
        };
        let expression = after_open[..close].trim();
        let (replacement, transformed_covariates) =
            formula_numeric_transform_term(transform, expression, numeric_covariates)?;
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

#[derive(Clone, Copy, Debug)]
struct FormulaNumericTransform {
    prefix: &'static str,
    label: &'static str,
    apply: fn(f64) -> f64,
}

const FORMULA_NUMERIC_TRANSFORMS: [FormulaNumericTransform; 7] = [
    FormulaNumericTransform {
        prefix: "poly(",
        label: "poly",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        prefix: "scale(",
        label: "scale",
        apply: std::convert::identity,
    },
    FormulaNumericTransform {
        prefix: "log10(",
        label: "log10",
        apply: f64::log10,
    },
    FormulaNumericTransform {
        prefix: "log2(",
        label: "log2",
        apply: f64::log2,
    },
    FormulaNumericTransform {
        prefix: "sqrt(",
        label: "sqrt",
        apply: f64::sqrt,
    },
    FormulaNumericTransform {
        prefix: "log(",
        label: "log",
        apply: f64::ln,
    },
    FormulaNumericTransform {
        prefix: "I(",
        label: "I",
        apply: std::convert::identity,
    },
];

fn next_formula_numeric_transform(rhs: &str) -> Option<(usize, FormulaNumericTransform)> {
    FORMULA_NUMERIC_TRANSFORMS
        .iter()
        .filter_map(|transform| rhs.find(transform.prefix).map(|idx| (idx, *transform)))
        .min_by_key(|(idx, _)| *idx)
}

fn formula_numeric_transform_term(
    transform: FormulaNumericTransform,
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    match transform.label {
        "I" => {
            let (name, values) = formula_numeric_power_term(expression, numeric_covariates)?;
            Ok((name.clone(), vec![(name, values)]))
        }
        "poly" => formula_numeric_raw_poly_term(expression, numeric_covariates),
        "scale" => {
            let (name, values) = formula_numeric_scale_term(expression, numeric_covariates)?;
            Ok((name.clone(), vec![(name, values)]))
        }
        _ => {
            let (name, values) =
                formula_numeric_function_term(transform, expression, numeric_covariates)?;
            Ok((name.clone(), vec![(name, values)]))
        }
    }
}

fn formula_numeric_function_term(
    transform: FormulaNumericTransform,
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    validate_formula_variable(expression)?;
    let covariate = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == expression)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{expression}' is not present in supplied design metadata"
            ),
        })?;
    let mut values = Vec::with_capacity(covariate.values.len());
    for (idx, value) in covariate.values.iter().copied().enumerate() {
        let transformed = (transform.apply)(value);
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
    Ok((format!("{}_{}", expression, transform.label), values))
}

fn formula_numeric_scale_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    let Some(numeric_name) = arguments.first().map(|value| value.trim()) else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' must provide a numeric covariate"
            ),
        });
    };
    validate_formula_variable(numeric_name)?;
    let mut center = FormulaScaleOption::Auto;
    let mut scale = FormulaScaleOption::Auto;
    for argument in arguments.iter().skip(1) {
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if let Some(value) = normalized.strip_prefix("center=") {
            center = parse_formula_scale_option(value, "center", expression)?;
            continue;
        }
        if let Some(value) = normalized.strip_prefix("scale=") {
            scale = parse_formula_scale_option(value, "scale", expression)?;
            continue;
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let covariate = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == numeric_name)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' is not present in supplied design metadata"
            ),
        })?;
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
    match value {
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

fn formula_numeric_raw_poly_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    if arguments.len() < 3 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' must provide numeric, degree, and raw=TRUE"
            ),
        });
    }
    let numeric_name = arguments[0].trim();
    validate_formula_variable(numeric_name)?;
    let mut degree_text = None;
    let mut has_raw_true = false;
    for (idx, argument) in arguments.iter().enumerate().skip(1) {
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if normalized == "raw=true" || normalized == "raw=t" {
            has_raw_true = true;
            continue;
        }
        if let Some(named_degree) = normalized.strip_prefix("degree=") {
            degree_text = Some(named_degree.to_string());
            continue;
        }
        if idx == 1 && degree_text.is_none() {
            degree_text = Some(argument.trim().to_string());
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
    let degree = degree_text
        .parse::<i32>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula raw polynomial degree '{degree_text}' must be an integer"),
        })?;
    if !(1..=16).contains(&degree) {
        return Err(DeseqError::InvalidOptions {
            reason: "formula raw polynomial degrees must be integers from 1 through 16".to_string(),
        });
    }
    if !has_raw_true {
        return Err(DeseqError::InvalidOptions {
            reason: "formula poly() transforms require raw=TRUE".to_string(),
        });
    }
    let covariate = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == numeric_name)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' is not present in supplied design metadata"
            ),
        })?;
    let mut replacement_terms = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for exponent in 1..=degree {
        let name = format!("{numeric_name}_poly_{exponent}");
        let mut values = Vec::with_capacity(covariate.values.len());
        for (idx, value) in covariate.values.iter().copied().enumerate() {
            let transformed = value.powi(exponent);
            if !transformed.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'poly({expression})' produced non-finite value at sample {idx}"
                    ),
                });
            }
            values.push(transformed);
        }
        replacement_terms.push(name.clone());
        derived.push((name, values));
    }
    Ok((format!("({})", replacement_terms.join(" + ")), derived))
}

fn split_formula_transform_arguments(expression: &str) -> Result<Vec<String>, DeseqError> {
    split_formula_top_level(expression, ',')
}

fn formula_numeric_power_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let Some((numeric_name, exponent_text)) = expression.split_once('^') else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' is not supported"),
        });
    };
    let numeric_name = numeric_name.trim();
    validate_formula_variable(numeric_name)?;
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
    let covariate = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == numeric_name)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' is not present in supplied design metadata"
            ),
        })?;
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
