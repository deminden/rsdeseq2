use crate::errors::{invalid_dimensions, DeseqError};
use crate::math::qr::{matrix_rank, DEFAULT_RANK_TOLERANCE};
use crate::matrix::RowMajorMatrix;

/// A design matrix generated outside the Rust core or by the Rust design helpers.
///
/// Rows are samples and columns are coefficients. The Rust formula helper covers
/// a primitive DESeq2-style subset; callers can still provide arbitrary matrices
/// directly when they need richer formula semantics.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignMatrix {
    matrix: RowMajorMatrix<f64>,
    coefficient_names: Option<Vec<String>>,
}

/// Expanded and standard design matrices for one factor.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorDesign {
    /// Expanded design with an intercept and one indicator column per level.
    pub expanded_design: DesignMatrix,
    /// Treatment-style standard design with an intercept and non-reference levels.
    pub standard_design: DesignMatrix,
    /// Expanded coefficient columns mapped to each standard-design coefficient.
    pub coefficient_groups: Vec<Vec<usize>>,
    /// Levels in first-observed order, with the reference first when present.
    pub levels: Vec<String>,
}

/// Caller-supplied metadata for one additive factor in an expanded design.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorSpec<'a> {
    /// Factor name used to build coefficient names.
    pub factor: &'a str,
    /// Per-sample factor levels in count-column order.
    pub sample_levels: &'a [String],
    /// Reference level for treatment-style reported coefficients.
    pub reference: &'a str,
}

/// Caller-supplied metadata for one additive numeric covariate.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedNumericSpec<'a> {
    /// Covariate name used as the coefficient name.
    pub name: &'a str,
    /// Per-sample finite numeric values in count-column order.
    pub values: &'a [f64],
}

/// Caller-supplied metadata for one additive factor-by-factor interaction.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorInteractionSpec<'a> {
    /// Left factor name, matching one [`ExpandedFactorSpec::factor`].
    pub left_factor: &'a str,
    /// Right factor name, matching one [`ExpandedFactorSpec::factor`].
    pub right_factor: &'a str,
}

/// Caller-supplied metadata for one factor-by-numeric interaction.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorNumericInteractionSpec<'a> {
    /// Factor name, matching one [`ExpandedFactorSpec::factor`].
    pub factor: &'a str,
    /// Numeric covariate name, matching one [`ExpandedNumericSpec::name`].
    pub numeric: &'a str,
}

/// Caller-supplied metadata for one numeric-by-numeric interaction.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedNumericInteractionSpec<'a> {
    /// Left numeric covariate name, matching one [`ExpandedNumericSpec::name`].
    pub left_numeric: &'a str,
    /// Right numeric covariate name, matching one [`ExpandedNumericSpec::name`].
    pub right_numeric: &'a str,
}

/// Expanded and standard design matrices for additive factor terms.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveFactorDesign {
    /// Expanded design with an intercept and one indicator column per factor level.
    pub expanded_design: DesignMatrix,
    /// Treatment-style standard design with an intercept and non-reference levels.
    pub standard_design: DesignMatrix,
    /// Expanded coefficient columns mapped to each standard-design coefficient.
    pub coefficient_groups: Vec<Vec<usize>>,
    /// Levels for each factor, in input factor order.
    pub factor_levels: Vec<Vec<String>>,
    /// Numeric covariate names in input order.
    pub numeric_covariates: Vec<String>,
    /// Factor-by-factor interaction names in input order.
    pub interactions: Vec<String>,
    /// Factor-by-numeric interaction names in input order.
    pub factor_numeric_interactions: Vec<String>,
    /// Numeric-by-numeric interaction names in input order.
    pub numeric_interactions: Vec<String>,
    /// Formula-only higher-order interaction names in input order.
    pub higher_order_interactions: Vec<String>,
}

/// Expanded formula design plus per-sample log-scale formula offsets.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFormulaDesignWithOffsets {
    /// Parsed expanded and reported design surfaces.
    pub design: ExpandedAdditiveFactorDesign,
    /// Per-sample offset values from supported `offset(numeric)` terms.
    pub offsets: Vec<f64>,
}

impl DesignMatrix {
    /// Create an intercept-only design matrix with one all-ones coefficient.
    pub fn intercept_only(n_samples: usize) -> Result<Self, DeseqError> {
        Self::from_row_major(
            n_samples,
            1,
            vec![1.0; n_samples],
            Some(vec!["Intercept".to_string()]),
        )
    }

    /// Create a design matrix from row-major values.
    pub fn from_row_major(
        n_samples: usize,
        n_coefficients: usize,
        values: Vec<f64>,
        coefficient_names: Option<Vec<String>>,
    ) -> Result<Self, DeseqError> {
        if let Some(names) = &coefficient_names {
            if names.len() != n_coefficients {
                return Err(invalid_dimensions(
                    "design coefficient names",
                    n_coefficients,
                    names.len(),
                ));
            }
        }
        let matrix = RowMajorMatrix::from_row_major(n_samples, n_coefficients, values)?;
        matrix.validate_finite("design matrix")?;
        Ok(Self {
            matrix,
            coefficient_names,
        })
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.matrix.n_rows()
    }

    /// Number of coefficients.
    pub fn n_coefficients(&self) -> usize {
        self.matrix.n_cols()
    }

    /// Reusable sample-index span.
    pub fn sample_indices(&self) -> core::range::Range<usize> {
        self.matrix.row_indices()
    }

    /// Reusable coefficient-index span.
    pub fn coefficient_indices(&self) -> core::range::Range<usize> {
        self.matrix.col_indices()
    }

    /// Return a contiguous sample-row block in row-major order.
    ///
    /// The range accepts both legacy range syntax (`1..3`) and the newer
    /// `core::range` types.
    pub fn sample_rows<R: core::ops::RangeBounds<usize>>(
        &self,
        samples: R,
    ) -> Result<&[f64], DeseqError> {
        self.matrix.rows(samples)
    }

    /// Matrix values.
    pub fn matrix(&self) -> &RowMajorMatrix<f64> {
        &self.matrix
    }

    /// Optional coefficient names.
    pub fn coefficient_names(&self) -> Option<&[String]> {
        self.coefficient_names.as_deref()
    }

    /// Resolve a coefficient name to its zero-based column index.
    pub fn coefficient_index(&self, name: &str) -> Result<usize, DeseqError> {
        let names = self
            .coefficient_names()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "coefficient names are required to resolve coefficient names".to_string(),
            })?;
        names
            .iter()
            .position(|candidate| candidate == name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("coefficient '{name}' is not present in coefficient names"),
            })
    }

    /// Numerical rank using the default deterministic tolerance.
    pub fn rank(&self) -> Result<usize, DeseqError> {
        self.rank_with_tolerance(DEFAULT_RANK_TOLERANCE)
    }

    /// Numerical rank using a caller-supplied absolute tolerance.
    pub fn rank_with_tolerance(&self, tolerance: f64) -> Result<usize, DeseqError> {
        matrix_rank(
            self.matrix.as_slice(),
            self.n_samples(),
            self.n_coefficients(),
            tolerance,
        )
    }

    /// Whether the design has full column rank using the default tolerance.
    pub fn is_full_rank(&self) -> Result<bool, DeseqError> {
        self.is_full_rank_with_tolerance(DEFAULT_RANK_TOLERANCE)
    }

    /// Whether the design has full column rank using a caller-supplied tolerance.
    pub fn is_full_rank_with_tolerance(&self, tolerance: f64) -> Result<bool, DeseqError> {
        Ok(self.rank_with_tolerance(tolerance)? == self.n_coefficients())
    }

    /// Return an error if the design is not full column rank.
    pub fn validate_full_rank(&self, context: &str) -> Result<(), DeseqError> {
        self.validate_full_rank_with_tolerance(context, DEFAULT_RANK_TOLERANCE)
    }

    /// Return an error if the design is not full column rank.
    pub fn validate_full_rank_with_tolerance(
        &self,
        context: &str,
        tolerance: f64,
    ) -> Result<(), DeseqError> {
        let rank = self.rank_with_tolerance(tolerance)?;
        if rank == self.n_coefficients() {
            return Ok(());
        }
        let zero_columns = self.zero_columns(tolerance);
        let reason = if zero_columns.is_empty() {
            format!(
                "{context} model matrix is not full rank: rank {rank}, columns {}",
                self.n_coefficients()
            )
        } else {
            format!(
                "{context} model matrix is not full rank: rank {rank}, columns {}; zero columns: {}",
                self.n_coefficients(),
                zero_columns.join(", ")
            )
        };
        Err(DeseqError::InvalidOptions { reason })
    }

    fn zero_columns(&self, tolerance: f64) -> Vec<String> {
        let mut columns = Vec::new();
        for col in 0..self.n_coefficients() {
            let is_zero = (0..self.n_samples())
                .all(|row| self.matrix.get(row, col).copied().unwrap_or(0.0).abs() <= tolerance);
            if is_zero {
                columns.push(
                    self.coefficient_names()
                        .and_then(|names| names.get(col))
                        .cloned()
                        .unwrap_or_else(|| col.to_string()),
                );
            }
        }
        columns
    }
}

/// Build a primitive expanded one-factor design and matching standard design.
///
/// The expanded matrix has columns `Intercept`, then one indicator per factor
/// level. The standard matrix has `Intercept`, then one treatment-style
/// indicator for every non-reference level. This helper owns only the matrix
/// construction step; callers still choose how to use the returned
/// `coefficient_groups` in beta-prior fitting and result assembly.
pub fn expanded_factor_design<S: AsRef<str>>(
    factor: &str,
    sample_levels: &[S],
    reference: &str,
) -> Result<ExpandedFactorDesign, DeseqError> {
    validate_factor_design_inputs(factor, sample_levels, reference)?;
    let levels = ordered_levels(sample_levels, reference);
    let n_samples = sample_levels.len();

    let mut expanded_names = Vec::with_capacity(levels.len() + 1);
    expanded_names.push("Intercept".to_string());
    expanded_names.extend(
        levels
            .iter()
            .map(|level| expanded_factor_coefficient_name(factor, level)),
    );
    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for level in sample_levels {
        let level = level.as_ref();
        expanded_values.push(1.0);
        for candidate in &levels {
            expanded_values.push((level == candidate) as u8 as f64);
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names),
    )?;

    let non_reference_levels = levels
        .iter()
        .filter(|level| level.as_str() != reference)
        .collect::<Vec<_>>();
    let mut standard_names = Vec::with_capacity(non_reference_levels.len() + 1);
    standard_names.push("Intercept".to_string());
    standard_names.extend(
        non_reference_levels
            .iter()
            .map(|level| standard_factor_coefficient_name(factor, level, reference)),
    );
    let mut standard_values = Vec::with_capacity(n_samples * standard_names.len());
    for level in sample_levels {
        let level = level.as_ref();
        standard_values.push(1.0);
        for candidate in &non_reference_levels {
            standard_values.push((level == candidate.as_str()) as u8 as f64);
        }
    }
    let standard_design = DesignMatrix::from_row_major(
        n_samples,
        standard_names.len(),
        standard_values,
        Some(standard_names),
    )?;

    let mut coefficient_groups = Vec::with_capacity(standard_design.n_coefficients());
    coefficient_groups.push(vec![0]);
    for level in non_reference_levels {
        let expanded_column = levels
            .iter()
            .position(|candidate| candidate == level)
            .map(|idx| idx + 1)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("factor level '{}' is not present in expanded levels", level),
            })?;
        coefficient_groups.push(vec![expanded_column]);
    }

    Ok(ExpandedFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        levels,
    })
}

/// Build an expanded additive-factor design and matching standard design.
///
/// This covers primitive `~ factor1 + factor2 + ...` model-matrix
/// construction for categorical terms. The expanded matrix has one intercept
/// plus one indicator per level of each factor. The standard matrix has one
/// intercept plus one treatment-style indicator per non-reference level of
/// each factor. Interactions, nested terms, continuous covariates, and formula
/// parsing remain caller or wrapper responsibilities.
pub fn expanded_additive_factor_design(
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design(factors, &[])
}

/// Build an expanded additive design with categorical factors and numeric covariates.
///
/// Categorical terms use expanded indicator columns and treatment-style
/// reported columns as in [`expanded_additive_factor_design`]. Numeric
/// covariates are included unchanged in both expanded and standard designs and
/// are mapped one-to-one in the coefficient groups. Interactions, splines,
/// nested terms, transformed variables, and formula parsing remain caller or
/// wrapper responsibilities.
pub fn expanded_additive_design(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design_with_interactions(factors, numeric_covariates, &[])
}

/// Build an expanded additive design with categorical factors, numerics, and factor interactions.
///
/// Factor-by-factor interactions add all level-pair products to the expanded
/// design and non-reference-by-non-reference treatment-style products to the
/// standard design. The interaction terms reference factor names already
/// supplied in `factors`.
pub fn expanded_additive_design_with_interactions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedFactorInteractionSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design_with_all_interactions(
        factors,
        numeric_covariates,
        interactions,
        &[],
        &[],
    )
}

/// Build an expanded additive design with all primitive pairwise interaction types.
pub fn expanded_additive_design_with_all_interactions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    factor_interactions: &[ExpandedFactorInteractionSpec<'_>],
    factor_numeric_interactions: &[ExpandedFactorNumericInteractionSpec<'_>],
    numeric_interactions: &[ExpandedNumericInteractionSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    let n_samples = additive_design_sample_count(factors, numeric_covariates)?;
    if !factors.is_empty() {
        validate_additive_factor_specs(factors)?;
    }
    validate_numeric_covariate_specs(numeric_covariates, n_samples)?;
    let interaction_indices = resolve_interaction_indices(factors, factor_interactions)?;
    let factor_numeric_indices = resolve_factor_numeric_interaction_indices(
        factors,
        numeric_covariates,
        factor_numeric_interactions,
    )?;
    let numeric_interaction_indices =
        resolve_numeric_interaction_indices(numeric_covariates, numeric_interactions)?;
    let factor_levels = factors
        .iter()
        .map(|factor| ordered_levels(factor.sample_levels, factor.reference))
        .collect::<Vec<_>>();

    let mut expanded_names = vec!["Intercept".to_string()];
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        expanded_names.extend(
            levels
                .iter()
                .map(|level| expanded_factor_coefficient_name(factor.factor, level)),
        );
    }
    expanded_names.extend(
        numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for (left_idx, right_idx) in &interaction_indices {
        let left = &factors[*left_idx];
        let right = &factors[*right_idx];
        for left_level in &factor_levels[*left_idx] {
            for right_level in &factor_levels[*right_idx] {
                expanded_names.push(interaction_coefficient_name(
                    left.factor,
                    left_level,
                    right.factor,
                    right_level,
                ));
            }
        }
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let factor = &factors[*factor_idx];
        let numeric = &numeric_covariates[*numeric_idx];
        for level in &factor_levels[*factor_idx] {
            expanded_names.push(factor_numeric_interaction_coefficient_name(
                factor.factor,
                level,
                numeric.name,
            ));
        }
    }
    for (left_idx, right_idx) in &numeric_interaction_indices {
        expanded_names.push(numeric_interaction_coefficient_name(
            numeric_covariates[*left_idx].name,
            numeric_covariates[*right_idx].name,
        ));
    }
    validate_unique_coefficient_names(&expanded_names, "expanded additive design")?;
    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for sample in 0..n_samples {
        expanded_values.push(1.0);
        for (factor, levels) in factors.iter().zip(&factor_levels) {
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels {
                expanded_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in numeric_covariates {
            expanded_values.push(covariate.values[sample]);
        }
        for (left_idx, right_idx) in &interaction_indices {
            let left = &factors[*left_idx];
            let right = &factors[*right_idx];
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in &factor_levels[*left_idx] {
                for right_level in &factor_levels[*right_idx] {
                    expanded_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for (factor_idx, numeric_idx) in &factor_numeric_indices {
            let factor = &factors[*factor_idx];
            let numeric = &numeric_covariates[*numeric_idx];
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in &factor_levels[*factor_idx] {
                expanded_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for (left_idx, right_idx) in &numeric_interaction_indices {
            expanded_values.push(
                numeric_covariates[*left_idx].values[sample]
                    * numeric_covariates[*right_idx].values[sample],
            );
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names),
    )?;

    let mut standard_names = vec!["Intercept".to_string()];
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        standard_names.extend(
            levels
                .iter()
                .filter(|level| level.as_str() != factor.reference)
                .map(|level| {
                    standard_factor_coefficient_name(factor.factor, level, factor.reference)
                }),
        );
    }
    standard_names.extend(
        numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for (left_idx, right_idx) in &interaction_indices {
        let left = &factors[*left_idx];
        let right = &factors[*right_idx];
        for left_level in factor_levels[*left_idx]
            .iter()
            .filter(|level| level.as_str() != left.reference)
        {
            for right_level in factor_levels[*right_idx]
                .iter()
                .filter(|level| level.as_str() != right.reference)
            {
                standard_names.push(standard_interaction_coefficient_name(
                    left.factor,
                    left_level,
                    left.reference,
                    right.factor,
                    right_level,
                    right.reference,
                ));
            }
        }
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let factor = &factors[*factor_idx];
        let numeric = &numeric_covariates[*numeric_idx];
        for level in factor_levels[*factor_idx]
            .iter()
            .filter(|level| level.as_str() != factor.reference)
        {
            standard_names.push(standard_factor_numeric_interaction_coefficient_name(
                factor.factor,
                level,
                factor.reference,
                numeric.name,
            ));
        }
    }
    for (left_idx, right_idx) in &numeric_interaction_indices {
        standard_names.push(numeric_interaction_coefficient_name(
            numeric_covariates[*left_idx].name,
            numeric_covariates[*right_idx].name,
        ));
    }
    validate_unique_coefficient_names(&standard_names, "standard additive design")?;
    let mut standard_values = Vec::with_capacity(n_samples * standard_names.len());
    for sample in 0..n_samples {
        standard_values.push(1.0);
        for (factor, levels) in factors.iter().zip(&factor_levels) {
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels
                .iter()
                .filter(|level| level.as_str() != factor.reference)
            {
                standard_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in numeric_covariates {
            standard_values.push(covariate.values[sample]);
        }
        for (left_idx, right_idx) in &interaction_indices {
            let left = &factors[*left_idx];
            let right = &factors[*right_idx];
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in factor_levels[*left_idx]
                .iter()
                .filter(|level| level.as_str() != left.reference)
            {
                for right_level in factor_levels[*right_idx]
                    .iter()
                    .filter(|level| level.as_str() != right.reference)
                {
                    standard_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for (factor_idx, numeric_idx) in &factor_numeric_indices {
            let factor = &factors[*factor_idx];
            let numeric = &numeric_covariates[*numeric_idx];
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in factor_levels[*factor_idx]
                .iter()
                .filter(|level| level.as_str() != factor.reference)
            {
                standard_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for (left_idx, right_idx) in &numeric_interaction_indices {
            standard_values.push(
                numeric_covariates[*left_idx].values[sample]
                    * numeric_covariates[*right_idx].values[sample],
            );
        }
    }
    let standard_design = DesignMatrix::from_row_major(
        n_samples,
        standard_names.len(),
        standard_values,
        Some(standard_names),
    )?;

    let mut coefficient_groups = Vec::with_capacity(standard_design.n_coefficients());
    coefficient_groups.push(vec![0]);
    let mut expanded_offset = 1;
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        for (level_idx, level) in levels.iter().enumerate() {
            if level.as_str() != factor.reference {
                coefficient_groups.push(vec![expanded_offset + level_idx]);
            }
        }
        expanded_offset += levels.len();
    }
    for idx in 0..numeric_covariates.len() {
        coefficient_groups.push(vec![expanded_offset + idx]);
    }
    expanded_offset += numeric_covariates.len();
    for (left_idx, right_idx) in &interaction_indices {
        let left_levels = &factor_levels[*left_idx];
        let right_levels = &factor_levels[*right_idx];
        for (left_level_idx, left_level) in left_levels.iter().enumerate() {
            if left_level.as_str() == factors[*left_idx].reference {
                continue;
            }
            for (right_level_idx, right_level) in right_levels.iter().enumerate() {
                if right_level.as_str() == factors[*right_idx].reference {
                    continue;
                }
                let expanded_column =
                    expanded_offset + left_level_idx * right_levels.len() + right_level_idx;
                coefficient_groups.push(vec![expanded_column]);
            }
        }
        expanded_offset += left_levels.len() * right_levels.len();
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let levels = &factor_levels[*factor_idx];
        let factor = &factors[*factor_idx];
        for (level_idx, level) in levels.iter().enumerate() {
            if level.as_str() != factor.reference {
                coefficient_groups.push(vec![expanded_offset + level_idx]);
            }
        }
        let _ = numeric_idx;
        expanded_offset += levels.len();
    }
    for idx in 0..numeric_interaction_indices.len() {
        coefficient_groups.push(vec![expanded_offset + idx]);
    }

    Ok(ExpandedAdditiveFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        factor_levels,
        numeric_covariates: numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string())
            .collect(),
        interactions: factor_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.left_factor, interaction.right_factor))
            .collect(),
        factor_numeric_interactions: factor_numeric_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.factor, interaction.numeric))
            .collect(),
        numeric_interactions: numeric_interactions
            .iter()
            .map(|interaction| {
                format!("{}:{}", interaction.left_numeric, interaction.right_numeric)
            })
            .collect(),
        higher_order_interactions: Vec::new(),
    })
}

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

fn expand_parenthesized_formula_terms(rhs: &str) -> Result<String, DeseqError> {
    let signed_terms = split_formula_signed_terms(rhs)?;
    let mut expanded = Vec::new();
    for (sign, term) in signed_terms {
        for expanded_term in expand_parenthesized_formula_term(&term)? {
            if sign < 0 {
                expanded.push(format!("- {expanded_term}"));
            } else {
                expanded.push(expanded_term);
            }
        }
    }
    Ok(join_formula_terms(&expanded))
}

fn join_formula_terms(terms: &[String]) -> String {
    let mut joined = String::new();
    for term in terms {
        if joined.is_empty() {
            joined.push_str(term);
        } else if term.trim_start().starts_with('-') {
            joined.push(' ');
            joined.push_str(term);
        } else {
            joined.push_str(" + ");
            joined.push_str(term);
        }
    }
    joined
}

fn split_formula_signed_terms(rhs: &str) -> Result<Vec<(i8, String)>, DeseqError> {
    let mut terms = Vec::new();
    let mut depth = 0_i32;
    let mut sign = 1_i8;
    let mut start = 0_usize;
    let mut saw_term = false;
    for (idx, character) in rhs.char_indices() {
        match character {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return Err(DeseqError::InvalidOptions {
                        reason: "formula parentheses are unbalanced".to_string(),
                    });
                }
            }
            '+' | '-' if depth == 0 => {
                let term = rhs[start..idx].trim();
                if term.is_empty() {
                    if saw_term || character == '+' {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula contains an empty term".to_string(),
                        });
                    }
                    sign = -1;
                    start = idx + character.len_utf8();
                    continue;
                }
                terms.push((sign, term.to_string()));
                saw_term = true;
                sign = if character == '-' { -1 } else { 1 };
                start = idx + character.len_utf8();
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula parentheses are unbalanced".to_string(),
        });
    }
    let term = rhs[start..].trim();
    if term.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula contains an empty term".to_string(),
        });
    }
    terms.push((sign, term.to_string()));
    Ok(terms)
}

fn expand_parenthesized_formula_term(term: &str) -> Result<Vec<String>, DeseqError> {
    let term = strip_formula_outer_parentheses(term.trim())?;
    if term.contains('-') {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' contains unsupported nested subtraction"),
        });
    }
    let additive_pieces = split_formula_top_level(term, '+')?;
    if additive_pieces.len() > 1 {
        return split_formula_additive_group(term);
    }
    for delimiter in ['*', ':', '/'] {
        let pieces = split_formula_top_level(term, delimiter)?;
        if pieces.len() > 1 {
            return expand_parenthesized_formula_operator(&pieces, delimiter);
        }
    }
    if term.contains('(') || term.contains(')') {
        return split_formula_additive_group(term);
    }
    Ok(vec![term.to_string()])
}

fn expand_parenthesized_formula_operator(
    pieces: &[String],
    delimiter: char,
) -> Result<Vec<String>, DeseqError> {
    let mut alternatives = Vec::with_capacity(pieces.len());
    for piece in pieces {
        alternatives.push(split_formula_additive_group(piece)?);
    }
    match delimiter {
        '*' => expand_parenthesized_star(&alternatives),
        ':' => Ok(formula_alternative_products(&alternatives)
            .into_iter()
            .map(|product| product.join(":"))
            .collect()),
        '/' => {
            let mut terms = Vec::new();
            for product in formula_alternative_products(&alternatives) {
                for prefix_len in 1..=product.len() {
                    push_unique_formula_term(&mut terms, product[..prefix_len].join(":"));
                }
            }
            Ok(terms)
        }
        _ => unreachable!("unsupported formula delimiter"),
    }
}

fn expand_parenthesized_star(alternatives: &[Vec<String>]) -> Result<Vec<String>, DeseqError> {
    let mut terms = Vec::new();
    for group in alternatives {
        for term in group {
            push_unique_formula_term(&mut terms, term.clone());
        }
    }
    for order in 2..=alternatives.len() {
        for group_subset in formula_group_combinations(alternatives, order) {
            for product in formula_alternative_products(&group_subset) {
                push_unique_formula_term(&mut terms, product.join(":"));
            }
        }
    }
    Ok(terms)
}

fn split_formula_additive_group(term: &str) -> Result<Vec<String>, DeseqError> {
    let stripped = strip_formula_outer_parentheses(term.trim())?;
    let pieces = split_formula_top_level(stripped, '+')?;
    if pieces.len() == 1 {
        if stripped.contains('-') {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula term '{term}' contains unsupported nested subtraction"),
            });
        }
        if stripped.contains('(') || stripped.contains(')') {
            return expand_parenthesized_formula_term(stripped);
        }
        return Ok(vec![stripped.to_string()]);
    }
    let mut terms = Vec::new();
    for piece in pieces {
        if piece.contains('-') {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula group '{term}' contains unsupported nested subtraction"),
            });
        }
        for expanded in expand_parenthesized_formula_term(&piece)? {
            push_unique_formula_term(&mut terms, expanded);
        }
    }
    Ok(terms)
}

fn strip_formula_outer_parentheses(term: &str) -> Result<&str, DeseqError> {
    let mut stripped = term.trim();
    loop {
        if !(stripped.starts_with('(') && stripped.ends_with(')')) {
            return Ok(stripped);
        }
        let mut depth = 0_i32;
        let mut encloses_whole_term = true;
        for (idx, character) in stripped.char_indices() {
            match character {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth < 0 {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula parentheses are unbalanced".to_string(),
                        });
                    }
                    if depth == 0 && idx + character.len_utf8() < stripped.len() {
                        encloses_whole_term = false;
                        break;
                    }
                }
                _ => {}
            }
        }
        if depth != 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "formula parentheses are unbalanced".to_string(),
            });
        }
        if !encloses_whole_term {
            return Ok(stripped);
        }
        stripped = stripped[1..stripped.len() - 1].trim();
    }
}

fn split_formula_top_level(term: &str, delimiter: char) -> Result<Vec<String>, DeseqError> {
    let mut pieces = Vec::new();
    let mut depth = 0_i32;
    let mut start = 0_usize;
    for (idx, character) in term.char_indices() {
        match character {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return Err(DeseqError::InvalidOptions {
                        reason: "formula parentheses are unbalanced".to_string(),
                    });
                }
            }
            _ if character == delimiter && depth == 0 => {
                let piece = term[start..idx].trim();
                if piece.is_empty() {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!("formula term '{term}' contains an empty component"),
                    });
                }
                pieces.push(piece.to_string());
                start = idx + character.len_utf8();
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula parentheses are unbalanced".to_string(),
        });
    }
    let piece = term[start..].trim();
    if piece.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' contains an empty component"),
        });
    }
    pieces.push(piece.to_string());
    Ok(pieces)
}

fn formula_alternative_products(alternatives: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut products: Vec<Vec<String>> = vec![Vec::new()];
    for group in alternatives {
        let mut next = Vec::new();
        for prefix in &products {
            for term in group {
                let mut product = prefix.clone();
                product.push(term.clone());
                next.push(product);
            }
        }
        products = next;
    }
    products
}

fn formula_group_combinations(groups: &[Vec<String>], order: usize) -> Vec<Vec<Vec<String>>> {
    fn push_group_combinations(
        groups: &[Vec<String>],
        order: usize,
        start: usize,
        current: &mut Vec<Vec<String>>,
        output: &mut Vec<Vec<Vec<String>>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=groups.len() - remaining {
            current.push(groups[idx].clone());
            push_group_combinations(groups, order, idx + 1, current, output);
            current.pop();
        }
    }
    if order == 0 || order > groups.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_group_combinations(groups, order, 0, &mut Vec::new(), &mut output);
    output
}

fn push_unique_formula_term(terms: &mut Vec<String>, term: String) {
    if !terms.iter().any(|candidate| candidate == &term) {
        terms.push(term);
    }
}

fn remove_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if term.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula subtraction must be followed by a term".to_string(),
        });
    }
    if term == "1" {
        state.has_intercept = false;
        return Ok(());
    }
    if term == "0" || term == "-1" {
        return Ok(());
    }
    if term.contains('^') || term.contains('(') || term.contains(')') {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' is not supported by the primitive parser"),
        });
    }
    if term.contains('/') {
        remove_nested_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if term.contains('*') {
        remove_star_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if term.contains(':') {
        remove_interaction_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    remove_main_formula_term(term, factors, numeric_covariates, state)
}

fn add_star_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '*')?;
    for piece in &pieces {
        add_main_formula_term(
            piece,
            factors,
            numeric_covariates,
            &mut state.selected_factors,
            &mut state.selected_numeric_covariates,
        )?;
    }
    for left in 0..pieces.len() {
        for right in (left + 1)..pieces.len() {
            add_pairwise_formula_interaction(
                pieces[left],
                pieces[right],
                term,
                factors,
                numeric_covariates,
                state,
            )?;
        }
    }
    for order in 3..=pieces.len() {
        for subset in formula_piece_combinations(&pieces, order) {
            add_higher_order_formula_interaction(
                &subset,
                term,
                factors,
                numeric_covariates,
                state,
            )?;
        }
    }
    Ok(())
}

fn add_nested_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '/')?;
    add_main_formula_term(
        pieces[0],
        factors,
        numeric_covariates,
        &mut state.selected_factors,
        &mut state.selected_numeric_covariates,
    )?;
    add_pairwise_formula_interaction(
        pieces[0],
        pieces[1],
        term,
        factors,
        numeric_covariates,
        state,
    )?;
    for prefix_len in 3..=pieces.len() {
        add_higher_order_formula_interaction(
            &pieces[..prefix_len],
            term,
            factors,
            numeric_covariates,
            state,
        )?;
    }
    Ok(())
}

fn remove_star_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '*')?;
    for piece in &pieces {
        remove_main_formula_term(piece, factors, numeric_covariates, state)?;
    }
    for left in 0..pieces.len() {
        for right in (left + 1)..pieces.len() {
            remove_pairwise_formula_interaction(
                pieces[left],
                pieces[right],
                factors,
                numeric_covariates,
                state,
            )?;
        }
    }
    for order in 3..=pieces.len() {
        for subset in formula_piece_combinations(&pieces, order) {
            remove_higher_order_formula_interaction(&subset, factors, numeric_covariates, state)?;
        }
    }
    Ok(())
}

fn remove_nested_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '/')?;
    remove_main_formula_term(pieces[0], factors, numeric_covariates, state)?;
    remove_pairwise_formula_interaction(pieces[0], pieces[1], factors, numeric_covariates, state)?;
    for prefix_len in 3..=pieces.len() {
        remove_higher_order_formula_interaction(
            &pieces[..prefix_len],
            factors,
            numeric_covariates,
            state,
        )?;
    }
    Ok(())
}

fn add_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    selected_factors: &mut Vec<ExpandedFactorSpec<'a>>,
    selected_numeric_covariates: &mut Vec<ExpandedNumericSpec<'a>>,
) -> Result<(), DeseqError> {
    validate_formula_variable(term)?;
    if let Some(factor) = factors.iter().find(|candidate| candidate.factor == term) {
        if selected_factors
            .iter()
            .any(|candidate| candidate.factor == term)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula main effect '{term}' appears more than once"),
            });
        }
        selected_factors.push(factor.clone());
        return Ok(());
    }
    if let Some(covariate) = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == term)
    {
        if selected_numeric_covariates
            .iter()
            .any(|candidate| candidate.name == term)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula main effect '{term}' appears more than once"),
            });
        }
        selected_numeric_covariates.push(covariate.clone());
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{term}' is not present in supplied design metadata"),
    })
}

fn remove_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(term)?;
    if factors.iter().any(|candidate| candidate.factor == term) {
        state
            .selected_factors
            .retain(|candidate| candidate.factor != term);
        return Ok(());
    }
    if numeric_covariates
        .iter()
        .any(|candidate| candidate.name == term)
    {
        state
            .selected_numeric_covariates
            .retain(|candidate| candidate.name != term);
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{term}' is not present in supplied design metadata"),
    })
}

fn add_interaction_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, ':')?;
    if pieces.len() == 2 {
        add_pairwise_formula_interaction(
            pieces[0],
            pieces[1],
            term,
            factors,
            numeric_covariates,
            state,
        )
    } else {
        add_higher_order_formula_interaction(&pieces, term, factors, numeric_covariates, state)
    }
}

fn add_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    display_term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(left)?;
    validate_formula_variable(right)?;
    if left == right {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{display_term}' cannot use one variable twice"),
        });
    }

    let left_factor = factors.iter().find(|candidate| candidate.factor == left);
    let right_factor = factors.iter().find(|candidate| candidate.factor == right);
    let left_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == left);
    let right_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == right);

    match (left_factor, right_factor, left_numeric, right_numeric) {
        (Some(left_factor), Some(right_factor), None, None) => {
            if state.factor_interactions.iter().any(|interaction| {
                (interaction.left_factor == left_factor.factor
                    && interaction.right_factor == right_factor.factor)
                    || (interaction.left_factor == right_factor.factor
                        && interaction.right_factor == left_factor.factor)
            }) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula interaction '{display_term}' appears more than once"),
                });
            }
            state
                .factor_interactions
                .push(ExpandedFactorInteractionSpec {
                    left_factor: left_factor.factor,
                    right_factor: right_factor.factor,
                });
        }
        (Some(factor), None, None, Some(numeric)) => {
            add_factor_numeric_formula_interaction(
                factor.factor,
                numeric.name,
                &mut state.factor_numeric_interactions,
                display_term,
            )?;
        }
        (None, Some(factor), Some(numeric), None) => {
            add_factor_numeric_formula_interaction(
                factor.factor,
                numeric.name,
                &mut state.factor_numeric_interactions,
                display_term,
            )?;
        }
        (None, None, Some(left_numeric), Some(right_numeric)) => {
            if state.numeric_interactions.iter().any(|interaction| {
                (interaction.left_numeric == left_numeric.name
                    && interaction.right_numeric == right_numeric.name)
                    || (interaction.left_numeric == right_numeric.name
                        && interaction.right_numeric == left_numeric.name)
            }) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula interaction '{display_term}' appears more than once"),
                });
            }
            state
                .numeric_interactions
                .push(ExpandedNumericInteractionSpec {
                    left_numeric: left_numeric.name,
                    right_numeric: right_numeric.name,
                });
        }
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{display_term}' references variables missing from supplied design metadata"
                ),
            });
        }
    }
    Ok(())
}

fn remove_interaction_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, ':')?;
    if pieces.len() == 2 {
        remove_pairwise_formula_interaction(
            pieces[0],
            pieces[1],
            factors,
            numeric_covariates,
            state,
        )
    } else {
        remove_higher_order_formula_interaction(&pieces, factors, numeric_covariates, state)
    }
}

fn remove_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(left)?;
    validate_formula_variable(right)?;
    if left == right {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{left}:{right}' cannot use one variable twice"),
        });
    }
    let left_variable = resolve_formula_variable(left, factors, numeric_covariates)?;
    let right_variable = resolve_formula_variable(right, factors, numeric_covariates)?;
    match (&left_variable, &right_variable) {
        (FormulaVariableRef::Factor(left_factor), FormulaVariableRef::Factor(right_factor)) => {
            state.factor_interactions.retain(|interaction| {
                !((interaction.left_factor == *left_factor
                    && interaction.right_factor == *right_factor)
                    || (interaction.left_factor == *right_factor
                        && interaction.right_factor == *left_factor))
            });
        }
        (FormulaVariableRef::Factor(factor), FormulaVariableRef::Numeric(numeric))
        | (FormulaVariableRef::Numeric(numeric), FormulaVariableRef::Factor(factor)) => {
            state.factor_numeric_interactions.retain(|interaction| {
                !(interaction.factor == *factor && interaction.numeric == *numeric)
            });
        }
        (FormulaVariableRef::Numeric(left_numeric), FormulaVariableRef::Numeric(right_numeric)) => {
            state.numeric_interactions.retain(|interaction| {
                !((interaction.left_numeric == *left_numeric
                    && interaction.right_numeric == *right_numeric)
                    || (interaction.left_numeric == *right_numeric
                        && interaction.right_numeric == *left_numeric))
            });
        }
    }
    Ok(())
}

fn add_higher_order_formula_interaction<'a>(
    pieces: &[&str],
    display_term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if pieces.len() < 3 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula interaction '{display_term}' must use at least three variables"
            ),
        });
    }
    let mut variables = Vec::with_capacity(pieces.len());
    for piece in pieces {
        validate_formula_variable(piece)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == *piece)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{display_term}' cannot use one variable twice"
                ),
            });
        }
        variables.push(resolve_formula_variable(
            piece,
            factors,
            numeric_covariates,
        )?);
    }
    if state
        .higher_order_interactions
        .iter()
        .any(|interaction| same_formula_variable_set(&interaction.variables, &variables))
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{display_term}' appears more than once"),
        });
    }
    state
        .higher_order_interactions
        .push(FormulaHigherOrderInteractionSpec { variables });
    Ok(())
}

fn remove_higher_order_formula_interaction<'a>(
    pieces: &[&str],
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if pieces.len() < 3 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula interaction removal must use at least three variables".to_string(),
        });
    }
    let mut variables = Vec::with_capacity(pieces.len());
    for piece in pieces {
        validate_formula_variable(piece)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == *piece)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{}' cannot use one variable twice",
                    pieces.join(":")
                ),
            });
        }
        variables.push(resolve_formula_variable(
            piece,
            factors,
            numeric_covariates,
        )?);
    }
    state
        .higher_order_interactions
        .retain(|interaction| !same_formula_variable_set(&interaction.variables, &variables));
    Ok(())
}

fn resolve_formula_variable<'a>(
    variable: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<FormulaVariableRef<'a>, DeseqError> {
    if let Some(factor) = factors
        .iter()
        .find(|candidate| candidate.factor == variable)
    {
        return Ok(FormulaVariableRef::Factor(factor.factor));
    }
    if let Some(covariate) = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == variable)
    {
        return Ok(FormulaVariableRef::Numeric(covariate.name));
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{variable}' is not present in supplied design metadata"),
    })
}

fn same_formula_variable_set(
    left: &[FormulaVariableRef<'_>],
    right: &[FormulaVariableRef<'_>],
) -> bool {
    left.len() == right.len()
        && left.iter().all(|left_variable| {
            right
                .iter()
                .any(|right_variable| right_variable == left_variable)
        })
}

fn expanded_formula_design_from_state<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    let used_factors = formula_used_factors(state, factors)?;
    let used_numeric_covariates = formula_used_numeric_covariates(state, numeric_covariates)?;
    let n_samples = formula_design_sample_count(
        state,
        &used_factors,
        &used_numeric_covariates,
        factors,
        numeric_covariates,
    )?;
    validate_formula_used_factor_specs(&used_factors, n_samples)?;
    validate_numeric_covariate_specs(&used_numeric_covariates, n_samples)?;

    let factor_levels = used_factors
        .iter()
        .map(|factor| ordered_levels(factor.sample_levels, factor.reference))
        .collect::<Vec<_>>();

    let mut expanded_names = Vec::new();
    if state.has_intercept {
        expanded_names.push("Intercept".to_string());
    }
    for factor in &state.selected_factors {
        let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
        expanded_names.extend(
            levels
                .iter()
                .map(|level| expanded_factor_coefficient_name(factor.factor, level)),
        );
    }
    expanded_names.extend(
        state
            .selected_numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for interaction in &state.factor_interactions {
        let left_levels =
            formula_factor_levels(interaction.left_factor, &used_factors, &factor_levels)?;
        let right_levels =
            formula_factor_levels(interaction.right_factor, &used_factors, &factor_levels)?;
        for left_level in left_levels {
            for right_level in right_levels {
                expanded_names.push(interaction_coefficient_name(
                    interaction.left_factor,
                    left_level,
                    interaction.right_factor,
                    right_level,
                ));
            }
        }
    }
    for interaction in &state.factor_numeric_interactions {
        let levels = formula_factor_levels(interaction.factor, &used_factors, &factor_levels)?;
        for level in levels {
            expanded_names.push(factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                interaction.numeric,
            ));
        }
    }
    for interaction in &state.numeric_interactions {
        expanded_names.push(numeric_interaction_coefficient_name(
            interaction.left_numeric,
            interaction.right_numeric,
        ));
    }
    for interaction in &state.higher_order_interactions {
        expanded_names.extend(formula_higher_order_expanded_names(
            interaction,
            &used_factors,
            &factor_levels,
        )?);
    }
    validate_unique_coefficient_names(&expanded_names, "expanded formula design")?;

    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for sample in 0..n_samples {
        if state.has_intercept {
            expanded_values.push(1.0);
        }
        for factor in &state.selected_factors {
            let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels {
                expanded_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in &state.selected_numeric_covariates {
            expanded_values.push(covariate.values[sample]);
        }
        for interaction in &state.factor_interactions {
            let left = formula_factor_spec(interaction.left_factor, &used_factors)?;
            let right = formula_factor_spec(interaction.right_factor, &used_factors)?;
            let left_levels =
                formula_factor_levels(interaction.left_factor, &used_factors, &factor_levels)?;
            let right_levels =
                formula_factor_levels(interaction.right_factor, &used_factors, &factor_levels)?;
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in left_levels {
                for right_level in right_levels {
                    expanded_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for interaction in &state.factor_numeric_interactions {
            let factor = formula_factor_spec(interaction.factor, &used_factors)?;
            let numeric = formula_numeric_spec(interaction.numeric, &used_numeric_covariates)?;
            let levels = formula_factor_levels(interaction.factor, &used_factors, &factor_levels)?;
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in levels {
                expanded_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for interaction in &state.numeric_interactions {
            expanded_values.push(
                formula_numeric_spec(interaction.left_numeric, &used_numeric_covariates)?.values
                    [sample]
                    * formula_numeric_spec(interaction.right_numeric, &used_numeric_covariates)?
                        .values[sample],
            );
        }
        for interaction in &state.higher_order_interactions {
            expanded_values.extend(formula_higher_order_expanded_values(
                interaction,
                &used_factors,
                &used_numeric_covariates,
                &factor_levels,
                sample,
            )?);
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names.clone()),
    )?;

    let mut standard_names = Vec::new();
    let mut standard_expanded_names = Vec::new();
    if state.has_intercept {
        standard_names.push("Intercept".to_string());
        standard_expanded_names.push("Intercept".to_string());
    }
    for (factor_idx, factor) in state.selected_factors.iter().enumerate() {
        let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
        let use_all_levels = !state.has_intercept && factor_idx == 0;
        for level in levels
            .iter()
            .filter(|level| use_all_levels || level.as_str() != factor.reference)
        {
            if use_all_levels {
                standard_names.push(expanded_factor_coefficient_name(factor.factor, level));
            } else {
                standard_names.push(standard_factor_coefficient_name(
                    factor.factor,
                    level,
                    factor.reference,
                ));
            }
            standard_expanded_names.push(expanded_factor_coefficient_name(factor.factor, level));
        }
    }
    for covariate in &state.selected_numeric_covariates {
        standard_names.push(covariate.name.to_string());
        standard_expanded_names.push(covariate.name.to_string());
    }
    for interaction in &state.factor_interactions {
        push_formula_factor_interaction_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_formula_factor_numeric_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    for interaction in &state.numeric_interactions {
        let name = numeric_interaction_coefficient_name(
            interaction.left_numeric,
            interaction.right_numeric,
        );
        standard_names.push(name.clone());
        standard_expanded_names.push(name);
    }
    for interaction in &state.higher_order_interactions {
        push_formula_higher_order_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    validate_unique_coefficient_names(&standard_names, "standard formula design")?;

    let mut standard_values = Vec::with_capacity(n_samples * standard_expanded_names.len());
    for sample in 0..n_samples {
        for expanded_name in &standard_expanded_names {
            let col = expanded_names
                .iter()
                .position(|candidate| candidate == expanded_name)
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula standard column '{expanded_name}' is not present in expanded design"
                    ),
                })?;
            standard_values.push(
                expanded_design
                    .matrix()
                    .get(sample, col)
                    .copied()
                    .unwrap_or(0.0),
            );
        }
    }
    let standard_design = DesignMatrix::from_row_major(
        n_samples,
        standard_names.len(),
        standard_values,
        Some(standard_names),
    )?;

    let mut coefficient_groups = Vec::with_capacity(standard_design.n_coefficients());
    for expanded_name in &standard_expanded_names {
        let col = expanded_names
            .iter()
            .position(|candidate| candidate == expanded_name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula standard column '{expanded_name}' is not present in expanded design"
                ),
            })?;
        coefficient_groups.push(vec![col]);
    }

    Ok(ExpandedAdditiveFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        factor_levels: state
            .selected_factors
            .iter()
            .map(|factor| ordered_levels(factor.sample_levels, factor.reference))
            .collect(),
        numeric_covariates: state
            .selected_numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string())
            .collect(),
        interactions: state
            .factor_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.left_factor, interaction.right_factor))
            .collect(),
        factor_numeric_interactions: state
            .factor_numeric_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.factor, interaction.numeric))
            .collect(),
        numeric_interactions: state
            .numeric_interactions
            .iter()
            .map(|interaction| {
                format!("{}:{}", interaction.left_numeric, interaction.right_numeric)
            })
            .collect(),
        higher_order_interactions: state
            .higher_order_interactions
            .iter()
            .map(formula_higher_order_interaction_name)
            .collect(),
    })
}

fn formula_used_factors<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    factors: &'a [ExpandedFactorSpec<'a>],
) -> Result<Vec<ExpandedFactorSpec<'a>>, DeseqError> {
    let mut used = Vec::new();
    for factor in &state.selected_factors {
        push_unique_factor_spec(&mut used, factor.clone())?;
    }
    for interaction in &state.factor_interactions {
        push_unique_factor_spec(
            &mut used,
            formula_factor_spec(interaction.left_factor, factors)?,
        )?;
        push_unique_factor_spec(
            &mut used,
            formula_factor_spec(interaction.right_factor, factors)?,
        )?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_unique_factor_spec(&mut used, formula_factor_spec(interaction.factor, factors)?)?;
    }
    for interaction in &state.higher_order_interactions {
        for variable in &interaction.variables {
            if let FormulaVariableRef::Factor(factor) = variable {
                push_unique_factor_spec(&mut used, formula_factor_spec(factor, factors)?)?;
            }
        }
    }
    Ok(used)
}

fn formula_used_numeric_covariates<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Vec<ExpandedNumericSpec<'a>>, DeseqError> {
    let mut used = Vec::new();
    for covariate in &state.selected_numeric_covariates {
        push_unique_numeric_spec(&mut used, covariate.clone())?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.numeric, numeric_covariates)?,
        )?;
    }
    for interaction in &state.numeric_interactions {
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.left_numeric, numeric_covariates)?,
        )?;
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.right_numeric, numeric_covariates)?,
        )?;
    }
    for interaction in &state.higher_order_interactions {
        for variable in &interaction.variables {
            if let FormulaVariableRef::Numeric(numeric) = variable {
                push_unique_numeric_spec(
                    &mut used,
                    formula_numeric_spec(numeric, numeric_covariates)?,
                )?;
            }
        }
    }
    Ok(used)
}

fn push_unique_factor_spec<'a>(
    used: &mut Vec<ExpandedFactorSpec<'a>>,
    factor: ExpandedFactorSpec<'a>,
) -> Result<(), DeseqError> {
    if used
        .iter()
        .any(|candidate| candidate.factor == factor.factor)
    {
        return Ok(());
    }
    used.push(factor);
    Ok(())
}

fn push_unique_numeric_spec<'a>(
    used: &mut Vec<ExpandedNumericSpec<'a>>,
    covariate: ExpandedNumericSpec<'a>,
) -> Result<(), DeseqError> {
    if used
        .iter()
        .any(|candidate| candidate.name == covariate.name)
    {
        return Ok(());
    }
    used.push(covariate);
    Ok(())
}

fn formula_design_sample_count(
    state: &ExpandedFormulaDesignState<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    used_numeric_covariates: &[ExpandedNumericSpec<'_>],
    supplied_factors: &[ExpandedFactorSpec<'_>],
    supplied_numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<usize, DeseqError> {
    if let Some(factor) = used_factors.first() {
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = used_numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "numeric covariate values must be non-empty".to_string(),
            });
        }
        return Ok(covariate.values.len());
    }
    if !state.has_intercept {
        return Err(DeseqError::InvalidOptions {
            reason: "formula without an intercept must include at least one design term"
                .to_string(),
        });
    }
    if let Some(factor) = supplied_factors.first() {
        if factor.sample_levels.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "intercept-only formula requires non-empty sample metadata".to_string(),
            });
        }
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = supplied_numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "intercept-only formula requires non-empty sample metadata".to_string(),
            });
        }
        return Ok(covariate.values.len());
    }
    Err(DeseqError::InvalidOptions {
        reason: "intercept-only formula requires supplied sample metadata".to_string(),
    })
}

fn validate_formula_used_factor_specs(
    factors: &[ExpandedFactorSpec<'_>],
    n_samples: usize,
) -> Result<(), DeseqError> {
    for (idx, factor) in factors.iter().enumerate() {
        validate_factor_design_inputs(factor.factor, factor.sample_levels, factor.reference)?;
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "formula factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        if factors[..idx]
            .iter()
            .any(|previous| previous.factor == factor.factor)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{}' appears more than once", factor.factor),
            });
        }
    }
    Ok(())
}

fn formula_factor_spec<'a>(
    factor: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
) -> Result<ExpandedFactorSpec<'a>, DeseqError> {
    factors
        .iter()
        .find(|candidate| candidate.factor == factor)
        .cloned()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{factor}' is not present in supplied design metadata"),
        })
}

fn formula_numeric_spec<'a>(
    numeric: &str,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedNumericSpec<'a>, DeseqError> {
    numeric_covariates
        .iter()
        .find(|candidate| candidate.name == numeric)
        .cloned()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric}' is not present in supplied design metadata"
            ),
        })
}

fn formula_factor_levels<'a>(
    factor: &str,
    factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
) -> Result<&'a [String], DeseqError> {
    let idx = factors
        .iter()
        .position(|candidate| candidate.factor == factor)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{factor}' is not present in supplied design metadata"),
        })?;
    Ok(&factor_levels[idx])
}

fn formula_has_factor_main_effect(state: &ExpandedFormulaDesignState<'_>, factor: &str) -> bool {
    state
        .selected_factors
        .iter()
        .any(|candidate| candidate.factor == factor)
}

fn formula_has_numeric_main_effect(state: &ExpandedFormulaDesignState<'_>, numeric: &str) -> bool {
    state
        .selected_numeric_covariates
        .iter()
        .any(|candidate| candidate.name == numeric)
}

fn push_formula_factor_interaction_standard_names(
    interaction: &ExpandedFactorInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let left = formula_factor_spec(interaction.left_factor, used_factors)?;
    let right = formula_factor_spec(interaction.right_factor, used_factors)?;
    let left_levels = formula_factor_levels(interaction.left_factor, used_factors, factor_levels)?;
    let right_levels =
        formula_factor_levels(interaction.right_factor, used_factors, factor_levels)?;
    let left_main = formula_has_factor_main_effect(state, interaction.left_factor);
    let right_main = formula_has_factor_main_effect(state, interaction.right_factor);
    let left_treatment = right_main;
    let right_treatment = left_main;
    for left_level in left_levels
        .iter()
        .filter(|level| !left_treatment || level.as_str() != left.reference)
    {
        for right_level in right_levels
            .iter()
            .filter(|level| !right_treatment || level.as_str() != right.reference)
        {
            standard_names.push(formula_factor_interaction_standard_name(
                FormulaFactorInteractionPiece {
                    factor: interaction.left_factor,
                    level: left_level,
                    reference: left.reference,
                    use_treatment_name: left_treatment,
                },
                FormulaFactorInteractionPiece {
                    factor: interaction.right_factor,
                    level: right_level,
                    reference: right.reference,
                    use_treatment_name: right_treatment,
                },
            ));
            standard_expanded_names.push(interaction_coefficient_name(
                interaction.left_factor,
                left_level,
                interaction.right_factor,
                right_level,
            ));
        }
    }
    Ok(())
}

fn push_formula_factor_numeric_standard_names(
    interaction: &ExpandedFactorNumericInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let factor = formula_factor_spec(interaction.factor, used_factors)?;
    let levels = formula_factor_levels(interaction.factor, used_factors, factor_levels)?;
    let numeric_main = formula_has_numeric_main_effect(state, interaction.numeric);
    let use_treatment_levels = numeric_main;
    for level in levels
        .iter()
        .filter(|level| !use_treatment_levels || level.as_str() != factor.reference)
    {
        if use_treatment_levels {
            standard_names.push(standard_factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                factor.reference,
                interaction.numeric,
            ));
        } else {
            standard_names.push(factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                interaction.numeric,
            ));
        }
        standard_expanded_names.push(factor_numeric_interaction_coefficient_name(
            interaction.factor,
            level,
            interaction.numeric,
        ));
    }
    Ok(())
}

fn formula_higher_order_expanded_names(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
) -> Result<Vec<String>, DeseqError> {
    formula_higher_order_names(interaction, used_factors, factor_levels, &|_| false)
}

fn formula_higher_order_expanded_values(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    used_numeric_covariates: &[ExpandedNumericSpec<'_>],
    factor_levels: &[Vec<String>],
    sample: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut values = Vec::new();
    let level_grid = formula_higher_order_level_grid(interaction, used_factors, factor_levels)?;
    for levels in level_grid {
        let mut value = 1.0;
        for (variable, level) in interaction.variables.iter().zip(levels) {
            match variable {
                FormulaVariableRef::Factor(factor) => {
                    let factor = formula_factor_spec(factor, used_factors)?;
                    let sample_level = factor.sample_levels[sample].as_str();
                    let Some(level) = level else {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula factor interaction level is missing".to_string(),
                        });
                    };
                    value *= (sample_level == level) as u8 as f64;
                }
                FormulaVariableRef::Numeric(numeric) => {
                    value *= formula_numeric_spec(numeric, used_numeric_covariates)?.values[sample];
                }
            }
        }
        values.push(value);
    }
    Ok(values)
}

fn push_formula_higher_order_standard_names(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let treatment_filter = |factor: &str| formula_has_factor_main_effect(state, factor);
    let standard =
        formula_higher_order_names(interaction, used_factors, factor_levels, &treatment_filter)?;
    let expanded_by_name = formula_higher_order_names_with_choice_filter(
        interaction,
        used_factors,
        factor_levels,
        &|_| false,
        &treatment_filter,
    )?;

    for (standard_name, expanded_name) in standard.into_iter().zip(expanded_by_name) {
        standard_names.push(standard_name);
        standard_expanded_names.push(expanded_name);
    }
    Ok(())
}

fn formula_higher_order_names<F>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    use_treatment_for_factor: &F,
) -> Result<Vec<String>, DeseqError>
where
    F: Fn(&str) -> bool,
{
    formula_higher_order_names_with_choice_filter(
        interaction,
        used_factors,
        factor_levels,
        use_treatment_for_factor,
        use_treatment_for_factor,
    )
}

fn formula_higher_order_names_with_choice_filter<F, G>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    use_treatment_for_factor: &F,
    use_treatment_choices_for_factor: &G,
) -> Result<Vec<String>, DeseqError>
where
    F: Fn(&str) -> bool,
    G: Fn(&str) -> bool,
{
    let mut names = Vec::new();
    let level_grid = formula_higher_order_level_grid_with_treatment(
        interaction,
        used_factors,
        factor_levels,
        use_treatment_choices_for_factor,
    )?;
    for levels in level_grid {
        let mut pieces = Vec::with_capacity(interaction.variables.len());
        for (variable, level) in interaction.variables.iter().zip(levels) {
            match variable {
                FormulaVariableRef::Factor(factor) => {
                    let factor_spec = formula_factor_spec(factor, used_factors)?;
                    let Some(level) = level else {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula factor interaction level is missing".to_string(),
                        });
                    };
                    if use_treatment_for_factor(factor) {
                        pieces.push(standard_factor_coefficient_name(
                            factor,
                            level,
                            factor_spec.reference,
                        ));
                    } else {
                        pieces.push(expanded_factor_coefficient_name(factor, level));
                    }
                }
                FormulaVariableRef::Numeric(numeric) => pieces.push((*numeric).to_string()),
            }
        }
        names.push(pieces.join(":"));
    }
    Ok(names)
}

fn formula_higher_order_level_grid<'a>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
) -> Result<Vec<Vec<Option<&'a str>>>, DeseqError> {
    formula_higher_order_level_grid_with_treatment(
        interaction,
        used_factors,
        factor_levels,
        &|_| false,
    )
}

fn formula_higher_order_level_grid_with_treatment<'a, F>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
    use_treatment_for_factor: &F,
) -> Result<Vec<Vec<Option<&'a str>>>, DeseqError>
where
    F: Fn(&str) -> bool,
{
    let mut grids: Vec<Vec<Option<&str>>> = vec![Vec::new()];
    for variable in &interaction.variables {
        let choices = match variable {
            FormulaVariableRef::Factor(factor) => {
                let factor_spec = formula_factor_spec(factor, used_factors)?;
                formula_factor_levels(factor, used_factors, factor_levels)?
                    .iter()
                    .filter(|level| {
                        !use_treatment_for_factor(factor) || level.as_str() != factor_spec.reference
                    })
                    .map(|level| Some(level.as_str()))
                    .collect::<Vec<_>>()
            }
            FormulaVariableRef::Numeric(_) => vec![None],
        };
        let mut next = Vec::with_capacity(grids.len() * choices.len());
        for prefix in &grids {
            for choice in &choices {
                let mut row = prefix.clone();
                row.push(*choice);
                next.push(row);
            }
        }
        grids = next;
    }
    Ok(grids)
}

fn formula_higher_order_interaction_name(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
) -> String {
    interaction
        .variables
        .iter()
        .map(FormulaVariableRef::name)
        .collect::<Vec<_>>()
        .join(":")
}

struct FormulaFactorInteractionPiece<'a> {
    factor: &'a str,
    level: &'a str,
    reference: &'a str,
    use_treatment_name: bool,
}

fn formula_factor_interaction_standard_name(
    left: FormulaFactorInteractionPiece<'_>,
    right: FormulaFactorInteractionPiece<'_>,
) -> String {
    let left = if left.use_treatment_name {
        standard_factor_coefficient_name(left.factor, left.level, left.reference)
    } else {
        expanded_factor_coefficient_name(left.factor, left.level)
    };
    let right = if right.use_treatment_name {
        standard_factor_coefficient_name(right.factor, right.level, right.reference)
    } else {
        expanded_factor_coefficient_name(right.factor, right.level)
    };
    format!("{left}:{right}")
}

fn add_factor_numeric_formula_interaction<'a>(
    factor: &'a str,
    numeric: &'a str,
    factor_numeric_interactions: &mut Vec<ExpandedFactorNumericInteractionSpec<'a>>,
    term: &str,
) -> Result<(), DeseqError> {
    if factor_numeric_interactions
        .iter()
        .any(|interaction| interaction.factor == factor && interaction.numeric == numeric)
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{term}' appears more than once"),
        });
    }
    factor_numeric_interactions.push(ExpandedFactorNumericInteractionSpec { factor, numeric });
    Ok(())
}

fn split_formula_pieces(term: &str, delimiter: char) -> Result<Vec<&str>, DeseqError> {
    let pieces = term.split(delimiter).map(str::trim).collect::<Vec<_>>();
    if pieces.len() < 2 || pieces.iter().any(|piece| piece.is_empty()) {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' must contain non-empty variables"),
        });
    }
    Ok(pieces)
}

fn formula_piece_combinations<'a>(pieces: &[&'a str], order: usize) -> Vec<Vec<&'a str>> {
    fn push_combinations<'a>(
        pieces: &[&'a str],
        order: usize,
        start: usize,
        current: &mut Vec<&'a str>,
        output: &mut Vec<Vec<&'a str>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=pieces.len() - remaining {
            current.push(pieces[idx]);
            push_combinations(pieces, order, idx + 1, current, output);
            current.pop();
        }
    }

    if order == 0 || order > pieces.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_combinations(pieces, order, 0, &mut Vec::new(), &mut output);
    output
}

fn validate_formula_variable(variable: &str) -> Result<(), DeseqError> {
    if variable.is_empty()
        || variable
            .chars()
            .any(|character| character.is_whitespace() || "+-*/:^()".contains(character))
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula variable '{variable}' is not a supported bare variable name"),
        });
    }
    Ok(())
}

fn validate_additive_factor_specs(factors: &[ExpandedFactorSpec<'_>]) -> Result<(), DeseqError> {
    if factors.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "at least one factor is required".to_string(),
        });
    }
    let n_samples = factors[0].sample_levels.len();
    for (idx, factor) in factors.iter().enumerate() {
        validate_factor_design_inputs(factor.factor, factor.sample_levels, factor.reference)?;
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "additive factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        if factors[..idx]
            .iter()
            .any(|previous| previous.factor == factor.factor)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{}' appears more than once", factor.factor),
            });
        }
    }
    Ok(())
}

fn additive_design_sample_count(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<usize, DeseqError> {
    if let Some(factor) = factors.first() {
        return Ok(factor.sample_levels.len());
    }
    let Some(covariate) = numeric_covariates.first() else {
        return Err(DeseqError::InvalidOptions {
            reason: "at least one factor or numeric covariate is required".to_string(),
        });
    };
    if covariate.values.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "numeric covariate values must be non-empty".to_string(),
        });
    }
    Ok(covariate.values.len())
}

fn validate_numeric_covariate_specs(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    n_samples: usize,
) -> Result<(), DeseqError> {
    for (idx, covariate) in numeric_covariates.iter().enumerate() {
        if covariate.name.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("numeric covariate {idx} name must be non-empty"),
            });
        }
        if covariate.values.len() != n_samples {
            return Err(invalid_dimensions(
                "additive numeric covariate values",
                n_samples,
                covariate.values.len(),
            ));
        }
        if numeric_covariates[..idx]
            .iter()
            .any(|previous| previous.name == covariate.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric covariate '{}' appears more than once",
                    covariate.name
                ),
            });
        }
        if !covariate.values.iter().all(|value| value.is_finite()) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric covariate '{}' contains non-finite values",
                    covariate.name
                ),
            });
        }
    }
    Ok(())
}

fn resolve_interaction_indices(
    factors: &[ExpandedFactorSpec<'_>],
    interactions: &[ExpandedFactorInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.left_factor.is_empty() || interaction.right_factor.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("interaction {idx} factor names must be non-empty"),
            });
        }
        if interaction.left_factor == interaction.right_factor {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "interaction {idx} cannot use factor '{}' twice",
                    interaction.left_factor
                ),
            });
        }
        let left = factor_index(factors, interaction.left_factor)?;
        let right = factor_index(factors, interaction.right_factor)?;
        let ordered = if left < right {
            (left, right)
        } else {
            (right, left)
        };
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "interaction '{}:{}' appears more than once",
                    interaction.left_factor, interaction.right_factor
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn resolve_factor_numeric_interaction_indices(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedFactorNumericInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.factor.is_empty() || interaction.numeric.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor-numeric interaction {idx} names must be non-empty"),
            });
        }
        let factor_idx = factor_index(factors, interaction.factor)?;
        let numeric_idx = numeric_index(numeric_covariates, interaction.numeric)?;
        let ordered = (factor_idx, numeric_idx);
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "factor-numeric interaction '{}:{}' appears more than once",
                    interaction.factor, interaction.numeric
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn resolve_numeric_interaction_indices(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedNumericInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.left_numeric.is_empty() || interaction.right_numeric.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("numeric interaction {idx} covariate names must be non-empty"),
            });
        }
        if interaction.left_numeric == interaction.right_numeric {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric interaction {idx} cannot use covariate '{}' twice",
                    interaction.left_numeric
                ),
            });
        }
        let left = numeric_index(numeric_covariates, interaction.left_numeric)?;
        let right = numeric_index(numeric_covariates, interaction.right_numeric)?;
        let ordered = if left < right {
            (left, right)
        } else {
            (right, left)
        };
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric interaction '{}:{}' appears more than once",
                    interaction.left_numeric, interaction.right_numeric
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn factor_index(factors: &[ExpandedFactorSpec<'_>], factor: &str) -> Result<usize, DeseqError> {
    factors
        .iter()
        .position(|candidate| candidate.factor == factor)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("interaction factor '{factor}' is not present in additive factors"),
        })
}

fn numeric_index(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    numeric: &str,
) -> Result<usize, DeseqError> {
    numeric_covariates
        .iter()
        .position(|candidate| candidate.name == numeric)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "interaction numeric covariate '{numeric}' is not present in numeric covariates"
            ),
        })
}

fn validate_unique_coefficient_names(names: &[String], context: &str) -> Result<(), DeseqError> {
    for (idx, name) in names.iter().enumerate() {
        if names[..idx].iter().any(|previous| previous == name) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} has duplicate coefficient name '{name}'"),
            });
        }
    }
    Ok(())
}

fn validate_factor_design_inputs<S: AsRef<str>>(
    factor: &str,
    sample_levels: &[S],
    reference: &str,
) -> Result<(), DeseqError> {
    if factor.is_empty() || reference.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "factor and reference level must be non-empty".to_string(),
        });
    }
    if sample_levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "sample levels must be non-empty".to_string(),
        });
    }
    let mut has_reference = false;
    for (idx, level) in sample_levels.iter().enumerate() {
        let level = level.as_ref();
        if level.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("sample level {idx} must be non-empty"),
            });
        }
        has_reference |= level == reference;
    }
    if !has_reference {
        return Err(DeseqError::InvalidOptions {
            reason: format!("reference level '{reference}' is not present in sample levels"),
        });
    }
    Ok(())
}

fn ordered_levels<S: AsRef<str>>(sample_levels: &[S], reference: &str) -> Vec<String> {
    let mut levels = vec![reference.to_string()];
    for level in sample_levels {
        let level = level.as_ref();
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.to_string());
        }
    }
    levels
}

fn expanded_factor_coefficient_name(factor: &str, level: &str) -> String {
    format!("{factor}{level}")
}

fn standard_factor_coefficient_name(factor: &str, level: &str, reference: &str) -> String {
    format!("{factor}_{level}_vs_{reference}")
}

fn interaction_coefficient_name(
    left_factor: &str,
    left_level: &str,
    right_factor: &str,
    right_level: &str,
) -> String {
    format!("{left_factor}{left_level}:{right_factor}{right_level}")
}

fn standard_interaction_coefficient_name(
    left_factor: &str,
    left_level: &str,
    left_reference: &str,
    right_factor: &str,
    right_level: &str,
    right_reference: &str,
) -> String {
    format!(
        "{left_factor}_{left_level}_vs_{left_reference}:{right_factor}_{right_level}_vs_{right_reference}"
    )
}

fn factor_numeric_interaction_coefficient_name(factor: &str, level: &str, numeric: &str) -> String {
    format!("{factor}{level}:{numeric}")
}

fn standard_factor_numeric_interaction_coefficient_name(
    factor: &str,
    level: &str,
    reference: &str,
    numeric: &str,
) -> String {
    format!("{factor}_{level}_vs_{reference}:{numeric}")
}

fn numeric_interaction_coefficient_name(left_numeric: &str, right_numeric: &str) -> String {
    format!("{left_numeric}:{right_numeric}")
}
