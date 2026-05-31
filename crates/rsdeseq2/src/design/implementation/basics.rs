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
