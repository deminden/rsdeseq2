pub(super) use super::super::DesignOrigin;
use super::*;
pub(super) use crate::errors::{Error, invalid_dimensions};
use crate::errors::{ErrorKind, missing_prerequisite};
pub(super) use crate::math::rank::{R_QR_TOLERANCE, matrix_rank, matrix_rank_r_qr};
pub(super) use crate::matrix::RowMajorMatrix;
use rcompat_numeric::matrix::try_alloc_matrix;
use std::sync::Arc;

/// A design matrix generated outside the Rust core or by the Rust design helpers.
///
/// Rows are samples and columns are coefficients. The Rust formula helper covers
/// a primitive DESeq2-style subset; callers can still provide arbitrary matrices
/// directly when they need richer formula semantics. Equality compares numeric
/// values and coefficient names; compare [`Self::origin`] separately when
/// formula-dependent result policy is relevant.
#[derive(Clone)]
pub struct DesignMatrix {
    values: Arc<DesignValues>,
    origin: DesignOrigin,
    // Immutable design products shared across row optimizers and clones.
    cox_products: std::sync::OnceLock<Option<std::sync::Arc<SharedDesignProducts>>>,
}

struct DesignValues {
    matrix: RowMajorMatrix<f64>,
    coefficient_names: Option<Vec<String>>,
}

#[doc(hidden)]
pub struct SharedDesignProducts {
    #[doc(hidden)]
    pub offsets: Vec<usize>,
    #[doc(hidden)]
    pub products: Vec<(usize, usize, f64)>,
}

impl PartialEq for DesignMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.values.matrix == other.values.matrix
            && self.values.coefficient_names == other.values.coefficient_names
    }
}

impl std::fmt::Debug for DesignMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DesignMatrix")
            .field("matrix", &self.values.matrix)
            .field("coefficient_names", &self.values.coefficient_names)
            .finish()
    }
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
    /// Optional R factor level order. When supplied, design columns follow this
    /// order instead of first-observed sample order, with `reference` still used
    /// as the treatment base.
    pub levels: Option<&'a [String]>,
}

/// Caller-supplied metadata for one additive numeric covariate.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedNumericSpec<'a> {
    /// Covariate name used as the coefficient name.
    pub name: &'a str,
    /// Per-sample finite numeric values in count-column order.
    pub values: &'a [f64],
}

/// Owned formula model frame for wrapper and object-style callers.
///
/// This bridges user-facing sample metadata to the borrowed primitive formula
/// design helpers. Factor references can be supplied explicitly; when omitted,
/// the first declared level is used as the treatment reference when `levels`
/// are present, otherwise the first observed level is used.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FormulaModelFrame {
    /// Categorical columns available to formula parsing.
    pub factors: Vec<FormulaFactorColumn>,
    /// Numeric columns available to formula parsing and supported transforms.
    pub numeric_covariates: Vec<FormulaNumericColumn>,
}

/// Resolved factor reference metadata for one model-frame factor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolvedFormulaFactorReference<'a> {
    /// Factor column name.
    pub factor: &'a str,
    /// Treatment reference/base level used by formula design construction.
    pub reference: &'a str,
    /// Optional declared R factor level order, if supplied.
    pub levels: Option<&'a [String]>,
}

/// Resolved numeric covariate metadata for one model-frame column.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedFormulaNumericCovariate<'a> {
    /// Numeric column name.
    pub name: &'a str,
    /// Per-sample finite numeric values.
    pub values: &'a [f64],
}

/// Owned categorical model-frame column.
#[derive(Clone, Debug, PartialEq)]
pub struct FormulaFactorColumn {
    /// Column name used in formula terms.
    pub name: String,
    /// Per-sample factor levels in count-column order.
    pub sample_levels: Vec<String>,
    /// Optional R factor level order. When supplied, every observed sample
    /// level must be present here and model-matrix columns follow this order.
    pub levels: Option<Vec<String>>,
    /// Optional treatment reference/base level. Defaults to the first declared
    /// level when [`Self::levels`] is present, otherwise the first observed
    /// sample level.
    pub reference: Option<String>,
}

/// Owned numeric model-frame column.
#[derive(Clone, Debug, PartialEq)]
pub struct FormulaNumericColumn {
    /// Column name used in formula terms.
    pub name: String,
    /// Per-sample finite numeric values in count-column order.
    pub values: Vec<f64>,
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
    /// Parsed expanded and reported design matrices.
    pub design: ExpandedAdditiveFactorDesign,
    /// Per-sample offset values from supported `offset(numeric)` terms.
    pub offsets: Vec<f64>,
    /// Formula-local model frame after supported transforms have been
    /// materialized.
    ///
    /// This includes the original supplied columns plus derived factor and
    /// numeric columns such as `relevel(condition, ref = "B")`,
    /// `factor(condition, levels = c(...))`, `droplevels(condition)`, and
    /// supported numeric transforms. Downstream result routing can use this to
    /// resolve character contrasts against the same labels and references that
    /// produced the design matrix.
    pub model_frame: FormulaModelFrame,
}

impl FormulaModelFrame {
    /// Validate this model frame using the same rules as formula design construction.
    pub fn validate(&self) -> Result<(), Error> {
        validate_formula_model_frame(self)
    }

    /// Number of samples represented by this model frame.
    pub fn n_samples(&self) -> Result<usize, Error> {
        formula_model_frame_sample_count(self)
    }

    /// Resolved treatment reference for an exact factor column name.
    ///
    /// The resolution order matches formula design construction: explicit
    /// `reference`, then the first declared level, then the first observed
    /// sample level. `Ok(None)` means the model frame is valid but has no
    /// factor with that exact name.
    pub fn resolved_factor_reference(&self, factor_name: &str) -> Result<Option<&str>, Error> {
        self.validate()?;
        self.factors
            .iter()
            .find(|factor| factor.name == factor_name)
            .map(formula_model_frame_factor_reference)
            .transpose()
    }

    /// Resolved treatment reference for a factor column name or R-cleaned alias.
    ///
    /// Exact factor names win first. If no exact factor exists, this accepts
    /// the same R-style cleanup used by formula coefficient aliases, so wrapper
    /// callers can resolve metadata from user-facing names that came through
    /// `make.names`-style paths. Ambiguous cleaned aliases return an error.
    pub fn resolved_factor_reference_by_alias(
        &self,
        factor_name: &str,
    ) -> Result<Option<&str>, Error> {
        self.validate()?;
        if let Some(reference) = self.resolved_factor_reference(factor_name)? {
            return Ok(Some(reference));
        }
        let matches = self
            .factors
            .iter()
            .filter(|factor| {
                r_like_name_candidates(&factor.name)
                    .into_iter()
                    .any(|candidate| candidate == factor_name)
            })
            .collect::<Vec<_>>();
        match matches.as_slice() {
            [factor] => Ok(Some(formula_model_frame_factor_reference(factor)?)),
            [] => Ok(None),
            _ => Err(Error::InvalidOptions {
                reason: format!(
                    "factor '{factor_name}' resolves ambiguously after R-style cleanup"
                ),
            }),
        }
    }

    /// Resolved reference metadata for all factor columns in model-frame order.
    pub fn resolved_factor_references(
        &self,
    ) -> Result<Vec<ResolvedFormulaFactorReference<'_>>, Error> {
        self.validate()?;
        self.factors
            .iter()
            .map(|factor| {
                Ok(ResolvedFormulaFactorReference {
                    factor: factor.name.as_str(),
                    reference: formula_model_frame_factor_reference(factor)?,
                    levels: factor.levels.as_deref(),
                })
            })
            .collect()
    }

    /// Resolved numeric covariate for an exact numeric column name.
    ///
    /// `Ok(None)` means the model frame is valid but has no numeric covariate
    /// with that exact name.
    pub fn resolved_numeric_covariate(
        &self,
        numeric_name: &str,
    ) -> Result<Option<ResolvedFormulaNumericCovariate<'_>>, Error> {
        self.validate()?;
        Ok(self
            .numeric_covariates
            .iter()
            .find(|covariate| covariate.name == numeric_name)
            .map(|covariate| ResolvedFormulaNumericCovariate {
                name: covariate.name.as_str(),
                values: covariate.values.as_slice(),
            }))
    }

    /// Resolved numeric covariate for a numeric column name or R-cleaned alias.
    ///
    /// Exact numeric names win first. If no exact numeric column exists, this
    /// accepts the same R-style cleanup used by formula numeric transform
    /// lookup. Ambiguous cleaned aliases return an error.
    pub fn resolved_numeric_covariate_by_alias(
        &self,
        numeric_name: &str,
    ) -> Result<Option<ResolvedFormulaNumericCovariate<'_>>, Error> {
        self.validate()?;
        if let Some(covariate) = self.resolved_numeric_covariate(numeric_name)? {
            return Ok(Some(covariate));
        }
        let matches = self
            .numeric_covariates
            .iter()
            .filter(|covariate| {
                r_like_name_candidates(&covariate.name)
                    .into_iter()
                    .any(|candidate| candidate == numeric_name)
            })
            .collect::<Vec<_>>();
        match matches.as_slice() {
            [covariate] => Ok(Some(ResolvedFormulaNumericCovariate {
                name: covariate.name.as_str(),
                values: covariate.values.as_slice(),
            })),
            [] => Ok(None),
            _ => Err(Error::InvalidOptions {
                reason: format!(
                    "numeric covariate '{numeric_name}' resolves ambiguously after R-style cleanup"
                ),
            }),
        }
    }

    /// Resolved metadata for all numeric columns in model-frame order.
    pub fn resolved_numeric_covariates(
        &self,
    ) -> Result<Vec<ResolvedFormulaNumericCovariate<'_>>, Error> {
        self.validate()?;
        Ok(self
            .numeric_covariates
            .iter()
            .map(|covariate| ResolvedFormulaNumericCovariate {
                name: covariate.name.as_str(),
                values: covariate.values.as_slice(),
            })
            .collect())
    }
}

impl DesignMatrix {
    #[doc(hidden)]
    pub fn cached_cox_products(
        &self,
        build: impl FnOnce() -> Option<SharedDesignProducts>,
    ) -> Option<std::sync::Arc<SharedDesignProducts>> {
        self.cox_products
            .get_or_init(|| build().map(std::sync::Arc::new))
            .clone()
    }

    /// Create an intercept-only design matrix with one all-ones coefficient.
    pub fn intercept_only(n_samples: usize) -> Result<Self, Error> {
        let mut values = try_alloc_matrix(n_samples, 1, "intercept-only design")?;
        values.resize(n_samples, 1.0);
        Self::from_row_major(n_samples, 1, values, Some(vec!["Intercept".to_string()]))
    }

    /// Create a design matrix from row-major values.
    pub fn from_row_major(
        n_samples: usize,
        n_coefficients: usize,
        values: Vec<f64>,
        coefficient_names: Option<Vec<String>>,
    ) -> Result<Self, Error> {
        if let Some(names) = &coefficient_names
            && names.len() != n_coefficients
        {
            return Err(invalid_dimensions(
                "design coefficient names",
                n_coefficients,
                names.len(),
            ));
        }
        let matrix = RowMajorMatrix::from_row_major(n_samples, n_coefficients, values)?;
        matrix.validate_finite("design matrix")?;
        Ok(Self {
            values: Arc::new(DesignValues {
                matrix,
                coefficient_names,
            }),
            origin: DesignOrigin::matrix(),
            cox_products: std::sync::OnceLock::new(),
        })
    }

    /// Attach explicit provenance without altering numerical design equality,
    /// coefficients, rank checks or cached numeric products.
    pub fn with_origin(mut self, origin: DesignOrigin) -> Self {
        self.origin = origin;
        self
    }

    /// Formula/matrix provenance used by formula-specific result policy.
    pub fn origin(&self) -> &DesignOrigin {
        &self.origin
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.values.matrix.n_rows()
    }

    /// Number of coefficients.
    pub fn n_coefficients(&self) -> usize {
        self.values.matrix.n_cols()
    }

    /// Reusable sample-index span.
    pub fn sample_indices(&self) -> core::range::Range<usize> {
        self.values.matrix.row_indices()
    }

    /// Reusable coefficient-index span.
    pub fn coefficient_indices(&self) -> core::range::Range<usize> {
        self.values.matrix.col_indices()
    }

    /// Return a contiguous sample-row block in row-major order.
    ///
    /// The range accepts both legacy range syntax (`1..3`) and the newer
    /// `core::range` types.
    pub fn sample_rows<R: core::ops::RangeBounds<usize>>(
        &self,
        samples: R,
    ) -> Result<&[f64], Error> {
        self.values.matrix.rows(samples)
    }

    /// Matrix values.
    pub fn matrix(&self) -> &RowMajorMatrix<f64> {
        &self.values.matrix
    }

    /// Optional coefficient names.
    pub fn coefficient_names(&self) -> Option<&[String]> {
        self.values.coefficient_names.as_deref()
    }

    /// Resolve a coefficient name to its zero-based column index.
    pub fn coefficient_index(&self, name: &str) -> Result<usize, Error> {
        let names = self.coefficient_names().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "coefficient name resolution",
                "coefficient names",
            )
        })?;
        names
            .iter()
            .position(|candidate| candidate == name)
            .ok_or_else(|| Error::CoefficientNotFound {
                reason: format!("coefficient '{name}' is not present in coefficient names"),
                candidates: vec![name.to_string()],
            })
    }

    /// Numerical rank as R's default `qr(x)$rank` computes it (LINPACK
    /// `dqrdc2`, tolerance 1e-7), the rank DESeq2's `checkFullRank` uses.
    pub fn rank(&self) -> Result<usize, Error> {
        matrix_rank_r_qr(
            self.values.matrix.as_slice(),
            self.n_samples(),
            self.n_coefficients(),
            R_QR_TOLERANCE,
        )
    }

    /// Numerical rank using a caller-supplied absolute tolerance.
    pub fn rank_with_tolerance(&self, tolerance: f64) -> Result<usize, Error> {
        matrix_rank(
            self.values.matrix.as_slice(),
            self.n_samples(),
            self.n_coefficients(),
            tolerance,
        )
    }

    /// Whether the design has full column rank using the default tolerance.
    pub fn is_full_rank(&self) -> Result<bool, Error> {
        Ok(self.rank()? == self.n_coefficients())
    }

    /// Whether the design has full column rank using a caller-supplied tolerance.
    pub fn is_full_rank_with_tolerance(&self, tolerance: f64) -> Result<bool, Error> {
        Ok(self.rank_with_tolerance(tolerance)? == self.n_coefficients())
    }

    /// Return an error if the design is not full column rank.
    pub fn validate_full_rank(&self, context: &str) -> Result<(), Error> {
        self.validate_rank_result(context, self.rank()?, 0.0)
    }

    /// Return an error if the design is not full column rank.
    pub fn validate_full_rank_with_tolerance(
        &self,
        context: &str,
        tolerance: f64,
    ) -> Result<(), Error> {
        self.validate_rank_result(context, self.rank_with_tolerance(tolerance)?, tolerance)
    }

    fn validate_rank_result(
        &self,
        context: &str,
        rank: usize,
        tolerance: f64,
    ) -> Result<(), Error> {
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
        Err(Error::InvalidOptions { reason })
    }

    fn zero_columns(&self, tolerance: f64) -> Vec<String> {
        let mut columns = Vec::new();
        for col in 0..self.n_coefficients() {
            let is_zero = (0..self.n_samples()).all(|row| {
                self.values
                    .matrix
                    .get(row, col)
                    .copied()
                    .unwrap_or(0.0)
                    .abs()
                    <= tolerance
            });
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
