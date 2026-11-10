//! Owned frozen dispersion curves. Imported locfit knots are evaluated without
//! refitting, foreign pointers, or an R interpreter.
use super::{validate_normalized_counts_and_max, vst_local_from_dispersion_grid};
use crate::design::expression::{self, Data, EvaluationOptions, Expression, Value};
use crate::{
    errors::DeseqError,
    matrix::{CountView, MatrixView, RowMajorMatrix},
};

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidDispersion {
        reason: reason.into(),
    }
}

fn grid_from_max(max: f64) -> Vec<f64> {
    let end = max.asinh();
    let delta = end / 999.0;
    (1..1000)
        .map(|i| (if i == 999 { end } else { i as f64 * delta }).sinh())
        .collect()
}

/// Build the frozen VST grid directly from borrowed counts and normalization
/// inputs. No full normalized matrix or gene-by-sample transport is allocated.
/// Normalization factors take precedence; every supplied factor is validated.
pub fn frozen_local_vst_grid_from_counts(
    counts: CountView<'_>,
    size_factors: Option<&[f64]>,
    factors: Option<MatrixView<'_, f64>>,
) -> Result<Vec<f64>, DeseqError> {
    if let Some(sf) = size_factors {
        crate::normalization::validate_size_factors(sf, counts.n_samples())?;
    }
    if let Some(nf) = factors {
        if nf.n_rows() != counts.n_genes() || nf.n_cols() != counts.n_samples() {
            return Err(invalid(
                "frozen local normalization factor shape differs from counts",
            ));
        }
    } else if size_factors.is_none() {
        return Err(invalid(
            "frozen local grid requires size or normalization factors",
        ));
    }
    let mut max = 0.0_f64;
    for sample in 0..counts.n_samples() {
        for gene in 0..counts.n_genes() {
            let factor = if let Some(nf) = factors {
                *nf.at(gene, sample)
            } else {
                size_factors.expect("validated size factors")[sample]
            };
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid(
                    "frozen local normalization factors must be finite and positive",
                ));
            }
            let normalized = f64::from(counts.at(gene, sample)) / factor;
            if !normalized.is_finite() {
                return Err(invalid("frozen local normalized counts must be finite"));
            }
            max = max.max(normalized);
        }
    }
    Ok(grid_from_max(max))
}

/// Saved one-dimensional quadratic locfit predictor. The local residual values
/// and first derivatives use cubic Hermite interpolation; outside the knot
/// interval they use the boundary tangent. A centered global quadratic is added
/// afterward. These are fitted coefficients, not training data for a refit.
#[derive(Clone, Debug, PartialEq)]
pub struct SavedLocalDispersionFit {
    knots: Vec<f64>,
    values: Vec<f64>,
    slopes: Vec<f64>,
    polynomial: [f64; 4],
}
impl SavedLocalDispersionFit {
    /// Import unordered log-mean knots, aligned residual values/slopes and
    /// `[center, intercept, slope, second_derivative]` global coefficients.
    pub fn new(
        knots: Vec<f64>,
        values: Vec<f64>,
        slopes: Vec<f64>,
        polynomial: [f64; 4],
    ) -> Result<Self, DeseqError> {
        if knots.len() < 2
            || values.len() != knots.len()
            || slopes.len() != knots.len()
            || knots
                .iter()
                .chain(&values)
                .chain(&slopes)
                .chain(&polynomial)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "saved local fit requires aligned finite knots, values, slopes and quadratic coefficients",
            ));
        }
        let mut order = (0..knots.len()).collect::<Vec<_>>();
        order.sort_by(|a, b| knots[*a].total_cmp(&knots[*b]));
        let reorder = |x: &[f64]| order.iter().map(|i| x[*i]).collect::<Vec<_>>();
        let knots = reorder(&knots);
        if knots.windows(2).any(|x| x[0] >= x[1]) {
            return Err(invalid("saved local fit knots must be distinct"));
        }
        Ok(Self {
            knots,
            values: reorder(&values),
            slopes: reorder(&slopes),
            polynomial,
        })
    }

    /// Predict positive finite dispersions at positive finite means, including
    /// outside the saved knot interval.
    pub fn evaluate(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        means
            .iter()
            .map(|mean| {
                if !mean.is_finite() || *mean <= 0.0 {
                    return Err(invalid("saved local fit means must be finite and positive"));
                }
                let x = mean.ln();
                let upper = self.knots.partition_point(|k| *k <= x);
                let residual = if upper == 0 || upper == self.knots.len() {
                    let i = if upper == 0 { 0 } else { upper - 1 };
                    self.values[i] + (x - self.knots[i]) * self.slopes[i]
                } else {
                    let i = upper - 1;
                    let width = self.knots[i + 1] - self.knots[i];
                    let t = (x - self.knots[i]) / width;
                    let u = 1.0 - t;
                    (1.0 + 2.0 * t) * u * u * self.values[i]
                        + t * t * (3.0 - 2.0 * t) * self.values[i + 1]
                        + t * u * u * width * self.slopes[i]
                        + t * t * (t - 1.0) * width * self.slopes[i + 1]
                };
                let dx = x - self.polynomial[0];
                let global = self.polynomial[1]
                    + dx * self.polynomial[2]
                    + (dx * dx / 2.0) * self.polynomial[3];
                let prediction = (residual + global).exp();
                if !prediction.is_finite() || prediction <= 0.0 {
                    return Err(invalid(
                        "saved local fit produced a non-finite or non-positive dispersion",
                    ));
                }
                Ok(prediction)
            })
            .collect()
    }
}

/// A frozen local curve represented only by owned numerical data or a supported
/// built-in expression. No executable frontend function is stored or called.
#[derive(Clone, Debug, PartialEq)]
pub enum FrozenLocalDispersion {
    /// Original saved coefficients, imported without refitting.
    Saved(SavedLocalDispersionFit),
    /// Restricted expression evaluated by the shared standalone expression engine.
    Expression {
        /// Name of the single mean argument.
        argument: String,
        /// Translated expression body.
        expression: Expression,
        /// Captured primitive constants, without promises or executable bindings.
        constants: Vec<(String, Value)>,
    },
}
impl FrozenLocalDispersion {
    /// Evaluate the frozen dispersion function without executing foreign code.
    pub fn evaluate(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        let values = match self {
            Self::Saved(fit) => fit.evaluate(means)?,
            Self::Expression {
                argument,
                expression,
                constants,
            } => {
                if argument.is_empty() || constants.iter().any(|(name, _)| name == argument) {
                    return Err(invalid(
                        "frozen dispersion argument must be nonempty and distinct from constants",
                    ));
                }
                let mut bindings = constants.clone();
                bindings.push((argument.clone(), Value::numbers(means.to_vec())));
                let result =
                    expression::evaluate(expression, &bindings, &EvaluationOptions::default())?;
                let Data::Number { values, .. } = result.data else {
                    return Err(invalid(
                        "frozen dispersion expression must return numeric values",
                    ));
                };
                if result.dimensions.is_some() {
                    return Err(invalid("frozen dispersion expression must return a vector"));
                }
                values
            }
        };
        if values.len() != means.len() || values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid(
                "frozen local dispersion function must return one finite positive value per mean",
            ));
        }
        Ok(values)
    }

    /// Construct and evaluate the complete DESeq2 sinh grid in Rust. Counts are
    /// already normalized by the shared normalization kernel.
    pub fn grid(
        &self,
        normalized: &RowMajorMatrix<f64>,
    ) -> Result<(Vec<f64>, Vec<f64>), DeseqError> {
        let max = validate_normalized_counts_and_max(normalized)?;
        let grid = grid_from_max(max);
        let dispersions = self.evaluate(&grid)?;
        Ok((grid, dispersions))
    }
}

/// Apply local VST using a saved predictor or restricted expression, with the
/// grid, prediction, integration and anchoring all in the standalone engine.
pub fn vst_frozen_local(
    normalized: &RowMajorMatrix<f64>,
    curve: &FrozenLocalDispersion,
    inverse_size_factor_mean: f64,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let (grid, dispersions) = curve.grid(normalized)?;
    vst_local_from_dispersion_grid(normalized, &grid, &dispersions, inverse_size_factor_mean)
}
