//! Request and output types for [`DeseqBuilder::run`](super::DeseqBuilder::run).
//!
//! An [`Analysis`] names one point on each independent axis of a DESeq2 run
//! (design input, test, dispersion source, mean model, result selection,
//! Cook's replacement and how far to run). Builder options such as thread
//! count, size-factor, IRLS and filtering settings stay on the builder.

use crate::contrasts::{ContrastSpec, FactorLevelContrast, ResultsContrast};
use crate::cooks::CooksReplacementOptions;
use crate::design::{DesignMatrix, FormulaModelFrame};
use crate::errors::DeseqError;
use crate::options::FitType;
use crate::results::DeseqResults;

use super::{CooksReplacementTestOutput, CountMatrix, DeseqFit};

/// Where the full design comes from.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum DesignInput<'a> {
    /// A design matrix.
    Matrix(&'a DesignMatrix),
    /// A formula evaluated against the builder's stored model frame
    /// ([`DeseqBuilder::model_frame`](super::DeseqBuilder::model_frame)).
    Formula(&'a str),
    /// No design. Valid only with [`Stage::SizeFactors`].
    None,
}

/// Which test to run.
#[derive(Clone, Copy, Debug, Default)]
#[non_exhaustive]
pub enum TestRequest<'a> {
    /// Wald test (DESeq2 `test = "Wald"`).
    #[default]
    Wald,
    /// Likelihood-ratio test against a reduced design. A formula reduced
    /// design requires a formula full design.
    Lrt {
        /// The reduced design.
        reduced: DesignInput<'a>,
    },
}

/// Where the dispersions come from.
#[derive(Clone, Copy, Debug, Default)]
#[non_exhaustive]
pub enum DispersionRequest<'a> {
    /// Estimate gene-wise, trend and MAP dispersions.
    #[default]
    Estimate,
    /// Use one caller-supplied dispersion per gene.
    Fixed(&'a [f64]),
}

/// Mean model used while estimating dispersions.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[non_exhaustive]
pub enum MeanModel {
    /// DESeq2's automatic choice between linear and GLM means.
    #[default]
    Auto,
    /// GLM fitted means, unless the builder already selects means
    /// automatically ([`DeseqBuilder::automatic_dispersion_means`](super::DeseqBuilder::automatic_dispersion_means)).
    GlmMu,
    /// Linear-model fitted means (`linearMu = TRUE`).
    LinearMu,
}

/// Which result column to report.
#[derive(Clone, Copy, Debug, Default)]
#[non_exhaustive]
pub enum ResultColumn<'a> {
    /// DESeq2's default: the last design coefficient.
    #[default]
    Default,
    /// A design coefficient by index.
    Coefficient(usize),
    /// A design coefficient by name (`results(name = ...)`).
    Name(&'a str),
    /// A numeric contrast vector over the design coefficients.
    Numeric(&'a [f64]),
    /// A primitive contrast specification.
    Spec(&'a ContrastSpec),
    /// A factor-level contrast with explicit per-sample levels.
    FactorLevel(FactorLevelContrast<'a>),
    /// A DESeq2 `results(contrast = ...)` request.
    Results {
        /// The contrast as given to `results()`.
        contrast: &'a ResultsContrast,
        /// Per-sample levels of the contrast factor, used for DESeq2's
        /// character-contrast all-zero rule.
        sample_levels: Option<&'a [String]>,
        /// A model frame to resolve character contrasts against, instead of
        /// the builder's stored one. Cannot be combined with `sample_levels`.
        model_frame: Option<&'a FormulaModelFrame>,
    },
}

/// Cook's distance outlier replacement.
#[derive(Clone, Copy, Debug, Default)]
#[non_exhaustive]
pub enum Replacement<'a> {
    /// No replacement refit.
    #[default]
    Off,
    /// Replace outlier counts and refit affected genes (`minReplicatesForReplace`).
    Cooks(&'a CooksReplacementOptions),
}

/// How far to run.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[non_exhaustive]
pub enum Stage {
    /// Size factors and base means.
    SizeFactors,
    /// Gene-wise dispersion estimates.
    GeneWise,
    /// The dispersion trend.
    Trend,
    /// MAP dispersions.
    Map,
    /// The full pipeline with result rows.
    #[default]
    Results,
}

/// One analysis request for [`DeseqBuilder::run`](super::DeseqBuilder::run).
///
/// [`Analysis::new`] matches `fit_with_results`: Wald test, estimated
/// dispersions, automatic means, the last coefficient and no replacement.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct Analysis<'a> {
    /// Raw counts.
    pub counts: &'a CountMatrix,
    /// The full design.
    pub design: DesignInput<'a>,
    /// The test.
    pub test: TestRequest<'a>,
    /// The dispersion source.
    pub dispersion: DispersionRequest<'a>,
    /// The mean model for dispersion estimation.
    pub means: MeanModel,
    /// The reported result column.
    pub result: ResultColumn<'a>,
    /// Cook's replacement.
    pub replacement: Replacement<'a>,
    /// How far to run.
    pub stop_after: Stage,
    /// Dispersion trend type for this run; `None` uses the builder's.
    pub fit_type: Option<FitType>,
}

impl<'a> Analysis<'a> {
    /// A default analysis of `counts` against a design matrix.
    pub fn new(counts: &'a CountMatrix, design: &'a DesignMatrix) -> Self {
        Self::with_design(counts, DesignInput::Matrix(design))
    }

    /// A default analysis of `counts` against a formula and the builder's model frame.
    pub fn formula(counts: &'a CountMatrix, formula: &'a str) -> Self {
        Self::with_design(counts, DesignInput::Formula(formula))
    }

    /// A size-factor-only analysis without a design.
    pub fn without_design(counts: &'a CountMatrix) -> Self {
        Self::with_design(counts, DesignInput::None).stop_after(Stage::SizeFactors)
    }

    fn with_design(counts: &'a CountMatrix, design: DesignInput<'a>) -> Self {
        Self {
            counts,
            design,
            test: TestRequest::Wald,
            dispersion: DispersionRequest::Estimate,
            means: MeanModel::Auto,
            result: ResultColumn::Default,
            replacement: Replacement::Off,
            stop_after: Stage::Results,
            fit_type: None,
        }
    }

    /// Run a likelihood-ratio test against a reduced design matrix.
    pub fn lrt(mut self, reduced: &'a DesignMatrix) -> Self {
        self.test = TestRequest::Lrt {
            reduced: DesignInput::Matrix(reduced),
        };
        self
    }

    /// Run a likelihood-ratio test against a reduced formula.
    pub fn lrt_formula(mut self, reduced: &'a str) -> Self {
        self.test = TestRequest::Lrt {
            reduced: DesignInput::Formula(reduced),
        };
        self
    }

    /// Set the test.
    pub fn test(mut self, test: TestRequest<'a>) -> Self {
        self.test = test;
        self
    }

    /// Use caller-supplied dispersions.
    pub fn fixed_dispersion(mut self, dispersions: &'a [f64]) -> Self {
        self.dispersion = DispersionRequest::Fixed(dispersions);
        self
    }

    /// Set the mean model for dispersion estimation.
    pub fn means(mut self, means: MeanModel) -> Self {
        self.means = means;
        self
    }

    /// Set the reported result column.
    pub fn result(mut self, result: ResultColumn<'a>) -> Self {
        self.result = result;
        self
    }

    /// Replace Cook's outliers and refit.
    pub fn cooks_replacement(mut self, options: &'a CooksReplacementOptions) -> Self {
        self.replacement = Replacement::Cooks(options);
        self
    }

    /// Stop after a stage.
    pub fn stop_after(mut self, stage: Stage) -> Self {
        self.stop_after = stage;
        self
    }

    /// Override the builder's dispersion trend type for this run.
    pub fn fit_type(mut self, fit_type: FitType) -> Self {
        self.fit_type = Some(fit_type);
        self
    }
}

/// The output of [`DeseqBuilder::run`](super::DeseqBuilder::run).
///
/// With Cook's replacement the original fit, the refit and the plan are kept
/// separately; [`RunOutput::into_fit_and_results`] merges them the way
/// DESeq2's `refitWithoutOutliers` does.
// The variants differ in size (a fit versus two fits and their results), but
// one value is returned per run and moved a few times, so boxing would only
// add an allocation and an extra dereference for callers.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
#[must_use]
#[non_exhaustive]
pub enum RunOutput {
    /// The fit state after an early [`Stage`].
    Fit(DeseqFit),
    /// A full run without replacement.
    Results {
        /// Final fit state.
        fit: DeseqFit,
        /// Result rows.
        results: DeseqResults,
    },
    /// A full run with Cook's replacement.
    Replacement(CooksReplacementTestOutput),
}

impl RunOutput {
    /// The final result rows, if the run reached [`Stage::Results`].
    pub fn results(&self) -> Option<&DeseqResults> {
        match self {
            Self::Fit(_) => None,
            Self::Results { results, .. } => Some(results),
            Self::Replacement(output) => Some(output.results()),
        }
    }

    /// The fit state of an early stage or of a run without replacement.
    /// With replacement this is the original (pre-replacement) fit.
    pub fn fit(&self) -> &DeseqFit {
        match self {
            Self::Fit(fit) | Self::Results { fit, .. } => fit,
            Self::Replacement(output) => match output {
                CooksReplacementTestOutput::Wald(output) => &output.original_fit,
                CooksReplacementTestOutput::Lrt(output) => &output.original_fit,
            },
        }
    }

    /// The fit (merged with replacement refits) and the final result rows.
    pub fn into_fit_and_results(self) -> Result<(DeseqFit, Option<DeseqResults>), DeseqError> {
        match self {
            Self::Fit(fit) => Ok((fit, None)),
            Self::Results { fit, results } => Ok((fit, Some(results))),
            Self::Replacement(output) => {
                let (fit, results, _plan) = output.into_merged_fit()?;
                Ok((fit, Some(results)))
            }
        }
    }
}
