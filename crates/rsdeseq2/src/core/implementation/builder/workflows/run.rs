//! The single `run` entry point (design note 4.2).
//!
//! `run` validates an [`Analysis`] once, resolves the design and the result
//! request to a canonical form, and dispatches to the leaf workflows.

use crate::errors::{ErrorKind, RunConflict, invalid_dimensions, missing_prerequisite};
use std::borrow::Cow;

use super::super::{
    Analysis, ContrastSpec, CooksReplacementLrtOutput, CooksReplacementOptions,
    CooksReplacementTestOutput, CooksReplacementWaldOutput, CountMatrix, DeseqBuilder, DeseqError,
    DeseqFit, DeseqResults, DesignInput, DesignMatrix, DispersionRequest, FactorLevelContrast,
    FitType, MeanModel, Replacement, ResultColumn, ResultsContrast, RunOutput, Stage, TestRequest,
    TestType, default_results_coefficient, factor_level_contrast_from_model_frame,
    factor_level_contrast_from_sample_levels, resolve_coefficient_index,
};

/// A result request after names, defaults and `results()` contrasts are resolved.
enum Column<'a> {
    Coefficient(usize),
    Numeric(&'a [f64]),
    Spec(Cow<'a, ContrastSpec>),
    FactorLevel(FactorLevelContrast<'a>),
}

fn invalid_run(conflict: RunConflict) -> DeseqError {
    DeseqError::InvalidRun { conflict }
}

fn unsupported(feature: &str) -> DeseqError {
    DeseqError::UnsupportedFeature {
        feature: format!("run: {feature}"),
    }
}

fn retain_formula_normalization_factors(
    output: &mut RunOutput,
    counts: &CountMatrix,
) -> Result<(), DeseqError> {
    let (fit, refit) = match output {
        RunOutput::Fit(fit) | RunOutput::Results { fit, .. } => (fit, None),
        RunOutput::BetaPrior(output) => (&mut output.fit, None),
        RunOutput::Replacement(CooksReplacementTestOutput::Wald(output)) => {
            (&mut output.original_fit, output.refit_fit.as_mut())
        }
        RunOutput::Replacement(CooksReplacementTestOutput::Lrt(output)) => {
            (&mut output.original_fit, output.refit_fit.as_mut())
        }
    };
    for fit in std::iter::once(fit).chain(refit) {
        fit.normalization_factors = Some(
            crate::normalization::normalization_factors_from_size_factors(
                counts,
                &fit.size_factors,
            )?,
        );
        fit.size_factors.fill(1.0);
    }
    Ok(())
}

impl DeseqBuilder {
    /// Run one analysis.
    ///
    /// This is the single entry point for the DESeq2 pipeline: every axis of
    /// the run (design, test, dispersion source, mean model, reported column,
    /// Cook's replacement and early stopping) is a field of [`Analysis`].
    /// Invalid combinations are rejected before any work is done.
    ///
    /// ```
    /// # use rsdeseq2::prelude::*;
    /// # fn main() -> Result<(), DeseqError> {
    /// let counts = CountMatrix::from_row_major_u32(
    ///     4,
    ///     6,
    ///     vec![
    ///         10, 12, 11, 30, 33, 29, //
    ///         5, 6, 7, 5, 7, 6, //
    ///         100, 90, 95, 120, 110, 130, //
    ///         40, 45, 38, 20, 22, 25,
    ///     ],
    /// )?;
    /// let design = DesignMatrix::from_row_major(
    ///     6,
    ///     2,
    ///     vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ///     Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    /// )?;
    /// let dispersions = [0.1; 4];
    /// let analysis = Analysis::new(&counts, &design)
    ///     .result(ResultColumn::Name("condition_B_vs_A"))
    ///     .fixed_dispersion(&dispersions);
    /// let output = DeseqBuilder::new().run(&analysis)?;
    /// assert_eq!(output.results().unwrap().rows.len(), 4);
    /// # Ok(())
    /// # }
    /// ```
    pub fn run(&self, analysis: &Analysis<'_>) -> Result<RunOutput, DeseqError> {
        self.in_pool(|| {
            validate_analysis(analysis)?;
            self.run_validated(analysis)
        })
    }

    fn run_validated(&self, analysis: &Analysis<'_>) -> Result<RunOutput, DeseqError> {
        let counts = analysis.counts;
        match analysis.design {
            DesignInput::None => self
                .leaf_size_factors_and_base_means(counts)
                .map(RunOutput::Fit),
            DesignInput::Matrix(design) => {
                let reduced = match analysis.test {
                    TestRequest::Wald => None,
                    TestRequest::Lrt {
                        reduced: DesignInput::Matrix(reduced),
                    } => Some(reduced),
                    TestRequest::Lrt { .. } => {
                        return Err(invalid_run(RunConflict::ReducedFormulaNeedsFullFormula));
                    }
                };
                self.run_matrix(analysis, design, reduced)
            }
            DesignInput::Formula(formula) if analysis.beta_prior => {
                self.run_beta_prior(analysis, formula)
            }
            DesignInput::Formula(formula) => {
                let full = self.formula_design_for_fitting(formula)?;
                let reduced = match analysis.test {
                    TestRequest::Wald => None,
                    TestRequest::Lrt {
                        reduced: DesignInput::Formula(reduced),
                    } => Some(Cow::Owned(self.reduced_formula_design(&full, reduced)?)),
                    TestRequest::Lrt {
                        reduced: DesignInput::Matrix(reduced),
                    } => Some(Cow::Borrowed(reduced)),
                    TestRequest::Lrt { .. } => {
                        return Err(invalid_run(RunConflict::ReducedDesignRequired));
                    }
                };
                let builder = self.with_formula_offsets(counts, &full)?;
                let mut output = builder.run_matrix(
                    analysis,
                    &full.design.standard_design,
                    reduced.as_deref(),
                )?;
                if builder.formula_sample_factors {
                    retain_formula_normalization_factors(&mut output, counts)?;
                }
                Ok(output)
            }
        }
    }

    /// `DESeq(betaPrior = TRUE)`: the usual run up to MAP dispersions on the
    /// standard design, then the expanded-model prior fit.
    fn run_beta_prior(
        &self,
        analysis: &Analysis<'_>,
        formula: &str,
    ) -> Result<RunOutput, DeseqError> {
        let full = self.formula_design_for_fitting(formula)?;
        let builder = self.with_formula_offsets(analysis.counts, &full)?;
        let design = &full.design.standard_design;
        let mut levels = Vec::new();
        let column = builder.resolve_column(analysis.result, design, &mut levels)?;
        let contrast = match column {
            Column::Coefficient(index) => {
                let mut unit = vec![0.0; design.n_coefficients()];
                *unit.get_mut(index).ok_or_else(|| {
                    invalid_dimensions(
                        "beta-prior coefficient index",
                        design.n_coefficients().saturating_sub(1),
                        index,
                    )
                })? = 1.0;
                unit
            }
            Column::Numeric(values) => values.to_vec(),
            Column::Spec(spec) => crate::contrasts::resolve_contrast(design, &spec)?,
            Column::FactorLevel(_) => {
                return Err(unsupported(
                    "a beta-prior run with a sample-level factor contrast",
                ));
            }
        };
        if contrast.len() != design.n_coefficients() {
            return Err(DeseqError::InvalidDimensions {
                context: "beta-prior contrast".to_string(),
                expected: design.n_coefficients(),
                actual: contrast.len(),
            });
        }
        let map = Analysis {
            stop_after: Stage::Map,
            result: ResultColumn::Default,
            beta_prior: false,
            ..*analysis
        };
        let fit = builder.run_matrix(&map, design, None)?.into_fit_only()?;
        let mut output = RunOutput::BetaPrior(builder.wald_beta_prior(
            analysis.counts,
            &full,
            fit,
            &contrast,
        )?);
        if builder.formula_sample_factors {
            retain_formula_normalization_factors(&mut output, analysis.counts)?;
        }
        Ok(output)
    }

    fn run_matrix(
        &self,
        analysis: &Analysis<'_>,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
    ) -> Result<RunOutput, DeseqError> {
        // A model frame supplied with a `results()` contrast replaces the
        // stored one for this run: it is validated against the design and
        // recorded on the fit, as the `_from_model_frame` methods do.
        if let ResultColumn::Results {
            model_frame: Some(frame),
            ..
        } = analysis.result
        {
            let framed = self.clone().try_model_frame(frame.clone())?;
            return framed.run_matrix_framed(analysis, design, reduced);
        }
        self.run_matrix_framed(analysis, design, reduced)
    }

    fn run_matrix_framed(
        &self,
        analysis: &Analysis<'_>,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
    ) -> Result<RunOutput, DeseqError> {
        let counts = analysis.counts;
        let mut builder = match (analysis.dispersion, analysis.means) {
            (DispersionRequest::Estimate, MeanModel::Auto) => {
                self.with_automatic_dispersion_means()
            }
            _ => Cow::Borrowed(self),
        };
        // A per-run trend type that differs from the builder's. Fixed
        // dispersions and the stages before the trend do not use it.
        let trend = analysis
            .fit_type
            .filter(|fit_type| *fit_type != self.fit_type)
            .filter(|_| matches!(analysis.dispersion, DispersionRequest::Estimate))
            .filter(|_| !matches!(analysis.stop_after, Stage::SizeFactors | Stage::GeneWise));
        let linear = analysis.means == MeanModel::LinearMu;
        if let Some(fit_type) = trend {
            // Dedicated leaves exist for these; everything else gets a builder
            // copy with the trend type set.
            let leaf = match (analysis.stop_after, fit_type) {
                (Stage::Trend, _) => {
                    Some(builder.run_trend_stage(counts, design, fit_type, linear))
                }
                (Stage::Map, FitType::Parametric) => Some(
                    if linear {
                        builder.leaf_map_dispersions_linear_mu_parametric(counts, design)
                    } else {
                        builder.leaf_map_dispersions_glm_mu_parametric(counts, design)
                    }
                    .map(RunOutput::Fit),
                ),
                _ => None,
            };
            if let Some(output) = leaf {
                return output;
            }
            if !(analysis.stop_after == Stage::Results
                && fit_type == FitType::Parametric
                && matches!(analysis.replacement, Replacement::Off))
            {
                builder = Cow::Owned(builder.into_owned().fit_type(fit_type));
            }
        }
        let parametric_leaves =
            trend == Some(FitType::Parametric) && matches!(analysis.replacement, Replacement::Off);
        if analysis.stop_after != Stage::Results {
            return builder.run_stage(counts, design, analysis.stop_after, analysis.means);
        }
        let mut levels = Vec::new();
        let column = self.resolve_column(analysis.result, design, &mut levels)?;
        let replacement = match analysis.replacement {
            Replacement::Off => None,
            Replacement::Cooks(options) => Some(options),
        };
        match (analysis.dispersion, analysis.means) {
            (DispersionRequest::Fixed(dispersions), _) => {
                builder.run_fixed(counts, design, reduced, dispersions, column, replacement)
            }
            (DispersionRequest::Estimate, MeanModel::LinearMu) => {
                if parametric_leaves {
                    return builder.run_linear_mu_parametric(counts, design, reduced, column);
                }
                builder.run_linear_mu(counts, design, reduced, column)
            }
            (DispersionRequest::Estimate, _) => {
                if parametric_leaves {
                    return builder.run_glm_mu_parametric(counts, design, reduced, column);
                }
                builder.run_glm_mu(counts, design, reduced, column, replacement)
            }
        }
    }

    fn run_stage(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        stage: Stage,
        means: MeanModel,
    ) -> Result<RunOutput, DeseqError> {
        let linear = means == MeanModel::LinearMu;
        match stage {
            Stage::SizeFactors => self.leaf_size_factors_and_base_means_with_design(counts, design),
            Stage::GeneWise if linear => self.leaf_gene_wise_dispersions_linear_mu(counts, design),
            Stage::GeneWise => self.leaf_gene_wise_dispersions_glm_mu(counts, design),
            Stage::Trend if linear => self.leaf_dispersion_trend_linear_mu(counts, design),
            Stage::Trend => self.leaf_dispersion_trend_glm_mu(counts, design),
            Stage::Map if linear => self.leaf_map_dispersions_linear_mu(counts, design),
            Stage::Map => self.leaf_map_dispersions_glm_mu(counts, design),
            Stage::Results => unreachable!("run_stage handles early stages only"),
        }
        .map(RunOutput::Fit)
    }

    fn resolve_column<'a>(
        &'a self,
        request: ResultColumn<'a>,
        design: &DesignMatrix,
        levels: &'a mut Vec<String>,
    ) -> Result<Column<'a>, DeseqError> {
        Ok(match request {
            ResultColumn::Default => Column::Coefficient(default_results_coefficient(design)?),
            ResultColumn::Coefficient(index) => Column::Coefficient(index),
            ResultColumn::Name(name) => {
                Column::Coefficient(resolve_coefficient_index(design, name)?)
            }
            ResultColumn::Numeric(contrast) => Column::Numeric(contrast),
            ResultColumn::Spec(spec) => Column::Spec(Cow::Borrowed(spec)),
            ResultColumn::FactorLevel(contrast) => Column::FactorLevel(contrast),
            ResultColumn::Results {
                contrast,
                sample_levels,
                model_frame,
            } => {
                let from_frame = match (model_frame, sample_levels) {
                    (Some(frame), None) => factor_level_contrast_from_model_frame(contrast, frame)?,
                    (Some(_), Some(_)) => {
                        return Err(invalid_run(RunConflict::DuplicateContrastMetadata));
                    }
                    (None, None) => self.model_frame_factor_level_contrast(contrast)?,
                    (None, Some(_)) => None,
                };
                if let Some(factor_contrast) = from_frame {
                    return Ok(Column::FactorLevel(factor_contrast));
                }
                match contrast {
                    ResultsContrast::Character {
                        factor,
                        numerator,
                        denominator,
                        reference,
                    } => {
                        let sample_levels = sample_levels.ok_or_else(|| {
                            missing_prerequisite(
                                ErrorKind::Options,
                                "character results contrast",
                                "sample levels for contrastAllZero",
                            )
                        })?;
                        levels.extend(sample_levels.iter().cloned());
                        Column::FactorLevel(factor_level_contrast_from_sample_levels(
                            factor,
                            numerator,
                            denominator,
                            reference.as_deref(),
                            levels,
                        )?)
                    }
                    ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
                        Column::Spec(Cow::Owned(contrast.as_contrast_spec()))
                    }
                }
            }
        })
    }

    fn run_trend_stage(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        fit_type: FitType,
        linear: bool,
    ) -> Result<RunOutput, DeseqError> {
        match (fit_type, linear) {
            (FitType::Parametric, true) => {
                self.leaf_parametric_dispersion_trend_linear_mu(counts, design)
            }
            (FitType::Parametric, false) => {
                self.leaf_parametric_dispersion_trend_glm_mu(counts, design)
            }
            (FitType::Mean, true) => self.leaf_mean_dispersion_trend_linear_mu(counts, design),
            (FitType::Mean, false) => self.leaf_mean_dispersion_trend_glm_mu(counts, design),
            (FitType::Local, true) => self.leaf_local_dispersion_trend_linear_mu(counts, design),
            (FitType::Local, false) => self.leaf_local_dispersion_trend_glm_mu(counts, design),
        }
        .map(RunOutput::Fit)
    }

    fn run_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
        column: Column<'_>,
    ) -> Result<RunOutput, DeseqError> {
        let d = design;
        results(match (reduced, column) {
            (None, Column::Coefficient(c)) => self.leaf_wald_glm_mu_parametric(counts, d, c),
            (None, Column::Numeric(v)) => self.leaf_wald_glm_mu_contrast_parametric(counts, d, v),
            (None, Column::Spec(s)) => {
                self.leaf_wald_glm_mu_contrast_spec_parametric(counts, d, &s)
            }
            (None, Column::FactorLevel(f)) => {
                self.leaf_wald_glm_mu_factor_level_contrast_parametric(counts, d, f)
            }
            (Some(r), Column::Coefficient(c)) => self.leaf_lrt_glm_mu_parametric(counts, d, r, c),
            (Some(r), Column::Numeric(v)) => {
                self.leaf_lrt_glm_mu_contrast_parametric(counts, d, r, v)
            }
            (Some(r), Column::Spec(s)) => {
                self.leaf_lrt_glm_mu_contrast_spec_parametric(counts, d, r, &s)
            }
            (Some(r), Column::FactorLevel(f)) => {
                self.leaf_lrt_glm_mu_factor_level_contrast_parametric(counts, d, r, f)
            }
        })
    }

    fn run_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
        column: Column<'_>,
    ) -> Result<RunOutput, DeseqError> {
        let d = design;
        results(match (reduced, column) {
            (None, Column::Coefficient(c)) => self.leaf_wald_linear_mu_parametric(counts, d, c),
            (None, Column::Numeric(v)) => {
                self.leaf_wald_linear_mu_contrast_parametric(counts, d, v)
            }
            (None, Column::Spec(s)) => {
                self.leaf_wald_linear_mu_contrast_spec_parametric(counts, d, &s)
            }
            (None, Column::FactorLevel(f)) => {
                self.leaf_wald_linear_mu_factor_level_contrast_parametric(counts, d, f)
            }
            (Some(r), Column::Coefficient(c)) => {
                self.leaf_lrt_linear_mu_parametric(counts, d, r, c)
            }
            (Some(r), Column::Numeric(v)) => {
                self.leaf_lrt_linear_mu_contrast_parametric(counts, d, r, v)
            }
            (Some(r), Column::Spec(s)) => {
                self.leaf_lrt_linear_mu_contrast_spec_parametric(counts, d, r, &s)
            }
            (Some(r), Column::FactorLevel(f)) => {
                self.leaf_lrt_linear_mu_factor_level_contrast_parametric(counts, d, r, f)
            }
        })
    }

    fn run_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
        column: Column<'_>,
        replacement: Option<&CooksReplacementOptions>,
    ) -> Result<RunOutput, DeseqError> {
        use CooksReplacementTestOutput::{Lrt, Wald};
        let (d, r) = (design, reduced);
        Ok(match (r, replacement, column) {
            (None, None, Column::Coefficient(c)) => results(self.leaf_wald_glm_mu(counts, d, c))?,
            (None, None, Column::Numeric(v)) => {
                results(self.leaf_wald_glm_mu_contrast(counts, d, v))?
            }
            (None, None, Column::Spec(s)) => {
                results(self.leaf_wald_glm_mu_contrast_spec(counts, d, &s))?
            }
            (None, None, Column::FactorLevel(f)) => {
                results(self.leaf_wald_glm_mu_factor_level_contrast(counts, d, f))?
            }
            (None, Some(o), Column::Coefficient(c)) => replaced(Wald(
                self.leaf_wald_glm_mu_with_cooks_replacement(counts, d, c, o)?,
            )),
            (None, Some(o), Column::Numeric(v)) => replaced(Wald(
                self.leaf_wald_glm_mu_contrast_with_cooks_replacement(counts, d, v, o)?,
            )),
            (None, Some(o), Column::Spec(s)) => replaced(Wald(
                self.leaf_wald_glm_mu_contrast_spec_with_cooks_replacement(counts, d, &s, o)?,
            )),
            (None, Some(o), Column::FactorLevel(f)) => replaced(Wald(
                self.leaf_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
                    counts, d, f, o,
                )?,
            )),
            (Some(r), None, Column::Coefficient(c)) => {
                results(self.leaf_lrt_glm_mu(counts, d, r, c))?
            }
            (Some(r), None, Column::Numeric(v)) => {
                results(self.leaf_lrt_glm_mu_contrast(counts, d, r, v))?
            }
            (Some(r), None, Column::Spec(s)) => {
                results(self.leaf_lrt_glm_mu_contrast_spec(counts, d, r, &s))?
            }
            (Some(r), None, Column::FactorLevel(f)) => {
                results(self.leaf_lrt_glm_mu_factor_level_contrast(counts, d, r, f))?
            }
            (Some(r), Some(o), Column::Coefficient(c)) => replaced(Lrt(
                self.leaf_lrt_glm_mu_with_cooks_replacement(counts, d, r, c, o)?
            )),
            (Some(r), Some(o), Column::Numeric(v)) => replaced(Lrt(
                self.leaf_lrt_glm_mu_contrast_with_cooks_replacement(counts, d, r, v, o)?
            )),
            (Some(r), Some(o), Column::Spec(s)) => {
                replaced(Lrt(self
                    .leaf_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                        counts, d, r, &s, o,
                    )?))
            }
            (Some(r), Some(o), Column::FactorLevel(f)) => replaced(Lrt(self
                .leaf_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
                    counts, d, r, f, o,
                )?)),
        })
    }

    fn run_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
        column: Column<'_>,
    ) -> Result<RunOutput, DeseqError> {
        let d = design;
        results(match (reduced, column) {
            (None, Column::Coefficient(c)) => self.leaf_wald_linear_mu(counts, d, c),
            (None, Column::Numeric(v)) => self.leaf_wald_linear_mu_contrast(counts, d, v),
            (None, Column::Spec(s)) => self.leaf_wald_linear_mu_contrast_spec(counts, d, &s),
            (None, Column::FactorLevel(f)) => {
                self.leaf_wald_linear_mu_factor_level_contrast(counts, d, f)
            }
            (Some(r), Column::Coefficient(c)) => self.leaf_lrt_linear_mu(counts, d, r, c),
            (Some(r), Column::Numeric(v)) => self.leaf_lrt_linear_mu_contrast(counts, d, r, v),
            (Some(r), Column::Spec(s)) => self.leaf_lrt_linear_mu_contrast_spec(counts, d, r, &s),
            (Some(r), Column::FactorLevel(f)) => {
                self.leaf_lrt_linear_mu_factor_level_contrast(counts, d, r, f)
            }
        })
    }

    fn run_fixed(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
        dispersions: &[f64],
        column: Column<'_>,
        replacement: Option<&CooksReplacementOptions>,
    ) -> Result<RunOutput, DeseqError> {
        use CooksReplacementTestOutput::{Lrt, Wald};
        let (d, x) = (design, dispersions);
        Ok(match (reduced, replacement, column) {
            (None, None, Column::Coefficient(c)) => {
                results(self.leaf_fixed_dispersion_wald(counts, d, x, c))?
            }
            (None, None, Column::Numeric(v)) => {
                results(self.leaf_fixed_dispersion_wald_contrast(counts, d, x, v))?
            }
            (None, None, Column::Spec(s)) => {
                results(self.leaf_fixed_dispersion_wald_contrast_spec(counts, d, x, &s))?
            }
            (None, None, Column::FactorLevel(f)) => {
                results(self.leaf_fixed_dispersion_wald_factor_level_contrast(counts, d, x, f))?
            }
            (None, Some(o), Column::Coefficient(c)) => replaced(Wald(
                self.leaf_fixed_dispersion_wald_with_cooks_replacement(counts, d, x, c, o)?,
            )),
            (None, Some(o), Column::Numeric(v)) => replaced(Wald(
                self.leaf_fixed_dispersion_wald_contrast_with_cooks_replacement(
                    counts, d, x, v, o,
                )?,
            )),
            (None, Some(o), Column::Spec(s)) => replaced(Wald(
                self.leaf_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
                    counts, d, x, &s, o,
                )?,
            )),
            (None, Some(o), Column::FactorLevel(f)) => replaced(Wald(
                self.leaf_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
                    counts, d, x, f, o,
                )?,
            )),
            (Some(r), None, Column::Coefficient(c)) => {
                results(self.leaf_fixed_dispersion_lrt(counts, d, r, x, c))?
            }
            (Some(r), None, Column::Numeric(v)) => {
                results(self.leaf_fixed_dispersion_lrt_contrast(counts, d, r, x, v))?
            }
            (Some(r), None, Column::Spec(s)) => {
                results(self.leaf_fixed_dispersion_lrt_contrast_spec(counts, d, r, x, &s))?
            }
            (Some(r), None, Column::FactorLevel(f)) => {
                results(self.leaf_fixed_dispersion_lrt_factor_level_contrast(counts, d, r, x, f))?
            }
            (Some(r), Some(o), Column::Coefficient(c)) => replaced(Lrt(
                self.leaf_fixed_dispersion_lrt_with_cooks_replacement(counts, d, r, x, c, o)?
            )),
            (Some(r), Some(o), Column::Numeric(v)) => replaced(Lrt(self
                .leaf_fixed_dispersion_lrt_contrast_with_cooks_replacement(
                    counts, d, r, x, v, o,
                )?)),
            (Some(r), Some(o), Column::Spec(s)) => replaced(Lrt(self
                .leaf_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
                    counts, d, r, x, &s, o,
                )?)),
            (Some(r), Some(o), Column::FactorLevel(f)) => replaced(Lrt(self
                .leaf_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
                    counts, d, r, x, f, o,
                )?)),
        })
    }
}

impl DeseqBuilder {
    /// The test a top-level `fit*` method runs: the builder's test type,
    /// with its stored reduced design for an LRT.
    pub(super) fn top_level_test(&self) -> Result<TestRequest<'_>, DeseqError> {
        Ok(match self.test {
            TestType::Wald => TestRequest::Wald,
            TestType::Lrt => TestRequest::Lrt {
                reduced: DesignInput::Matrix(self.reduced_design_for_top_level_lrt()?),
            },
        })
    }
}

fn shape_error(expected: &str, output: &RunOutput) -> DeseqError {
    let actual = match output {
        RunOutput::Fit(_) => "an early-stage fit",
        RunOutput::Results { .. } => "a fit with results",
        RunOutput::Replacement(CooksReplacementTestOutput::Wald(_)) => "a Wald replacement output",
        RunOutput::Replacement(CooksReplacementTestOutput::Lrt(_)) => "an LRT replacement output",
        RunOutput::BetaPrior(_) => "a beta-prior output",
    };
    DeseqError::InvalidOptions {
        reason: format!("run returned {actual} where {expected} was expected"),
    }
}

impl RunOutput {
    pub(in super::super) fn into_results_pair(
        self,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self {
            Self::Results { fit, results } => Ok((fit, results)),
            other => Err(shape_error("a fit with results", &other)),
        }
    }

    pub(in super::super) fn into_fit_only(self) -> Result<DeseqFit, DeseqError> {
        match self {
            Self::Fit(fit) | Self::Results { fit, .. } => Ok(fit),
            other => Err(shape_error("a fit", &other)),
        }
    }

    pub(in super::super) fn into_replacement(
        self,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self {
            Self::Replacement(output) => Ok(output),
            other => Err(shape_error("a replacement output", &other)),
        }
    }

    pub(in super::super) fn into_wald_replacement(
        self,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self {
            Self::Replacement(CooksReplacementTestOutput::Wald(output)) => Ok(output),
            other => Err(shape_error("a Wald replacement output", &other)),
        }
    }

    pub(in super::super) fn into_lrt_replacement(
        self,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        match self {
            Self::Replacement(CooksReplacementTestOutput::Lrt(output)) => Ok(output),
            other => Err(shape_error("an LRT replacement output", &other)),
        }
    }
}

fn replaced(output: CooksReplacementTestOutput) -> RunOutput {
    RunOutput::Replacement(output)
}

fn results(output: Result<(DeseqFit, DeseqResults), DeseqError>) -> Result<RunOutput, DeseqError> {
    output.map(|(fit, results)| RunOutput::Results { fit, results })
}

fn validate_analysis(analysis: &Analysis<'_>) -> Result<(), DeseqError> {
    let early = analysis.stop_after != Stage::Results;
    if matches!(analysis.design, DesignInput::None) && analysis.stop_after != Stage::SizeFactors {
        return Err(invalid_run(RunConflict::DesignRequired));
    }
    if matches!(
        analysis.test,
        TestRequest::Lrt {
            reduced: DesignInput::None
        }
    ) {
        return Err(invalid_run(RunConflict::ReducedDesignRequired));
    }
    if let DispersionRequest::Fixed(dispersions) = analysis.dispersion {
        if early {
            return Err(invalid_run(RunConflict::FixedDispersionEarlyStop));
        }
        if analysis.means != MeanModel::Auto {
            return Err(invalid_run(RunConflict::FixedDispersionMeanModel));
        }
        if dispersions.len() != analysis.counts.n_genes() {
            return Err(DeseqError::InvalidDimensions {
                context: "run fixed dispersions".to_string(),
                expected: analysis.counts.n_genes(),
                actual: dispersions.len(),
            });
        }
    }
    if analysis.beta_prior {
        if !matches!(analysis.design, DesignInput::Formula(_)) {
            return Err(invalid_run(RunConflict::BetaPriorNeedsFormula));
        }
        if !matches!(analysis.test, TestRequest::Wald) {
            return Err(invalid_run(RunConflict::BetaPriorNeedsWald));
        }
        if !matches!(analysis.dispersion, DispersionRequest::Estimate) {
            return Err(invalid_run(RunConflict::BetaPriorNeedsEstimatedDispersion));
        }
        if early {
            return Err(invalid_run(RunConflict::BetaPriorEarlyStop));
        }
        if !matches!(analysis.replacement, Replacement::Off) {
            return Err(unsupported("Cook's replacement with the beta prior"));
        }
    }
    if early {
        if !matches!(analysis.replacement, Replacement::Off) {
            return Err(invalid_run(RunConflict::ReplacementEarlyStop));
        }
        if !matches!(analysis.result, ResultColumn::Default) {
            return Err(invalid_run(RunConflict::ResultColumnEarlyStop));
        }
        if !matches!(analysis.test, TestRequest::Wald) {
            return Err(invalid_run(RunConflict::LrtEarlyStop));
        }
    }
    if matches!(
        (analysis.design, analysis.test),
        (
            DesignInput::Matrix(_),
            TestRequest::Lrt {
                reduced: DesignInput::Formula(_)
            }
        )
    ) {
        return Err(invalid_run(RunConflict::ReducedFormulaNeedsFullFormula));
    }
    if matches!(
        analysis.result,
        ResultColumn::Results {
            sample_levels: Some(_),
            model_frame: Some(_),
            ..
        }
    ) {
        return Err(invalid_run(RunConflict::DuplicateContrastMetadata));
    }
    if analysis.means == MeanModel::LinearMu && !matches!(analysis.replacement, Replacement::Off) {
        return Err(unsupported("Cook's replacement with linear-mu dispersions"));
    }
    Ok(())
}
