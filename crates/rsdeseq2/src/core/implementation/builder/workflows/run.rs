//! The single `run` entry point (design note 4.2).
//!
//! `run` validates an [`Analysis`] once, resolves the design and the result
//! request to a canonical form, and dispatches to the leaf workflows.

use std::borrow::Cow;

use super::super::{
    Analysis, ContrastSpec, CooksReplacementOptions, CooksReplacementTestOutput, CountMatrix,
    DeseqBuilder, DeseqError, DesignInput, DesignMatrix, DispersionRequest, FactorLevelContrast,
    MeanModel, Replacement, ResultColumn, ResultsContrast, RunOutput, Stage, TestRequest,
    default_results_coefficient, factor_level_contrast_from_model_frame,
    factor_level_contrast_from_sample_levels, resolve_coefficient_index,
};

/// A result request after names, defaults and `results()` contrasts are resolved.
enum Column<'a> {
    Coefficient(usize),
    Numeric(&'a [f64]),
    Spec(Cow<'a, ContrastSpec>),
    FactorLevel(FactorLevelContrast<'a>),
}

fn options_error(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("run: {reason}"),
    }
}

fn unsupported(feature: &str) -> DeseqError {
    DeseqError::UnsupportedFeature {
        feature: format!("run: {feature}"),
    }
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
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.run(analysis));
        }
        validate_analysis(analysis)?;
        let builder = match analysis.fit_type {
            Some(fit_type) if fit_type != self.fit_type => {
                Cow::Owned(self.clone().fit_type(fit_type))
            }
            _ => Cow::Borrowed(self),
        };
        builder.run_validated(analysis)
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
                        return Err(options_error(
                            "an LRT reduced formula needs a formula full design",
                        ));
                    }
                };
                self.run_matrix(analysis, design, reduced)
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
                        return Err(options_error("an LRT needs a reduced design"));
                    }
                };
                self.with_formula_offsets(counts, &full)?.run_matrix(
                    analysis,
                    &full.design.standard_design,
                    reduced.as_deref(),
                )
            }
        }
    }

    fn run_matrix(
        &self,
        analysis: &Analysis<'_>,
        design: &DesignMatrix,
        reduced: Option<&DesignMatrix>,
    ) -> Result<RunOutput, DeseqError> {
        let counts = analysis.counts;
        let builder = match (analysis.dispersion, analysis.means) {
            (DispersionRequest::Estimate, MeanModel::Auto) => {
                self.with_automatic_dispersion_means()
            }
            _ => Cow::Borrowed(self),
        };
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
                if replacement.is_some() {
                    return Err(unsupported("Cook's replacement with linear-mu dispersions"));
                }
                builder.run_linear_mu(counts, design, reduced, column)
            }
            (DispersionRequest::Estimate, _) => {
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
                        return Err(options_error(
                            "results contrast takes either sample levels or a model frame",
                        ));
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
                        let sample_levels = sample_levels.ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "character results contrast requires sample levels for contrastAllZero"
                                .to_string(),
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

fn replaced(output: CooksReplacementTestOutput) -> RunOutput {
    RunOutput::Replacement(output)
}

fn results(
    output: Result<(super::super::DeseqFit, super::super::DeseqResults), DeseqError>,
) -> Result<RunOutput, DeseqError> {
    output.map(|(fit, results)| RunOutput::Results { fit, results })
}

fn validate_analysis(analysis: &Analysis<'_>) -> Result<(), DeseqError> {
    let early = analysis.stop_after != Stage::Results;
    if matches!(analysis.design, DesignInput::None) && analysis.stop_after != Stage::SizeFactors {
        return Err(options_error(
            "a run without a design must stop after size factors",
        ));
    }
    if matches!(
        analysis.test,
        TestRequest::Lrt {
            reduced: DesignInput::None
        }
    ) {
        return Err(options_error("an LRT needs a reduced design"));
    }
    if let DispersionRequest::Fixed(dispersions) = analysis.dispersion {
        if early {
            return Err(options_error(
                "fixed dispersions cannot stop at an estimation stage",
            ));
        }
        if analysis.means != MeanModel::Auto {
            return Err(options_error("fixed dispersions do not use a mean model"));
        }
        if dispersions.len() != analysis.counts.n_genes() {
            return Err(DeseqError::InvalidDimensions {
                context: "run fixed dispersions".to_string(),
                expected: analysis.counts.n_genes(),
                actual: dispersions.len(),
            });
        }
    }
    if early {
        if !matches!(analysis.replacement, Replacement::Off) {
            return Err(options_error("Cook's replacement needs a full run"));
        }
        if !matches!(analysis.result, ResultColumn::Default) {
            return Err(options_error("a result column needs a full run"));
        }
        if !matches!(analysis.test, TestRequest::Wald) {
            return Err(options_error("an LRT needs a full run"));
        }
    }
    Ok(())
}
