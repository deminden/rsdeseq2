use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};
impl DeseqBuilder {
    pub(super) fn leaf_wald_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_pipeline_wald_coefficient(design, coefficient)?;
            let fit = self.leaf_map_dispersions_linear_mu(counts, design)?;
            self.attach_native_wald(counts, design, coefficient, fit)
        })
    }

    pub(super) fn leaf_wald_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let fit = self.leaf_map_dispersions_linear_mu(counts, design)?;
            self.attach_native_wald_contrast(counts, design, contrast, None, fit)
        })
    }

    pub(super) fn leaf_wald_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let (fit, mut results) =
                self.leaf_wald_linear_mu_contrast(counts, design, &numeric_contrast)?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_wald_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let fit = self.leaf_map_dispersions_linear_mu(counts, design)?;
            let (fit, mut results) = self.attach_native_wald_contrast(
                counts,
                design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_wald_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_pipeline_wald_coefficient(design, coefficient)?;
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, design)?;
            self.attach_native_wald(counts, design, coefficient, fit)
        })
    }

    pub(super) fn leaf_wald_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, design)?;
            self.attach_native_wald_contrast(counts, design, contrast, None, fit)
        })
    }

    pub(super) fn leaf_wald_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let (fit, mut results) =
                self.leaf_wald_linear_mu_contrast_parametric(counts, design, &numeric_contrast)?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_wald_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, design)?;
            let (fit, mut results) = self.attach_native_wald_contrast(
                counts,
                design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_wald_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_pipeline_wald_coefficient(design, coefficient)?;
            let fit = self.leaf_map_dispersions_glm_mu(counts, design)?;
            self.attach_native_wald(counts, design, coefficient, fit)
        })
    }

    pub(super) fn leaf_wald_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let fit = self.leaf_map_dispersions_glm_mu(counts, design)?;
            self.attach_native_wald_contrast(counts, design, contrast, None, fit)
        })
    }

    pub(super) fn leaf_wald_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(design, contrast)?;
            let (fit, mut results) =
                self.leaf_wald_glm_mu_contrast(counts, design, &numeric_contrast)?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_wald_glm_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let fit = raw_builder.leaf_map_dispersions_glm_mu(counts, design)?;
            let (fit, mut results) = self.attach_native_wald_contrast(
                counts,
                design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let cooks = fit.cooks.as_ref().ok_or_else(|| {
                missing_prerequisite(
                    ErrorKind::Options,
                    "factor-level Cook's filtering",
                    "Cook's distances",
                )
            })?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                Some(cooks),
                design,
                self.model_frame.as_ref(),
                Some(contrast),
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_linear_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_linear_mu(counts, full_design)?;
            self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
        })
    }

    pub(super) fn leaf_lrt_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_linear_mu(counts, full_design)?;
            self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
                None,
                fit,
            )
        })
    }

    pub(super) fn leaf_lrt_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let (fit, mut results) = self.leaf_lrt_linear_mu_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
            )?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let fit = self.leaf_map_dispersions_linear_mu(counts, full_design)?;
            let (fit, mut results) = self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, full_design)?;
            self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
        })
    }

    pub(super) fn leaf_lrt_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, full_design)?;
            self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
                None,
                fit,
            )
        })
    }

    pub(super) fn leaf_lrt_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let (fit, mut results) = self.leaf_lrt_linear_mu_contrast_parametric(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
            )?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let fit = self.leaf_map_dispersions_linear_mu_parametric(counts, full_design)?;
            let (fit, mut results) = self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_glm_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_glm_mu(counts, full_design)?;
            self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_glm_mu(counts, full_design)?;
            self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
                None,
                fit,
            )
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let (fit, mut results) = self.leaf_lrt_glm_mu_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
            )?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_glm_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let raw_builder = self
                .clone()
                .disable_cooks_cutoff()
                .disable_independent_filtering();
            let fit = raw_builder.leaf_map_dispersions_glm_mu(counts, full_design)?;
            let (fit, mut results) = self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let cooks_cutoff = resolve_cooks_cutoff(
                self.cooks_cutoff,
                full_design.n_samples(),
                full_design.n_coefficients(),
            )?;
            let cooks = fit.cooks.as_ref().ok_or_else(|| {
                missing_prerequisite(
                    ErrorKind::Options,
                    "factor-level Cook's filtering",
                    "Cook's distances",
                )
            })?;
            apply_cooks_cutoff_for_design(
                &mut results,
                cooks_cutoff,
                counts,
                Some(cooks),
                full_design,
                self.model_frame.as_ref(),
                Some(contrast),
            )?;
            apply_independent_filtering(&mut results, &self.independent_filtering_options)?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_glm_mu_parametric(counts, full_design)?;
            self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let fit = self.leaf_map_dispersions_glm_mu_parametric(counts, full_design)?;
            self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                contrast,
                None,
                fit,
            )
        })
    }

    pub(super) fn leaf_lrt_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            let numeric_contrast = resolve_contrast(full_design, contrast)?;
            let (fit, mut results) = self.leaf_lrt_glm_mu_contrast_parametric(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
            )?;
            results.set_resolved_contrast_metadata(
                contrast.result_name(),
                contrast.comparison(),
                &numeric_contrast,
            );
            Ok((fit, results))
        })
    }

    pub(super) fn leaf_lrt_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.in_pool(|| {
            validate_lrt_designs(full_design, reduced_design)?;
            let contrast_spec = match contrast.reference {
                Some(reference) => ContrastSpec::factor_level_with_reference(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                    reference,
                ),
                None => ContrastSpec::factor_level(
                    contrast.factor,
                    contrast.numerator,
                    contrast.denominator,
                ),
            };
            let numeric_contrast = resolve_contrast(full_design, &contrast_spec)?;
            let contrast_all_zero = contrast_all_zero_factor_levels(
                counts,
                contrast.sample_levels,
                contrast.numerator,
                contrast.denominator,
            )?;
            let fit = self.leaf_map_dispersions_glm_mu_parametric(counts, full_design)?;
            let (fit, mut results) = self.attach_native_lrt_contrast(
                counts,
                full_design,
                reduced_design,
                &numeric_contrast,
                Some(&contrast_all_zero),
                fit,
            )?;
            let (result_name, comparison) = factor_level_result_metadata(contrast);
            results.set_resolved_contrast_metadata(result_name, comparison, &numeric_contrast);
            Ok((fit, results))
        })
    }
}
