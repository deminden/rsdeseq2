//! Explicit-path `fit*` methods (GLM-mu, linear-mu, fixed dispersion,
//! `_parametric` and stage methods), each a one-line [`DeseqBuilder::run`]
//! request. The workflows themselves are the private `leaf_*` methods.

#![allow(clippy::too_many_arguments)]

use super::*;

impl DeseqBuilder {
    /// Run the current GLM-mu native Wald path with limited Cook's replacement refit.
    ///
    /// This is an explicitly scoped analogue of the `replaceOutliers` /
    /// `refitWithoutOutliers` part of DESeq2 for the currently implemented
    /// native GLM-mu Wald branch. It preserves the original size factors,
    /// builds replacement counts from original Cook's distances, reruns the
    /// implemented GLM-mu dispersion/MAP/Wald path on replacement counts, and
    /// merges only rows marked for refit. Beta priors are applied by the
    /// results-level beta-prior functions; Bioconductor object slots are
    /// outside the Rust API.
    pub fn fit_wald_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(coefficient))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the current GLM-mu native Wald contrast path with limited Cook's replacement refit.
    ///
    /// This mirrors [`Self::fit_wald_glm_mu_with_cooks_replacement`] for
    /// primitive numeric contrasts: original and replacement-count refits both
    /// use the native GLM-mu contrast path, then only rows marked by the
    /// replacement plan are merged into the final result table.
    pub fn fit_wald_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run native GLM-mu Wald replacement refit for a named primitive contrast specification.
    pub fn fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run native GLM-mu Wald replacement refit for a caller-supplied factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the current GLM-mu native LRT path with limited Cook's replacement refit.
    ///
    /// This mirrors the scoped Wald replacement-refit path for the implemented
    /// native GLM-mu LRT branch: first fit on original counts with Cook's
    /// filtering disabled, replace eligible Cook's outlier counts, rerun the
    /// native GLM-mu dispersion/MAP/LRT path on replacement counts with the
    /// original size factors, and merge only rows marked by the refit plan.
    /// Broader DESeq2 behavior for contrasts, beta priors, and Bioconductor
    /// object metadata remains future work.
    pub fn fit_lrt_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(coefficient))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run the current GLM-mu native LRT contrast path with limited Cook's replacement refit.
    pub fn fit_lrt_glm_mu_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run native GLM-mu LRT replacement refit for a caller-supplied factor-level contrast.
    pub fn fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run native GLM-mu LRT replacement refit for a named primitive contrast specification.
    pub fn fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run a supplied-dispersion likelihood-ratio test pipeline.
    ///
    /// This mirrors the DESeq2 `nbinomLRT` shape for primitive matrices when
    /// dispersions are already available. The full-model beta fields are
    /// exposed in result rows alongside the model-level LRT statistic and
    /// p-value.
    pub fn fit_fixed_dispersion_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the supplied-dispersion LRT path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_lrt_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Coefficient(coefficient))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run a supplied-dispersion likelihood-ratio test and report a numeric contrast.
    ///
    /// This keeps the LRT model comparison unchanged while reporting contrast
    /// estimates and standard errors from the full model in result rows.
    pub fn fit_fixed_dispersion_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the supplied-dispersion LRT numeric-contrast path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_lrt_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run a supplied-dispersion LRT and report a named full-model contrast.
    pub fn fit_fixed_dispersion_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run supplied-dispersion LRT replacement refit for a named primitive contrast specification.
    pub fn fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run a supplied-dispersion LRT and report a factor-level full-model contrast.
    ///
    /// This resolves DESeq2-shaped coefficient names from the full design
    /// matrix and applies character-style `contrastAllZero` handling from
    /// caller-supplied sample levels. As in DESeq2 LRT result tables, the
    /// all-zero cleanup only zeroes the displayed LFC; the model-comparison
    /// statistic and p-values remain unchanged.
    pub fn fit_fixed_dispersion_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run supplied-dispersion LRT replacement refit for a factor-level full-model contrast.
    pub fn fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_lrt_replacement()
    }

    /// Run a supplied-dispersion Wald pipeline for one coefficient.
    ///
    /// Analogue of `nbinomWaldTest` with caller-supplied dispersions. Contrasts
    /// use the `fit_fixed_dispersion_wald_contrast*` variants and beta priors
    /// the results-level beta-prior functions.
    /// All-zero rows are skipped during GLM fitting and expanded back as
    /// missing numeric values, matching DESeq2's `buildMatrixWithNARows`
    /// pattern.
    pub fn fit_fixed_dispersion_wald(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run a supplied-dispersion Wald pipeline for a primitive numeric contrast.
    ///
    /// The contrast vector must contain one finite value per design
    /// coefficient. This is the Rust primitive analogue of DESeq2's numeric
    /// contrast path after R has resolved formula terms and coefficient names.
    pub fn fit_fixed_dispersion_wald_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the supplied-dispersion Wald coefficient path with Cook's replacement refit.
    ///
    /// This extends the replacement/refit machinery to the fixed-dispersion
    /// core path. Dispersions remain caller-supplied for both the original fit
    /// and the replacement-count refit; no native dispersion re-estimation or
    /// beta-prior refit is performed here.
    pub fn fit_fixed_dispersion_wald_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Coefficient(coefficient))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the supplied-dispersion Wald numeric-contrast path with Cook's replacement refit.
    pub fn fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Numeric(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run supplied-dispersion Wald replacement refit for a named primitive contrast specification.
    pub fn fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Spec(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run a supplied-dispersion Wald pipeline for a named primitive contrast specification.
    ///
    /// This resolves coefficient names and DESeq2-style positive/negative
    /// coefficient-name lists to a numeric contrast before calling
    /// `fit_fixed_dispersion_wald_contrast`.
    pub fn fit_fixed_dispersion_wald_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run a supplied-dispersion Wald pipeline for a factor-level contrast.
    ///
    /// This resolves DESeq2-shaped coefficient names from the design matrix and
    /// applies DESeq2-style character contrast all-zero handling using
    /// caller-supplied sample levels. Model-frame callers can use
    /// [`Self::fit_fixed_dispersion_wald_results_contrast_from_model_frame`].
    pub fn fit_fixed_dispersion_wald_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run supplied-dispersion Wald replacement refit for a factor-level contrast.
    pub fn fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .fixed_dispersion(dispersions)
                .result(ResultColumn::FactorLevel(contrast))
                .cooks_replacement(replacement_options),
        )?
        .into_wald_replacement()
    }

    /// Run the parametric GLM-mu native dispersion path and then a Wald test.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native Wald path for a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_glm_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the current linear-mu gene-wise dispersion estimator.
    ///
    /// This is a narrow Phase 3 stage for designs where DESeq2's
    /// `linearMu=TRUE` branch is appropriate. It estimates size factors,
    /// normalized counts, base row metadata, linear fitted means, and
    /// fixed-mean gene-wise dispersions. Cox-Reid correction, iterative GLM
    /// mean refits, trend fitting, and MAP shrinkage remain future stages.
    pub fn fit_gene_wise_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .stop_after(Stage::GeneWise),
        )?
        .into_fit_only()
    }

    /// Run the current GLM-mu gene-wise dispersion estimator.
    ///
    /// This is the first non-`linearMu` foundation for
    /// `estimateDispersionsGeneEst`: rough/moments starts are followed by
    /// fixed-dispersion NB GLM mean fitting and fixed-mean dispersion
    /// optimization. Builder-supplied observation weights are preprocessed in
    /// the same design-aware stage used by fixed-dispersion Wald/LRT paths and
    /// then passed into both the fixed-dispersion mean fit and fixed-mean
    /// dispersion objective.
    pub fn fit_gene_wise_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .stop_after(Stage::GeneWise),
        )?
        .into_fit_only()
    }

    /// Run linear-mu gene-wise, selected trend, prior variance, and MAP dispersion stages.
    ///
    /// This fills final `dispersion` values using the builder's `fit_type`.
    /// Parametric, local, and mean trends are implemented. It follows the implemented
    /// subset of DESeq2's
    /// `estimateDispersionsMAP(type="DESeq2")`: no observation weights and
    /// deterministic prior-variance estimation, including the low-df
    /// histogram/KL branch.
    pub fn fit_map_dispersions_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .stop_after(Stage::Map),
        )?
        .into_fit_only()
    }

    /// Run linear-mu gene-wise, parametric trend, prior variance, and MAP dispersion stages.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_map_dispersions_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .stop_after(Stage::Map),
        )?
        .into_fit_only()
    }

    /// Run GLM-mu gene-wise, selected trend, prior variance, and MAP dispersion stages.
    ///
    /// This fills final `dispersion` values using the builder's `fit_type`.
    /// Parametric, local, and mean trends are implemented. Builder-supplied
    /// observation weights flow through the GLM-mu gene-wise, MAP, and native
    /// Wald stages after DESeq2-style preprocessing.
    pub fn fit_map_dispersions_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .stop_after(Stage::Map),
        )?
        .into_fit_only()
    }

    /// Run GLM-mu gene-wise, parametric trend, prior variance, and MAP dispersion stages.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_map_dispersions_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .stop_after(Stage::Map),
        )?
        .into_fit_only()
    }

    /// Run only the implemented initial normalization stages.
    pub fn fit_size_factors_and_base_means(
        &self,
        counts: &CountMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(&Analysis::without_design(counts))?.into_fit_only()
    }

    /// Run initial normalization stages with design-aware observation-weight checks.
    ///
    /// This is useful for parity checks against DESeq2's early metadata when a
    /// `weights` assay is present: raw weights are used for `baseMean` and
    /// `baseVar`, weights are row-normalized for fitting checks, design/rank
    /// failures are recorded in `weights_fail`, and failed rows are marked in
    /// `all_zero` for downstream skipping.
    pub fn fit_size_factors_and_base_means_with_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(&Analysis::new(counts, design).stop_after(Stage::SizeFactors))?
            .into_fit_only()
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the parametric trend.
    ///
    /// This mirrors the implemented subset of `estimateDispersionsGeneEst`
    /// followed by `estimateDispersionsFit(fitType="parametric")`. It fills
    /// `dispGeneEst`, `dispFit`, and the linear fitted mean matrix, but it does
    /// not yet estimate prior variance or MAP dispersions.
    pub fn fit_parametric_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the mean trend.
    ///
    /// This mirrors the implemented subset of `estimateDispersionsGeneEst`
    /// followed by `estimateDispersionsFit(fitType="mean")`. It fills
    /// `dispGeneEst`, a constant `dispFit` for non-all-zero rows, and the
    /// linear fitted mean matrix. Prior variance and MAP dispersions are the
    /// next stage (`fit_map_dispersions_linear_mu`).
    pub fn fit_mean_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Mean)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the parametric trend.
    ///
    /// This mirrors the current non-`linearMu` gene-wise branch
    /// followed by `estimateDispersionsFit(fitType="parametric")`.
    pub fn fit_parametric_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the mean trend.
    pub fn fit_mean_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Mean)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Local)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Local)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run the implemented linear-mu dispersion trend path selected by `fit_type`.
    ///
    /// Dispatches to the parametric, local or mean trend.
    pub fn fit_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run the implemented GLM-mu dispersion trend path selected by `fit_type`.
    ///
    /// Dispatches to the parametric, local or mean trend.
    pub fn fit_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .stop_after(Stage::Trend),
        )?
        .into_fit_only()
    }

    /// Run the implemented native dispersion path and then a Wald test.
    ///
    /// Analogue of `DESeq(..., test="Wald")` with linear fitted means: size
    /// factors, base means, linear-mu gene-wise dispersions, the selected
    /// parametric, local or mean trend, prior variance, MAP dispersions, the GLM
    /// fit, Cook's distances and Wald results for one coefficient. Contrasts use
    /// the `fit_wald_linear_mu_contrast*` variants, observation weights the
    /// GLM-mu branch (`fit_wald_glm_mu`), and beta priors the results-level
    /// beta-prior functions.
    pub fn fit_wald_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native dispersion path and then a Wald test for a numeric contrast.
    ///
    /// This is the primitive contrast companion to [`Self::fit_wald_linear_mu`].
    /// It reuses the implemented linear-mu dispersion/MAP path, then reports
    /// the requested contrast from the final fixed-dispersion GLM fit.
    pub fn fit_wald_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native Wald path for a named primitive contrast specification.
    pub fn fit_wald_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native Wald path for a caller-supplied factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric native dispersion path and then a Wald test.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_wald_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native Wald path for a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_wald_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native dispersion path and then a Wald test.
    ///
    /// This mirrors `fit_wald_linear_mu` but uses the GLM-mu mean/dispersion
    /// alternation before trend, MAP, fixed-dispersion GLM, Cook's distances,
    /// and selected-coefficient Wald results. Builder-supplied observation
    /// weights are supported for this branch.
    pub fn fit_wald_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native dispersion path and then a Wald test for a numeric contrast.
    ///
    /// This follows the same native dispersion/MAP and final GLM fitting path
    /// as [`Self::fit_wald_glm_mu`], then reports a primitive numeric
    /// contrast over the fitted coefficient vector. Higher-level formula and
    /// factor handling remains caller or wrapper responsibility.
    pub fn fit_wald_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native Wald path for a named primitive contrast specification.
    ///
    /// This resolves coefficient names, coefficient-name lists, or supported
    /// factor-level coefficient shapes to a numeric contrast before running
    /// the implemented native GLM-mu Wald contrast path.
    pub fn fit_wald_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native Wald path for a factor-level contrast.
    ///
    /// In addition to resolving the coefficient contrast, this applies
    /// DESeq2-style character `contrastAllZero` handling using the supplied
    /// per-sample factor levels.
    pub fn fit_wald_glm_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the implemented linear-mu native dispersion path and then an LRT.
    ///
    /// This is a limited native analogue of `nbinomLRT`: the full design is
    /// used for the currently implemented linear-mu dispersion/MAP stages, then
    /// full and reduced fixed-dispersion GLMs are fit using those final
    /// dispersions.
    pub fn fit_lrt_linear_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native LRT path and report a full-model numeric contrast.
    ///
    /// The likelihood-ratio statistic and p-values remain the full-vs-reduced
    /// model comparison; the result table's effect-size columns come from the
    /// supplied contrast over the full-model coefficients.
    pub fn fit_lrt_linear_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the linear-mu native LRT path for a factor-level full-model contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native dispersion path and then an LRT.
    ///
    /// This compatibility-named entry point keeps the original parametric-only
    /// behavior even if the builder's `fit_type` is set to another value.
    pub fn fit_lrt_linear_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native LRT path and report a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_linear_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native LRT path and report a named contrast.
    pub fn fit_lrt_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric linear-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::LinearMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the implemented GLM-mu native dispersion path and then an LRT.
    ///
    /// Builder-supplied observation weights are supported for the GLM-mu branch
    /// through the same preprocessing used by native Wald.
    pub fn fit_lrt_glm_mu(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native LRT path and report a full-model numeric contrast.
    ///
    /// The model comparison, deviance, and p-values come from the LRT. The
    /// result table's effect-size columns use the supplied contrast over the
    /// fitted full-model coefficients, matching DESeq2's result-table shape.
    pub fn fit_lrt_glm_mu_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the GLM-mu native LRT path for a factor-level full-model contrast.
    ///
    /// In addition to resolving the coefficient contrast, this applies
    /// DESeq2-style character `contrastAllZero` handling using the supplied
    /// per-sample factor levels. For LRT result tables, only the reported LFC
    /// is zeroed; the full-vs-reduced statistic and p-values are preserved.
    pub fn fit_lrt_glm_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native dispersion path and then an LRT.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_glm_mu_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Coefficient(coefficient)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native LRT path and report a numeric contrast.
    ///
    /// This compatibility-named entry point keeps parametric behavior even if
    /// the builder's `fit_type` is set to another value.
    pub fn fit_lrt_glm_mu_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Numeric(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native LRT path and report a named contrast.
    pub fn fit_lrt_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::Spec(contrast)),
        )?
        .into_results_pair()
    }

    /// Run the parametric GLM-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.run(
            &Analysis::new(counts, full_design)
                .lrt(reduced_design)
                .means(MeanModel::GlmMu)
                .fit_type(FitType::Parametric)
                .result(ResultColumn::FactorLevel(contrast)),
        )?
        .into_results_pair()
    }
}
