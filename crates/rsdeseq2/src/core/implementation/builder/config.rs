impl Default for DeseqBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl DeseqBuilder {
    /// Construct a builder with conservative DESeq2-like defaults.
    pub fn new() -> Self {
        Self {
            fit_type: FitType::default(),
            test: TestType::default(),
            size_factor_options: SizeFactorOptions::default(),
            normalization_factors: None,
            observation_weights: None,
            observation_weight_options: ObservationWeightOptions::default(),
            execution_mode: ExecutionMode::default(),
            threads: None,
            reduced_design: None,
            model_frame: None,
            irls_options: IrlsOptions::default(),
            gene_wise_dispersion_options: GeneWiseDispersionOptions::default(),
            wald_test_options: WaldTestOptions::default(),
            cooks_cutoff: CooksCutoff::default(),
            independent_filtering_options: IndependentFilteringOptions::default(),
        }
    }

    /// Set future dispersion fit type.
    pub fn fit_type(mut self, fit_type: FitType) -> Self {
        self.fit_type = fit_type;
        self
    }

    /// Set future test type.
    pub fn test(mut self, test: TestType) -> Self {
        self.test = test;
        self
    }

    /// Set size-factor method.
    pub fn size_factor_method(mut self, method: SizeFactorMethod) -> Self {
        self.size_factor_options.method = method;
        self
    }

    /// Set all size-factor options.
    pub fn size_factor_options(mut self, options: SizeFactorOptions) -> Self {
        self.size_factor_options = options;
        self
    }

    /// Set supplied geometric means for frozen size-factor estimation.
    pub fn geometric_means(mut self, geo_means: Vec<f64>) -> Self {
        self.size_factor_options.geo_means = Some(geo_means);
        self
    }

    /// Set caller-supplied size factors, bypassing size-factor estimation.
    pub fn size_factors(mut self, size_factors: Vec<f64>) -> Self {
        self.size_factor_options.supplied_size_factors = Some(size_factors);
        self
    }

    /// Set caller-supplied gene/sample normalization factors.
    ///
    /// As in DESeq2, these count-scale factors preempt size factors for
    /// normalized counts and fixed-dispersion GLM offsets.
    pub fn normalization_factors(mut self, normalization_factors: RowMajorMatrix<f64>) -> Self {
        self.normalization_factors = Some(normalization_factors);
        self
    }

    /// Set caller-supplied gene/sample observation weights.
    ///
    /// Initial no-design stages use these weights directly for DESeq2-style
    /// weighted base metadata. Design-aware stages first row-normalize and
    /// check them with `getAndCheckWeights`-style preprocessing.
    pub fn observation_weights(mut self, observation_weights: RowMajorMatrix<f64>) -> Self {
        self.observation_weights = Some(observation_weights);
        self
    }

    /// Set DESeq2-style observation-weight preprocessing options.
    ///
    /// The `weight_threshold` is also used by weighted Cox-Reid dispersion
    /// fitting, matching DESeq2's single `weightThreshold` argument.
    pub fn observation_weight_options(mut self, options: ObservationWeightOptions) -> Self {
        self.gene_wise_dispersion_options.weight_threshold = options.weight_threshold;
        self.observation_weight_options = options;
        self
    }

    /// Set zero-based control-gene row indices.
    pub fn control_genes(mut self, control_genes: Vec<usize>) -> Self {
        self.size_factor_options.control_genes = Some(ControlGenes::Indices(control_genes));
        self
    }

    /// Set a logical control-gene mask with one value per gene.
    pub fn control_gene_mask(mut self, control_gene_mask: Vec<bool>) -> Self {
        self.size_factor_options.control_genes = Some(ControlGenes::Mask(control_gene_mask));
        self
    }

    /// Set execution mode.
    pub fn execution_mode(mut self, mode: ExecutionMode) -> Self {
        self.execution_mode = mode;
        self
    }

    /// Set the desired worker thread count for future parallel stages.
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = Some(threads);
        self
    }

    /// Store a reduced design matrix for top-level LRT workflows.
    pub fn reduced_design(mut self, reduced_design: DesignMatrix) -> Self {
        self.reduced_design = Some(reduced_design);
        self
    }

    /// Store owned formula/model-frame metadata for object-style result routing.
    ///
    /// Character `results(contrast=...)` requests can use this metadata to
    /// infer the factor reference and per-sample levels when call sites do not
    /// pass explicit sample-level vectors.
    pub fn model_frame(mut self, model_frame: FormulaModelFrame) -> Self {
        self.model_frame = Some(model_frame);
        self
    }

    /// Build a supported expanded formula design from stored model-frame metadata.
    ///
    /// This is the object-style companion to
    /// [`crate::design::expanded_formula_design_from_model_frame`].
    pub fn expanded_formula_design(
        &self,
        formula: &str,
    ) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
        let model_frame = self.model_frame.as_ref().ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: "formula design construction requires builder model-frame metadata"
                    .to_string(),
            }
        })?;
        expanded_formula_design_from_model_frame(formula, model_frame)
    }

    /// Build a supported expanded formula design plus formula offsets from stored metadata.
    pub fn expanded_formula_design_with_offsets(
        &self,
        formula: &str,
    ) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
        let model_frame = self.model_frame.as_ref().ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: "formula design construction requires builder model-frame metadata"
                    .to_string(),
            }
        })?;
        expanded_formula_design_with_offsets_from_model_frame(formula, model_frame)
    }

    /// Set IRLS options for fixed-dispersion GLM fitting.
    pub fn irls_options(mut self, options: IrlsOptions) -> Self {
        self.irls_options = options;
        self
    }

    /// Set options for the current linear-mu gene-wise dispersion estimator.
    ///
    /// The `weight_threshold` is also used for observation-weight
    /// preprocessing, matching DESeq2's single `weightThreshold` argument.
    pub fn gene_wise_dispersion_options(mut self, options: GeneWiseDispersionOptions) -> Self {
        self.observation_weight_options.weight_threshold = options.weight_threshold;
        self.gene_wise_dispersion_options = options;
        self
    }

    /// Set Wald p-value options.
    pub fn wald_test_options(mut self, options: WaldTestOptions) -> Self {
        self.wald_test_options = options;
        self
    }

    /// Set DESeq2-style selected-coefficient LFC threshold testing.
    pub fn wald_lfc_threshold(mut self, threshold: f64, alternative: WaldAlternative) -> Self {
        self.wald_test_options.lfc_threshold = threshold;
        self.wald_test_options.alternative = alternative;
        self
    }

    /// Use DESeq2 `useT=TRUE` with residual degrees of freedom.
    pub fn wald_t_residual_degrees_of_freedom(mut self) -> Self {
        self.wald_test_options.pvalue_type =
            WaldTestOptions::t_residual_degrees_of_freedom().pvalue_type;
        self
    }

    /// Use DESeq2 `useT=TRUE` with one supplied degrees-of-freedom value.
    pub fn wald_t_degrees_of_freedom(mut self, degrees_of_freedom: f64) -> Self {
        self.wald_test_options.pvalue_type =
            WaldTestOptions::t_degrees_of_freedom(degrees_of_freedom).pvalue_type;
        self
    }

    /// Use DESeq2 `useT=TRUE` with one degrees-of-freedom value per gene.
    pub fn wald_t_per_gene_degrees_of_freedom(mut self, degrees_of_freedom: Vec<f64>) -> Self {
        self.wald_test_options.pvalue_type =
            WaldTestOptions::t_per_gene_degrees_of_freedom(degrees_of_freedom).pvalue_type;
        self
    }

    /// Set Cook's distance p-value filtering behavior for result rows.
    pub fn cooks_cutoff(mut self, cutoff: CooksCutoff) -> Self {
        self.cooks_cutoff = cutoff;
        self
    }

    /// Disable Cook's distance p-value filtering.
    pub fn disable_cooks_cutoff(mut self) -> Self {
        self.cooks_cutoff = CooksCutoff::Disabled;
        self
    }

    /// Use an explicit Cook's distance cutoff.
    pub fn cooks_cutoff_threshold(mut self, cutoff: f64) -> Self {
        self.cooks_cutoff = CooksCutoff::Threshold(cutoff);
        self
    }

    /// Set independent-filtering options for result-row assembly.
    pub fn independent_filtering_options(mut self, options: IndependentFilteringOptions) -> Self {
        self.independent_filtering_options = options;
        self
    }

    /// Disable independent filtering and use regular BH adjustment.
    pub fn disable_independent_filtering(mut self) -> Self {
        self.independent_filtering_options.enabled = false;
        self
    }

    /// Set the alpha used to select the independent-filtering threshold.
    pub fn independent_filtering_alpha(mut self, alpha: f64) -> Self {
        self.independent_filtering_options.alpha = alpha;
        self
    }

    /// Set an explicit independent-filtering theta grid.
    pub fn independent_filtering_theta(mut self, theta: Vec<f64>) -> Self {
        self.independent_filtering_options.theta = Some(theta);
        self
    }

    /// Current fit type option.
    pub fn current_fit_type(&self) -> FitType {
        self.fit_type
    }

    /// Current test option.
    pub fn current_test(&self) -> TestType {
        self.test
    }

    /// Current execution mode.
    pub fn current_execution_mode(&self) -> ExecutionMode {
        self.execution_mode
    }

    /// Requested thread count.
    pub fn requested_threads(&self) -> Option<usize> {
        self.threads
    }

    /// Current reduced design for top-level LRT workflows, if supplied.
    pub fn current_reduced_design(&self) -> Option<&DesignMatrix> {
        self.reduced_design.as_ref()
    }

    /// Current formula/model-frame metadata, if supplied.
    pub fn current_model_frame(&self) -> Option<&FormulaModelFrame> {
        self.model_frame.as_ref()
    }

    /// Current IRLS options.
    pub fn current_irls_options(&self) -> IrlsOptions {
        self.irls_options.clone()
    }

    /// Current gene-wise dispersion options.
    pub fn current_gene_wise_dispersion_options(&self) -> GeneWiseDispersionOptions {
        self.gene_wise_dispersion_options
    }

    /// Current Wald p-value options.
    pub fn current_wald_test_options(&self) -> &WaldTestOptions {
        &self.wald_test_options
    }

    /// Current Cook's cutoff option.
    pub fn current_cooks_cutoff(&self) -> CooksCutoff {
        self.cooks_cutoff
    }

    /// Current independent-filtering options.
    pub fn current_independent_filtering_options(&self) -> &IndependentFilteringOptions {
        &self.independent_filtering_options
    }

    /// Current size-factor options.
    pub fn current_size_factor_options(&self) -> &SizeFactorOptions {
        &self.size_factor_options
    }

    /// Current caller-supplied normalization factors, if any.
    pub fn current_normalization_factors(&self) -> Option<&RowMajorMatrix<f64>> {
        self.normalization_factors.as_ref()
    }

    /// Current caller-supplied observation weights, if any.
    pub fn current_observation_weights(&self) -> Option<&RowMajorMatrix<f64>> {
        self.observation_weights.as_ref()
    }

    /// Current observation-weight preprocessing options.
    pub fn current_observation_weight_options(&self) -> ObservationWeightOptions {
        self.observation_weight_options
    }

    fn model_frame_factor_level_contrast<'a>(
        &'a self,
        contrast: &'a ResultsContrast,
    ) -> Result<Option<FactorLevelContrast<'a>>, DeseqError> {
        match self.model_frame.as_ref() {
            Some(model_frame) => factor_level_contrast_from_model_frame(contrast, model_frame),
            None => Ok(None),
        }
    }

    fn model_frame_factor_level_contrast_for_coefficient<'a>(
        &'a self,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Option<FactorLevelContrast<'a>> {
        let model_frame = self.model_frame.as_ref()?;
        let [column] = model_frame.factors.as_slice() else {
            return None;
        };
        if column.sample_levels.len() != design.n_samples() {
            return None;
        }
        let mut observed_levels = Vec::<&str>::new();
        for level in &column.sample_levels {
            if !observed_levels.iter().any(|observed| *observed == level) {
                observed_levels.push(level);
            }
        }
        let [first_level, second_level] = observed_levels.as_slice() else {
            return None;
        };
        let reference = column
            .reference
            .as_deref()
            .or_else(|| {
                column
                    .levels
                    .as_ref()
                    .and_then(|levels| levels.first().map(String::as_str))
            })
            .unwrap_or(first_level);
        if reference != *first_level && reference != *second_level {
            return None;
        }
        let numerator = if reference == *first_level {
            *second_level
        } else {
            *first_level
        };
        let contrast = FactorLevelContrast {
            factor: &column.name,
            numerator,
            denominator: reference,
            reference: Some(reference),
            sample_levels: &column.sample_levels,
        };
        let contrast_spec =
            ContrastSpec::factor_level_with_reference(&column.name, numerator, reference, reference);
        let Ok(numeric_contrast) = resolve_contrast(design, &contrast_spec) else {
            return None;
        };
        if numeric_contrast
            .iter()
            .enumerate()
            .all(|(idx, value)| {
                if idx == coefficient {
                    (*value - 1.0).abs() <= f64::EPSILON
                } else {
                    value.abs() <= f64::EPSILON
                }
            })
        {
            Some(contrast)
        } else {
            None
        }
    }

    fn standard_design_from_formula_without_offsets(
        &self,
        formula: &str,
    ) -> Result<DesignMatrix, DeseqError> {
        Ok(self.formula_design_without_offsets(formula)?.design.standard_design)
    }

    fn formula_design_without_offsets(
        &self,
        formula: &str,
    ) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
        if formula_has_offset_terms(formula)? {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level formula workflows with formula offsets".to_string(),
            });
        }
        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        if !formula_design.offsets.is_empty() {
            return Err(DeseqError::UnsupportedFeature {
                feature: "top-level formula workflows with formula offsets".to_string(),
            });
        }
        Ok(formula_design)
    }
}
