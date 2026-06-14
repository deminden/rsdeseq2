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
    ///
    /// This setter intentionally preserves existing builder-style chaining and
    /// does not validate immediately. Wrapper and object-ingestion code should
    /// prefer [`Self::try_model_frame`], while formula-built helper paths store
    /// model frames that were already validated during formula construction.
    pub fn model_frame(mut self, model_frame: FormulaModelFrame) -> Self {
        self.model_frame = Some(model_frame);
        self
    }

    /// Validate and store owned formula/model-frame metadata.
    ///
    /// This checked companion to [`Self::model_frame`] is intended for
    /// wrapper/object ingestion paths where invalid metadata should fail before
    /// a later formula, contrast, or Cook's/refit route depends on it.
    pub fn try_model_frame(mut self, model_frame: FormulaModelFrame) -> Result<Self, DeseqError> {
        model_frame.validate()?;
        self.model_frame = Some(model_frame);
        Ok(self)
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
    ) -> Result<Option<FactorLevelContrast<'a>>, DeseqError> {
        let mut coefficient_contrast = vec![0.0; design.n_coefficients()];
        coefficient_contrast[coefficient] = 1.0;
        self.model_frame_factor_level_contrast_for_numeric_contrast(design, &coefficient_contrast)
    }

    fn model_frame_factor_level_contrast_for_numeric_contrast<'a>(
        &'a self,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<Option<FactorLevelContrast<'a>>, DeseqError> {
        if contrast.len() != design.n_coefficients() {
            return Ok(None);
        }
        let Some(model_frame) = self.model_frame.as_ref() else {
            return Ok(None);
        };
        model_frame.validate()?;
        let [column] = model_frame.factors.as_slice() else {
            return Ok(None);
        };
        if column.sample_levels.len() != design.n_samples() {
            return Ok(None);
        }
        let mut observed_levels = Vec::<&str>::new();
        for level in &column.sample_levels {
            if !observed_levels.iter().any(|observed| *observed == level) {
                observed_levels.push(level);
            }
        }
        let [first_level, second_level] = observed_levels.as_slice() else {
            return Ok(None);
        };
        let raw_reference = column
            .reference
            .as_deref()
            .or_else(|| {
                column
                    .levels
                    .as_ref()
                    .and_then(|levels| levels.first().map(String::as_str))
            })
            .unwrap_or(first_level);
        let reference =
            match resolve_observed_level_alias_for_builder(raw_reference, &observed_levels) {
                Some(reference) => reference,
                None => return Ok(None),
            };
        if reference != *first_level && reference != *second_level {
            return Ok(None);
        }
        let numerator = if reference == *first_level {
            *second_level
        } else {
            *first_level
        };
        let factor_contrast = FactorLevelContrast {
            factor: &column.name,
            numerator,
            denominator: reference,
            reference: Some(reference),
            sample_levels: &column.sample_levels,
        };
        let contrast_spec =
            ContrastSpec::factor_level_with_reference(&column.name, numerator, reference, reference);
        let Ok(numeric_contrast) = resolve_contrast(design, &contrast_spec) else {
            return Ok(None);
        };
        if numeric_contrasts_match(&numeric_contrast, contrast) {
            Ok(Some(factor_contrast))
        } else if numeric_contrasts_match_with_scale(&numeric_contrast, contrast, -1.0) {
            Ok(Some(FactorLevelContrast {
                factor: &column.name,
                numerator: reference,
                denominator: numerator,
                reference: Some(reference),
                sample_levels: &column.sample_levels,
            }))
        } else {
            Ok(None)
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

    fn with_formula_offsets(
        &self,
        counts: &CountMatrix,
        formula_design: &ExpandedFormulaDesignWithOffsets,
    ) -> Result<Self, DeseqError> {
        if formula_design.offsets.is_empty() {
            return Ok(self.clone().model_frame(formula_design.model_frame.clone()));
        }
        if formula_design.offsets.len() != counts.n_samples() {
            return Err(invalid_dimensions(
                "formula offsets",
                counts.n_samples(),
                formula_design.offsets.len(),
            ));
        }
        let offset_multipliers = formula_design
            .offsets
            .iter()
            .copied()
            .enumerate()
            .map(|(sample, offset)| {
                let multiplier = offset.exp();
                if multiplier.is_finite() && multiplier > 0.0 {
                    Ok(multiplier)
                } else {
                    Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula offset exponentiated to non-finite factor at sample {sample}"
                        ),
                    })
                }
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut builder = self.clone().model_frame(formula_design.model_frame.clone());
        let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
        match &self.normalization_factors {
            Some(factors) => {
                validate_normalization_factors(counts, factors)?;
                for gene in 0..counts.n_genes() {
                    let row = factors.row(gene)?;
                    for (sample, multiplier) in offset_multipliers.iter().copied().enumerate() {
                        values.push(row[sample] * multiplier);
                    }
                }
            }
            None => {
                let control_gene_indices = self
                    .size_factor_options
                    .control_genes
                    .as_ref()
                    .map(|control_genes| control_genes.to_indices(counts.n_genes()))
                    .transpose()?;
                let size_factors = match &self.size_factor_options.supplied_size_factors {
                    Some(size_factors) => size_factors.clone(),
                    None => estimate_size_factors_with_options(
                        counts,
                        self.size_factor_options.method,
                        self.size_factor_options.geo_means.as_deref(),
                        control_gene_indices.as_deref(),
                    )?,
                };
                normalized_counts(counts, &size_factors)?;
                for _gene in 0..counts.n_genes() {
                    for (sample, multiplier) in offset_multipliers.iter().copied().enumerate() {
                        values.push(size_factors[sample] * multiplier);
                    }
                }
                builder.size_factor_options.supplied_size_factors = Some(size_factors);
            }
        }
        builder.normalization_factors = Some(RowMajorMatrix::from_row_major(
            counts.n_genes(),
            counts.n_samples(),
            values,
        )?);
        Ok(builder)
    }
}

fn numeric_contrasts_match(expected: &[f64], observed: &[f64]) -> bool {
    numeric_contrasts_match_with_scale(expected, observed, 1.0)
}

fn numeric_contrasts_match_with_scale(expected: &[f64], observed: &[f64], scale: f64) -> bool {
    expected
        .iter()
        .zip(observed.iter())
        .all(|(expected, observed)| (scale * *expected - *observed).abs() <= f64::EPSILON)
}

fn resolve_observed_level_alias_for_builder<'a>(
    requested: &str,
    observed_levels: &[&'a str],
) -> Option<&'a str> {
    let exact = observed_levels
        .iter()
        .copied()
        .filter(|level| *level == requested)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [level] => return Some(*level),
        [] => {}
        _ => return None,
    }

    let matches = observed_levels
        .iter()
        .copied()
        .filter(|level| {
            design_name_candidates(level)
                .into_iter()
                .any(|candidate| candidate == requested)
        })
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [level] => Some(*level),
        _ => None,
    }
}

fn factor_level_contrast_from_sample_levels<'a>(
    factor: &'a str,
    numerator: &str,
    denominator: &str,
    reference: Option<&str>,
    sample_levels: &'a [String],
) -> Result<FactorLevelContrast<'a>, DeseqError> {
    let mut observed_levels = Vec::<&str>::new();
    for level in sample_levels {
        if !observed_levels.iter().any(|observed| *observed == level) {
            observed_levels.push(level);
        }
    }
    let numerator =
        resolve_sample_level_alias_for_builder(factor, numerator, "numerator", &observed_levels)?;
    let denominator = resolve_sample_level_alias_for_builder(
        factor,
        denominator,
        "denominator",
        &observed_levels,
    )?;
    if numerator == denominator {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor '{factor}' numerator and denominator resolve to the same level '{numerator}'"
            ),
        });
    }
    let reference = reference
        .map(|reference| {
            resolve_sample_level_alias_for_builder(factor, reference, "reference", &observed_levels)
        })
        .transpose()?;
    Ok(FactorLevelContrast {
        factor,
        numerator,
        denominator,
        reference,
        sample_levels,
    })
}

fn resolve_sample_level_alias_for_builder<'a>(
    factor: &str,
    requested: &str,
    role: &str,
    observed_levels: &[&'a str],
) -> Result<&'a str, DeseqError> {
    let exact = observed_levels
        .iter()
        .copied()
        .filter(|level| *level == requested)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [level] => return Ok(*level),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "factor '{factor}' {role} level '{requested}' appears more than once"
                ),
            });
        }
    }

    let matches = observed_levels
        .iter()
        .copied()
        .filter(|level| {
            design_name_candidates(level)
                .into_iter()
                .any(|candidate| candidate == requested)
        })
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [level] => Ok(*level),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' does not contain {role} level '{requested}'"),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor '{factor}' {role} level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}
