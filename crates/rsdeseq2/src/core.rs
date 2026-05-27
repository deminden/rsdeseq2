use crate::all_zero::{
    expand_gene_values_with_fill_rows, expand_gene_values_with_nan_rows,
    expand_matrix_with_nan_rows, mask_all_zero_values_with_nan_rows, nonzero_gene_indices,
};
use crate::contrasts::{
    contrast_all_zero_factor_levels, contrast_all_zero_numeric, resolve_contrast, ContrastSpec,
    FactorLevelContrast,
};
use crate::cooks::{
    calculate_cooks_distance, prepare_cooks_replacement_refit, CooksOutput, CooksRefitPlan,
    CooksReplacementOptions,
};
use crate::design::DesignMatrix;
use crate::dispersion::{
    estimate_dispersion_prior_variance, estimate_gene_wise_dispersions_glm_mu,
    estimate_gene_wise_dispersions_linear_mu, estimate_map_dispersions, fit_local_dispersion_trend,
    fit_mean_dispersion_trend, fit_parametric_dispersion_trend, DispersionTrendFit,
    GeneWiseDispersionInput, GeneWiseDispersionOptions, LocalDispersionTrendOptions,
    MapDispersionInput, MapDispersionOptions, MeanDispersionTrendOptions,
    ParametricDispersionTrendOptions,
};
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights,
    fit_intercept_only_fixed_dispersion_with_normalization_factors,
    fit_intercept_only_fixed_dispersion_with_weights, lrt_test,
    preprocess_observation_weights_with_options, wald_test_coefficient_with_options,
    wald_test_contrast_with_options, IrlsOptions, LrtOutput, NbinomGlmFit,
    ObservationWeightOptions, WaldAlternative, WaldContrastOutput, WaldOutput, WaldTestOptions,
};
use crate::independent_filtering::{apply_independent_filtering, IndependentFilteringOptions};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{
    base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
    estimate_size_factors_with_options, normalized_counts, normalized_counts_with_factors,
    validate_normalization_factors,
};
use crate::options::{
    ControlGenes, CooksCutoff, ExecutionMode, FitType, SizeFactorMethod, SizeFactorOptions,
    TestType,
};
use crate::results::{
    apply_cooks_cutoff, build_lrt_contrast_results, build_lrt_results, build_wald_contrast_results,
    build_wald_results_from_wald, resolve_cooks_cutoff, DeseqResultRow, DeseqResults,
};
use crate::transform::{
    fast_vst_eligible_count as count_fast_vst_eligible_rows,
    fast_vst_subset as build_fast_vst_subset, norm_transform,
    rlog_frozen_with_normalization_factors, rlog_frozen_with_size_factors,
    rlog_with_estimated_prior_and_normalization_factors,
    rlog_with_estimated_prior_and_size_factors,
    vst_with_dispersion_trend_and_normalization_factors,
    vst_with_dispersion_trend_and_size_factors, FastVstSubset, RlogMetadata, RlogOffsetMode,
    RlogOutput,
};

/// DESeq2's default number of genes used to fit the fast `vst()` trend.
pub const DEFAULT_FAST_VST_NSUB: usize = 1000;

/// Row-major genes x samples count matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct CountMatrix {
    n_genes: usize,
    n_samples: usize,
    counts: Vec<u32>,
    gene_names: Option<Vec<String>>,
    sample_names: Option<Vec<String>>,
}

impl CountMatrix {
    /// Build a count matrix from row-major `u32` values.
    pub fn from_row_major_u32(
        n_genes: usize,
        n_samples: usize,
        counts: Vec<u32>,
    ) -> Result<Self, DeseqError> {
        Self::from_row_major_u32_with_names(n_genes, n_samples, counts, None, None)
    }

    /// Build a count matrix from row-major `u64` values.
    ///
    /// The initial storage type is `u32`; values larger than `u32::MAX` are
    /// rejected instead of silently truncating.
    pub fn from_row_major_u64(
        n_genes: usize,
        n_samples: usize,
        counts: Vec<u64>,
    ) -> Result<Self, DeseqError> {
        let converted = counts
            .into_iter()
            .enumerate()
            .map(|(idx, value)| {
                u32::try_from(value).map_err(|_| DeseqError::InvalidCounts {
                    reason: format!("count at index {idx} exceeds u32::MAX"),
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_row_major_u32(n_genes, n_samples, converted)
    }

    /// Build a count matrix with optional row and column names.
    pub fn from_row_major_u32_with_names(
        n_genes: usize,
        n_samples: usize,
        counts: Vec<u32>,
        gene_names: Option<Vec<String>>,
        sample_names: Option<Vec<String>>,
    ) -> Result<Self, DeseqError> {
        if n_genes == 0 {
            return Err(DeseqError::InvalidCounts {
                reason: "count matrix must contain at least one gene".to_string(),
            });
        }
        if n_samples == 0 {
            return Err(DeseqError::InvalidCounts {
                reason: "count matrix must contain at least one sample".to_string(),
            });
        }
        let expected =
            n_genes
                .checked_mul(n_samples)
                .ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "count matrix shape overflow".to_string(),
                    expected: usize::MAX,
                    actual: counts.len(),
                })?;
        if counts.len() != expected {
            return Err(invalid_dimensions(
                "count matrix values",
                expected,
                counts.len(),
            ));
        }
        if let Some(names) = &gene_names {
            if names.len() != n_genes {
                return Err(invalid_dimensions("gene names", n_genes, names.len()));
            }
        }
        if let Some(names) = &sample_names {
            if names.len() != n_samples {
                return Err(invalid_dimensions("sample names", n_samples, names.len()));
            }
        }
        Ok(Self {
            n_genes,
            n_samples,
            counts,
            gene_names,
            sample_names,
        })
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        self.n_genes
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Count values in row-major order.
    pub fn as_slice(&self) -> &[u32] {
        &self.counts
    }

    /// Return a gene row.
    pub fn row(&self, gene: usize) -> Result<&[u32], DeseqError> {
        if gene >= self.n_genes {
            return Err(DeseqError::InvalidDimensions {
                context: "gene index".to_string(),
                expected: self.n_genes.saturating_sub(1),
                actual: gene,
            });
        }
        Ok(self.row_values(gene))
    }

    pub(crate) fn row_values(&self, gene: usize) -> &[u32] {
        let start = gene * self.n_samples;
        &self.counts[start..start + self.n_samples]
    }

    /// Whether all samples are zero for a gene.
    pub fn is_all_zero_gene(&self, gene: usize) -> Result<bool, DeseqError> {
        Ok(self.row(gene)?.iter().all(|count| *count == 0))
    }

    /// Per-gene all-zero flags.
    pub fn all_zero_flags(&self) -> Vec<bool> {
        (0..self.n_genes)
            .map(|gene| self.row_values(gene).iter().all(|count| *count == 0))
            .collect()
    }

    /// Optional gene names.
    pub fn gene_names(&self) -> Option<&[String]> {
        self.gene_names.as_deref()
    }

    /// Optional sample names.
    pub fn sample_names(&self) -> Option<&[String]> {
        self.sample_names.as_deref()
    }

    /// Return a basic summary used by fit-state output.
    pub fn summary(&self) -> CountsSummary {
        let all_zero_genes = (0..self.n_genes)
            .filter(|gene| self.row_values(*gene).iter().all(|count| *count == 0))
            .count();
        CountsSummary {
            n_genes: self.n_genes,
            n_samples: self.n_samples,
            all_zero_genes,
        }
    }

    /// Select gene rows in caller-provided order while preserving sample names.
    pub fn select_rows(&self, gene_indices: &[usize]) -> Result<Self, DeseqError> {
        if gene_indices.is_empty() {
            return Err(DeseqError::InvalidDimensions {
                context: "selected count rows".to_string(),
                expected: 1,
                actual: 0,
            });
        }
        let mut values = Vec::with_capacity(gene_indices.len() * self.n_samples);
        for gene in gene_indices {
            values.extend_from_slice(self.row(*gene)?);
        }
        let gene_names = self.gene_names().map(|names| {
            gene_indices
                .iter()
                .map(|gene| names[*gene].clone())
                .collect::<Vec<_>>()
        });
        let sample_names = self.sample_names().map(<[String]>::to_vec);
        CountMatrix::from_row_major_u32_with_names(
            gene_indices.len(),
            self.n_samples,
            values,
            gene_names,
            sample_names,
        )
    }
}

/// Basic count matrix summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CountsSummary {
    /// Number of genes.
    pub n_genes: usize,
    /// Number of samples.
    pub n_samples: usize,
    /// Number of all-zero genes.
    pub all_zero_genes: usize,
}

/// Inspectable fit state for all implemented and future DESeq2 stages.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqFit {
    /// Count matrix summary.
    pub counts_summary: CountsSummary,
    /// Optional design matrix information.
    pub design: Option<DesignMatrix>,
    /// Optional reduced design matrix information for likelihood-ratio tests.
    pub reduced_design: Option<DesignMatrix>,
    /// Estimated or supplied size factors.
    pub size_factors: Vec<f64>,
    /// Optional gene/sample normalization factors.
    pub normalization_factors: Option<RowMajorMatrix<f64>>,
    /// Optional DESeq2-style observation weights used by builder stages.
    ///
    /// When a design is available these are row-normalized weights from
    /// `getAndCheckWeights`-style preprocessing. `baseMean` and `baseVar` still
    /// use the raw caller-supplied weights first, matching DESeq2's
    /// `getBaseMeansAndVariances` ordering.
    pub observation_weights: Option<RowMajorMatrix<f64>>,
    /// Rows that failed DESeq2-style observation-weight design checks.
    pub weights_fail: Option<Vec<bool>>,
    /// Rank of the unweighted design during observation-weight preprocessing.
    pub weights_design_rank: Option<usize>,
    /// Per-gene base means.
    pub base_mean: Vec<f64>,
    /// Per-gene sample variance of normalized counts, matching DESeq2 `baseVar`.
    pub base_var: Vec<f64>,
    /// Per-gene all-zero flags, matching DESeq2 `allZero`.
    ///
    /// For design-aware observation-weight stages, rows with `weightsFail` are
    /// also marked true, matching DESeq2's downstream skip behavior.
    pub all_zero: Vec<bool>,
    /// Gene-wise dispersion estimates.
    pub disp_gene_est: Option<Vec<f64>>,
    /// Gene-wise dispersion iteration counts, matching DESeq2 `dispGeneIter`.
    pub disp_gene_iter: Option<Vec<usize>>,
    /// Fitted dispersion trend values.
    pub disp_fit: Option<Vec<f64>>,
    /// Fitted dispersion trend object used by implemented VST dispatch.
    pub dispersion_trend: Option<DispersionTrendFit>,
    /// MAP dispersion estimates before outlier replacement.
    pub disp_map: Option<Vec<f64>>,
    /// Final dispersion estimates.
    pub dispersion: Option<Vec<f64>>,
    /// MAP dispersion iteration counts.
    pub disp_iter: Option<Vec<usize>>,
    /// MAP dispersion outlier flags.
    pub disp_outlier: Option<Vec<bool>>,
    /// Dispersion prior variance.
    pub disp_prior_var: Option<f64>,
    /// Dispersion convergence flags.
    pub dispersion_converged: Option<Vec<bool>>,
    /// GLM beta estimates.
    pub beta: Option<RowMajorMatrix<f64>>,
    /// GLM beta standard errors.
    pub beta_se: Option<RowMajorMatrix<f64>>,
    /// Per-gene GLM beta covariance matrices on log2 scale.
    ///
    /// Stored as genes x `(n_coefficients * n_coefficients)`, with each gene
    /// row containing a row-major coefficient covariance matrix.
    pub beta_covariance: Option<RowMajorMatrix<f64>>,
    /// GLM beta convergence flags.
    pub beta_converged: Option<Vec<bool>>,
    /// GLM beta iteration counts.
    pub beta_iter: Option<Vec<usize>>,
    /// Per-gene fitted-model log likelihoods from the full GLM.
    pub log_like: Option<Vec<f64>>,
    /// Per-gene full-model deviance, matching DESeq2's `-2 * logLike` field.
    pub full_deviance: Option<Vec<f64>>,
    /// Per-gene reduced-model log likelihoods for LRT pipelines.
    pub reduced_log_like: Option<Vec<f64>>,
    /// Per-gene reduced-model beta convergence flags for LRT pipelines.
    pub reduced_beta_converged: Option<Vec<bool>>,
    /// Per-gene reduced-model beta iteration counts for LRT pipelines.
    pub reduced_beta_iter: Option<Vec<usize>>,
    /// Reduced-model fitted mean matrix for LRT pipelines.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub reduced_mu: Option<RowMajorMatrix<f64>>,
    /// Reduced-model hat diagonal matrix for LRT pipelines.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub reduced_hat_diagonal: Option<RowMajorMatrix<f64>>,
    /// Fitted mean matrix.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub mu: Option<RowMajorMatrix<f64>>,
    /// Cook's distance matrix.
    pub cooks: Option<RowMajorMatrix<f64>>,
    /// Per-gene maximum Cook's distance over eligible samples.
    pub max_cooks: Option<Vec<Option<f64>>>,
    /// Hat diagonal matrix.
    ///
    /// Kept as matrix-valued fit state rather than `mcols(dds)`-style row
    /// metadata.
    pub hat_diagonal: Option<RowMajorMatrix<f64>>,
    /// Wald output.
    pub wald: Option<WaldOutput>,
    /// LRT output.
    pub lrt: Option<LrtOutput>,
}

/// Builder for implemented and future DESeq2 workflow stages.
#[derive(Clone, Debug)]
pub struct DeseqBuilder {
    fit_type: FitType,
    test: TestType,
    size_factor_options: SizeFactorOptions,
    normalization_factors: Option<RowMajorMatrix<f64>>,
    observation_weights: Option<RowMajorMatrix<f64>>,
    observation_weight_options: ObservationWeightOptions,
    execution_mode: ExecutionMode,
    threads: Option<usize>,
    reduced_design: Option<DesignMatrix>,
    irls_options: IrlsOptions,
    gene_wise_dispersion_options: GeneWiseDispersionOptions,
    wald_test_options: WaldTestOptions,
    cooks_cutoff: CooksCutoff,
    independent_filtering_options: IndependentFilteringOptions,
}

impl DeseqFit {
    /// Reconstruct DESeq2-style normalized counts for the count matrix used to create this fit.
    pub fn normalized_counts(
        &self,
        counts: &CountMatrix,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "normalized counts")?;
        match &self.normalization_factors {
            Some(factors) => normalized_counts_with_factors(counts, factors),
            None => normalized_counts(counts, &self.size_factors),
        }
    }

    /// Apply DESeq2's `normTransform` to the count matrix used to create this fit.
    pub fn norm_transform(&self, counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
        let normalized = self.normalized_counts(counts)?;
        norm_transform(&normalized)
    }

    /// Apply VST using this fit's stored dispersion trend.
    ///
    /// This is a lower-level analogue of DESeq2's
    /// `getVarianceStabilizedData`: the fit must already include a fitted
    /// dispersion trend from one of the implemented trend stages.
    pub fn variance_stabilizing_transform(
        &self,
        counts: &CountMatrix,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "variance-stabilizing transform")?;
        let trend_fit =
            self.dispersion_trend
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "a fitted dispersion trend is required before VST".to_string(),
                })?;
        self.variance_stabilizing_transform_with_trend(counts, trend_fit)
    }

    /// Apply VST using a caller-supplied fitted dispersion trend.
    ///
    /// This is useful for DESeq2's fast `vst()` shape, where normalization and
    /// full-count reconstruction come from the full fit but the trend may have
    /// been estimated on a deterministic subset.
    pub fn variance_stabilizing_transform_with_trend(
        &self,
        counts: &CountMatrix,
        trend_fit: &DispersionTrendFit,
    ) -> Result<RowMajorMatrix<f64>, DeseqError> {
        validate_fit_counts_shape(self, counts, "variance-stabilizing transform")?;
        let normalized = self.normalized_counts(counts)?;
        match &self.normalization_factors {
            Some(factors) => {
                vst_with_dispersion_trend_and_normalization_factors(&normalized, trend_fit, factors)
            }
            None => vst_with_dispersion_trend_and_size_factors(
                &normalized,
                trend_fit,
                &self.size_factors,
            ),
        }
    }

    /// Short alias for [`DeseqFit::variance_stabilizing_transform`].
    pub fn vst(&self, counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
        self.variance_stabilizing_transform(counts)
    }

    /// Apply the implemented regularized-log sample-effect transform.
    ///
    /// The fit must already contain fitted trend dispersions (`dispFit`) and
    /// final gene-wise dispersions (`dispersion`). The rlog sample-effect prior
    /// is estimated from this fit's normalized counts, `baseMean`, and
    /// `dispFit`, then the sample-effect ridge GLM is fit with this fit's
    /// size factors or normalization factors.
    pub fn regularized_log_transform(
        &self,
        counts: &CountMatrix,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_options(counts, IrlsOptions::default())
    }

    /// Apply the implemented regularized-log transform with explicit IRLS options.
    pub fn regularized_log_transform_with_options(
        &self,
        counts: &CountMatrix,
        options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        validate_fit_counts_shape(self, counts, "regularized-log transform")?;
        if self.all_zero.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog allZero rows",
                counts.n_genes(),
                self.all_zero.len(),
            ));
        }
        let disp_fit = self
            .disp_fit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted dispersion trend values are required before rlog".to_string(),
            })?;
        let dispersions =
            self.dispersion
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "final dispersions are required before rlog".to_string(),
                })?;
        let nonzero_rows = nonzero_gene_indices(&self.all_zero);
        if nonzero_rows.len() != counts.n_genes() {
            if nonzero_rows.is_empty() {
                let transformed =
                    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 0.0)?;
                let offset_mode = match &self.normalization_factors {
                    Some(_) => RlogOffsetMode::NormalizationFactors,
                    None => RlogOffsetMode::SizeFactors,
                };
                return Ok(RlogOutput {
                    transformed,
                    intercept: vec![0.0; counts.n_genes()],
                    sample_prior_variance: 0.0,
                    offset_mode,
                });
            }
            let compact_counts = compact_counts(counts, &nonzero_rows)?;
            let compact_base_mean = compact_f64_values(&self.base_mean, &nonzero_rows)?;
            let compact_disp_fit = compact_f64_values(disp_fit, &nonzero_rows)?;
            let compact_dispersions = compact_f64_values(dispersions, &nonzero_rows)?;
            let compact_output = match &self.normalization_factors {
                Some(factors) => {
                    let compact_factors = compact_matrix_rows(factors, &nonzero_rows)?;
                    rlog_with_estimated_prior_and_normalization_factors(
                        &compact_counts,
                        &compact_factors,
                        &compact_base_mean,
                        &compact_disp_fit,
                        &compact_dispersions,
                        options,
                    )?
                }
                None => rlog_with_estimated_prior_and_size_factors(
                    &compact_counts,
                    &self.size_factors,
                    &compact_base_mean,
                    &compact_disp_fit,
                    &compact_dispersions,
                    options,
                )?,
            };
            return expand_rlog_output_with_all_zero_rows(
                compact_output,
                &self.all_zero,
                counts.n_samples(),
            );
        }
        match &self.normalization_factors {
            Some(factors) => rlog_with_estimated_prior_and_normalization_factors(
                counts,
                factors,
                &self.base_mean,
                disp_fit,
                dispersions,
                options,
            ),
            None => rlog_with_estimated_prior_and_size_factors(
                counts,
                &self.size_factors,
                &self.base_mean,
                disp_fit,
                dispersions,
                options,
            ),
        }
    }

    /// Short alias for [`DeseqFit::regularized_log_transform`].
    pub fn rlog(&self, counts: &CountMatrix) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform(counts)
    }

    /// Apply rlog with supplied frozen intercepts and sample-effect prior variance.
    ///
    /// This uses the fit state's final dispersions and size-factor or
    /// normalization-factor offsets, while treating `frozen_intercept` as the
    /// fixed log2 intercept surface. It is the fit-level building block for
    /// frozen rlog reuse.
    pub fn regularized_log_transform_with_frozen_intercept(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_frozen_intercept_and_options(
            counts,
            frozen_intercept,
            sample_prior_variance,
            IrlsOptions::default(),
        )
    }

    /// Apply frozen-intercept rlog with explicit IRLS options.
    pub fn regularized_log_transform_with_frozen_intercept_and_options(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
        options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        validate_fit_counts_shape(self, counts, "frozen regularized-log transform")?;
        if frozen_intercept.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog frozen intercept rows",
                counts.n_genes(),
                frozen_intercept.len(),
            ));
        }
        if self.all_zero.len() != counts.n_genes() {
            return Err(invalid_dimensions(
                "rlog allZero rows",
                counts.n_genes(),
                self.all_zero.len(),
            ));
        }
        let dispersions =
            self.dispersion
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "final dispersions are required before frozen rlog".to_string(),
                })?;
        let nonzero_rows = nonzero_gene_indices(&self.all_zero);
        if nonzero_rows.len() != counts.n_genes() {
            if nonzero_rows.is_empty() {
                let transformed =
                    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 0.0)?;
                let offset_mode = match &self.normalization_factors {
                    Some(_) => RlogOffsetMode::NormalizationFactors,
                    None => RlogOffsetMode::SizeFactors,
                };
                return Ok(RlogOutput {
                    transformed,
                    intercept: frozen_intercept.to_vec(),
                    sample_prior_variance,
                    offset_mode,
                });
            }
            let compact_counts = compact_counts(counts, &nonzero_rows)?;
            let compact_dispersions = compact_f64_values(dispersions, &nonzero_rows)?;
            let compact_intercept = compact_f64_values(frozen_intercept, &nonzero_rows)?;
            let compact_output = match &self.normalization_factors {
                Some(factors) => {
                    let compact_factors = compact_matrix_rows(factors, &nonzero_rows)?;
                    rlog_frozen_with_normalization_factors(
                        &compact_counts,
                        &compact_factors,
                        &compact_dispersions,
                        sample_prior_variance,
                        &compact_intercept,
                        options,
                    )?
                }
                None => rlog_frozen_with_size_factors(
                    &compact_counts,
                    &self.size_factors,
                    &compact_dispersions,
                    sample_prior_variance,
                    &compact_intercept,
                    options,
                )?,
            };
            return expand_frozen_rlog_output_with_all_zero_rows(
                compact_output,
                &self.all_zero,
                counts.n_samples(),
                frozen_intercept,
            );
        }
        match &self.normalization_factors {
            Some(factors) => rlog_frozen_with_normalization_factors(
                counts,
                factors,
                dispersions,
                sample_prior_variance,
                frozen_intercept,
                options,
            ),
            None => rlog_frozen_with_size_factors(
                counts,
                &self.size_factors,
                dispersions,
                sample_prior_variance,
                frozen_intercept,
                options,
            ),
        }
    }

    /// Short alias for [`DeseqFit::regularized_log_transform_with_frozen_intercept`].
    pub fn frozen_rlog(
        &self,
        counts: &CountMatrix,
        frozen_intercept: &[f64],
        sample_prior_variance: f64,
    ) -> Result<RlogOutput, DeseqError> {
        self.regularized_log_transform_with_frozen_intercept(
            counts,
            frozen_intercept,
            sample_prior_variance,
        )
    }

    /// Number of rows eligible for DESeq2's fast `vst()` subset.
    ///
    /// This uses the stored `baseMean` vector and the DESeq2 `baseMean > 5`
    /// rule, with the same finite-value validation as subset construction.
    pub fn fast_vst_eligible_count(&self) -> Result<usize, DeseqError> {
        count_fast_vst_eligible_rows(&self.base_mean)
    }

    /// Build the row-aligned subset used by DESeq2's fast `vst()` trend fit.
    ///
    /// The returned bundle includes raw counts, normalized counts, optional
    /// normalization factors, optional observation weights, and original row
    /// indices, all selected from this fit's `baseMean` vector.
    pub fn fast_vst_subset(
        &self,
        counts: &CountMatrix,
        nsub: usize,
    ) -> Result<FastVstSubset, DeseqError> {
        validate_fit_counts_shape(self, counts, "fast VST subset")?;
        let normalized = self.normalized_counts(counts)?;
        build_fast_vst_subset(
            counts,
            &normalized,
            &self.base_mean,
            nsub,
            self.normalization_factors.as_ref(),
            self.observation_weights.as_ref(),
        )
    }
}

fn validate_fit_counts_shape(
    fit: &DeseqFit,
    counts: &CountMatrix,
    context: &str,
) -> Result<(), DeseqError> {
    let expected = fit.counts_summary.n_genes * fit.counts_summary.n_samples;
    let actual = counts.n_genes() * counts.n_samples();
    if fit.counts_summary.n_genes != counts.n_genes()
        || fit.counts_summary.n_samples != counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected,
            actual,
        });
    }
    Ok(())
}

struct NormalizationStages {
    size_factors: Vec<f64>,
    base_mean: Vec<f64>,
    base_var: Vec<f64>,
    all_zero: Vec<bool>,
    normalized: RowMajorMatrix<f64>,
    normalization_factors: Option<RowMajorMatrix<f64>>,
    observation_weights: Option<RowMajorMatrix<f64>>,
    weights_fail: Option<Vec<bool>>,
    weights_design_rank: Option<usize>,
}

impl NormalizationStages {
    fn into_base_fit_input(self) -> BaseFitInput {
        BaseFitInput {
            size_factors: self.size_factors,
            normalization_factors: self.normalization_factors,
            observation_weights: self.observation_weights,
            weights_fail: self.weights_fail,
            weights_design_rank: self.weights_design_rank,
            base_mean: self.base_mean,
            base_var: self.base_var,
            all_zero: self.all_zero,
        }
    }
}

struct BaseFitInput {
    size_factors: Vec<f64>,
    normalization_factors: Option<RowMajorMatrix<f64>>,
    observation_weights: Option<RowMajorMatrix<f64>>,
    weights_fail: Option<Vec<bool>>,
    weights_design_rank: Option<usize>,
    base_mean: Vec<f64>,
    base_var: Vec<f64>,
    all_zero: Vec<bool>,
}

struct WeightedBaseMetadata {
    base_mean: Vec<f64>,
    base_var: Vec<f64>,
    observation_weights: Option<RowMajorMatrix<f64>>,
    weights_fail: Option<Vec<bool>>,
    weights_design_rank: Option<usize>,
}

struct WaldPipelineInput<'a> {
    counts: &'a CountMatrix,
    design: &'a DesignMatrix,
    size_factors: &'a [f64],
    normalization_factors: Option<&'a RowMajorMatrix<f64>>,
    observation_weights: Option<&'a RowMajorMatrix<f64>>,
    normalized: &'a RowMajorMatrix<f64>,
    base_mean: &'a [f64],
    all_zero: &'a [bool],
    dispersions: &'a [f64],
    coefficient: usize,
}

struct LrtPipelineInput<'a> {
    counts: &'a CountMatrix,
    full_design: &'a DesignMatrix,
    reduced_design: &'a DesignMatrix,
    size_factors: &'a [f64],
    normalization_factors: Option<&'a RowMajorMatrix<f64>>,
    observation_weights: Option<&'a RowMajorMatrix<f64>>,
    normalized: &'a RowMajorMatrix<f64>,
    base_mean: &'a [f64],
    all_zero: &'a [bool],
    dispersions: &'a [f64],
    coefficient: usize,
}

#[derive(Clone, Copy)]
struct FixedDispersionGlmInput<'a> {
    counts: &'a CountMatrix,
    design: &'a DesignMatrix,
    size_factors: &'a [f64],
    normalization_factors: Option<&'a RowMajorMatrix<f64>>,
    observation_weights: Option<&'a RowMajorMatrix<f64>>,
    all_zero: &'a [bool],
    dispersions: &'a [f64],
}

struct FixedDispersionGlmOutput {
    glm_fit: NbinomGlmFit,
    expanded_dispersions: Vec<f64>,
}

struct WaldPipelineOutput {
    glm_fit: NbinomGlmFit,
    wald: WaldOutput,
    cooks: CooksOutput,
    results: DeseqResults,
    expanded_dispersions: Vec<f64>,
}

struct WaldContrastPipelineOutput {
    glm_fit: NbinomGlmFit,
    wald_contrast: WaldContrastOutput,
    cooks: CooksOutput,
    results: DeseqResults,
    expanded_dispersions: Vec<f64>,
}

struct LrtPipelineOutput {
    full_fit: NbinomGlmFit,
    reduced_fit: NbinomGlmFit,
    lrt: LrtOutput,
    cooks: CooksOutput,
    results: DeseqResults,
    expanded_dispersions: Vec<f64>,
}

/// Output from the limited native Wald replacement-refit path.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksReplacementWaldOutput {
    /// Original fit on the caller-supplied counts, before replacement.
    pub original_fit: DeseqFit,
    /// Original result rows before replacement/refit.
    pub original_results: DeseqResults,
    /// Replacement/refit planning metadata.
    pub refit_plan: CooksRefitPlan,
    /// Refit on replacement counts, when any non-all-zero replacement row exists.
    pub refit_fit: Option<DeseqFit>,
    /// Result rows from the replacement-count refit, before merge.
    pub refit_results: Option<DeseqResults>,
    /// Final merged result rows after replacing only `refitRows`.
    pub results: DeseqResults,
}

/// Output from the limited native LRT replacement-refit path.
#[derive(Clone, Debug, PartialEq)]
pub struct CooksReplacementLrtOutput {
    /// Original fit on the caller-supplied counts, before replacement.
    pub original_fit: DeseqFit,
    /// Original result rows before replacement/refit.
    pub original_results: DeseqResults,
    /// Replacement/refit planning metadata.
    pub refit_plan: CooksRefitPlan,
    /// Refit on replacement counts, when any non-all-zero replacement row exists.
    pub refit_fit: Option<DeseqFit>,
    /// Result rows from the replacement-count refit, before merge.
    pub refit_results: Option<DeseqResults>,
    /// Final merged result rows after replacing only `refitRows`.
    pub results: DeseqResults,
}

/// Output from a top-level Cook's replacement-refit workflow selected by `test`.
#[derive(Clone, Debug, PartialEq)]
pub enum CooksReplacementTestOutput {
    /// Wald replacement-refit output.
    Wald(CooksReplacementWaldOutput),
    /// LRT replacement-refit output.
    Lrt(CooksReplacementLrtOutput),
}

impl CooksReplacementTestOutput {
    /// Final merged result rows after replacement/refit.
    pub fn results(&self) -> &DeseqResults {
        match self {
            Self::Wald(output) => &output.results,
            Self::Lrt(output) => &output.results,
        }
    }

    /// Original fit before replacement/refit.
    pub fn original_fit(&self) -> &DeseqFit {
        match self {
            Self::Wald(output) => &output.original_fit,
            Self::Lrt(output) => &output.original_fit,
        }
    }

    /// Replacement/refit planning metadata.
    pub fn refit_plan(&self) -> &CooksRefitPlan {
        match self {
            Self::Wald(output) => &output.refit_plan,
            Self::Lrt(output) => &output.refit_plan,
        }
    }
}

/// Output from the current fast-VST GLM-mu helper.
#[derive(Clone, Debug, PartialEq)]
pub struct FastVstGlmMuOutput {
    /// Full count matrix transformed with the subset-fitted dispersion trend.
    pub transformed: RowMajorMatrix<f64>,
    /// Fit state for the deterministic fast-VST subset used to estimate the trend.
    pub subset_fit: DeseqFit,
    /// Row-aligned subset bundle with original row indices and optional factors.
    pub subset: FastVstSubset,
}

/// Metadata summary for the explicit fast-VST GLM-mu helper.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastVstGlmMuMetadata {
    /// Number of rows in the transformed full matrix.
    pub transformed_rows: usize,
    /// Number of columns in the transformed full matrix.
    pub transformed_cols: usize,
    /// Number of rows in the deterministic fast subset.
    pub fast_subset_rows: usize,
    /// Number of samples in the deterministic fast subset.
    pub fast_subset_cols: usize,
    /// Original zero-based row indices selected for the fast subset.
    pub fast_subset_indices: Vec<usize>,
    /// Number of rows used to fit the subset trend.
    pub trend_fit_rows: usize,
    /// Number of samples in the subset trend fit.
    pub trend_fit_cols: usize,
    /// Stable fit-type label for the fitted dispersion trend.
    pub trend_fit_type: Option<&'static str>,
}

/// Output from the automatic GLM-mu VST helper.
#[derive(Clone, Debug, PartialEq)]
pub struct VstGlmMuOutput {
    /// Full count matrix after variance stabilization.
    pub transformed: RowMajorMatrix<f64>,
    /// Fit state that supplied the dispersion trend.
    pub trend_fit: DeseqFit,
    /// Source of the fitted dispersion trend used for this transform.
    pub trend_source: VstTrendSource,
    /// Fast-VST subset diagnostics when the fast path was used.
    pub fast_subset: Option<FastVstSubset>,
}

/// Output from the GLM-mu rlog builder helper with retained fit state.
#[derive(Clone, Debug, PartialEq)]
pub struct RlogGlmMuOutput {
    /// Regularized-log output matrix and prior metadata.
    pub rlog: RlogOutput,
    /// Fit state that supplied `baseMean`, `dispFit`, and final dispersions.
    pub fit: DeseqFit,
    /// Stable design mode used by the builder helper.
    pub design_mode: RlogDesignMode,
}

/// Output from a builder-level frozen-rlog reuse workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenRlogGlmMuOutput {
    /// Initial rlog fit that supplied frozen intercepts and prior variance.
    pub source_rlog: RlogOutput,
    /// Frozen-intercept rlog transform fit from the same dispersion state.
    pub frozen_rlog: RlogOutput,
    /// Fit state that supplied final dispersions and offsets.
    pub fit: DeseqFit,
    /// Stable design mode used by the builder helper.
    pub design_mode: RlogDesignMode,
}

/// Metadata summary for the GLM-mu rlog builder helper.
#[derive(Clone, Debug, PartialEq)]
pub struct RlogGlmMuMetadata {
    /// Metadata from the rlog transform itself.
    pub rlog: RlogMetadata,
    /// Stable design mode label.
    pub design_mode: &'static str,
    /// Number of rows in the fit state that supplied dispersion inputs.
    pub fit_rows: usize,
    /// Number of samples in the fit state that supplied dispersion inputs.
    pub fit_cols: usize,
    /// Stable fit-type label for the fitted dispersion trend.
    pub trend_fit_type: Option<&'static str>,
}

/// Metadata summary for a builder-level frozen-rlog reuse workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenRlogGlmMuMetadata {
    /// Metadata from the initial rlog fit.
    pub source_rlog: RlogMetadata,
    /// Metadata from the frozen-intercept transform.
    pub frozen_rlog: RlogMetadata,
    /// Stable design mode label.
    pub design_mode: &'static str,
    /// Number of rows in the fit state that supplied dispersion inputs.
    pub fit_rows: usize,
    /// Number of samples in the fit state that supplied dispersion inputs.
    pub fit_cols: usize,
    /// Stable fit-type label for the fitted dispersion trend.
    pub trend_fit_type: Option<&'static str>,
}

/// Design mode used by a GLM-mu rlog builder helper.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RlogDesignMode {
    /// Caller-supplied design-aware fit.
    DesignAware,
    /// Intercept-only blind fit.
    Blind,
}

/// Metadata summary for the automatic GLM-mu VST helper.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VstGlmMuMetadata {
    /// Stable source label for the fitted trend.
    pub trend_source: &'static str,
    /// Requested fast-subset size considered by automatic VST.
    pub nsub: usize,
    /// Number of rows passing the `baseMean > 5` fast-VST eligibility rule.
    pub eligible_rows: usize,
    /// Whether the deterministic fast subset supplied the fitted trend.
    pub used_fast_subset: bool,
    /// Stable reason label when the full-data trend path was selected.
    pub full_data_reason: Option<&'static str>,
    /// Number of rows in the transformed full matrix.
    pub transformed_rows: usize,
    /// Number of columns in the transformed full matrix.
    pub transformed_cols: usize,
    /// Number of rows used to fit the trend.
    pub trend_fit_rows: usize,
    /// Number of samples in the trend fit.
    pub trend_fit_cols: usize,
    /// Stable fit-type label for the fitted dispersion trend.
    pub trend_fit_type: Option<&'static str>,
    /// Number of rows in the fast subset, when that path was used.
    pub fast_subset_rows: Option<usize>,
    /// Original zero-based row indices in the fast subset, when that path was used.
    pub fast_subset_indices: Option<Vec<usize>>,
}

/// Source of the fitted dispersion trend used by automatic VST.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VstTrendSource {
    /// The trend was fit on DESeq2's deterministic fast-VST subset.
    FastSubset {
        /// Requested fast-subset size.
        nsub: usize,
        /// Number of rows passing the `baseMean > 5` eligibility rule.
        eligible_rows: usize,
    },
    /// The trend was fit on the full count matrix.
    FullData {
        /// Requested fast-subset size that could not be satisfied.
        nsub: usize,
        /// Number of rows passing the `baseMean > 5` eligibility rule.
        eligible_rows: usize,
        /// Reason the full-data trend path was selected.
        reason: VstFullDataReason,
    },
}

/// Reason automatic VST selected the full-data trend path.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VstFullDataReason {
    /// Fewer than `nsub` rows passed the fast-VST eligibility rule.
    InsufficientEligibleRows,
}

impl FastVstGlmMuOutput {
    /// DESeq2-shaped metadata view for explicit fast-VST diagnostics.
    pub fn metadata(&self) -> FastVstGlmMuMetadata {
        FastVstGlmMuMetadata {
            transformed_rows: self.transformed.n_rows(),
            transformed_cols: self.transformed.n_cols(),
            fast_subset_rows: self.subset.counts.n_genes(),
            fast_subset_cols: self.subset.counts.n_samples(),
            fast_subset_indices: self.subset.row_indices.clone(),
            trend_fit_rows: self.subset_fit.counts_summary.n_genes,
            trend_fit_cols: self.subset_fit.counts_summary.n_samples,
            trend_fit_type: self
                .subset_fit
                .dispersion_trend
                .as_ref()
                .map(|trend| trend.fit_type_label()),
        }
    }
}

impl VstGlmMuOutput {
    /// DESeq2-shaped metadata view for automatic VST diagnostics.
    pub fn metadata(&self) -> VstGlmMuMetadata {
        VstGlmMuMetadata {
            trend_source: self.trend_source.as_str(),
            nsub: self.trend_source.nsub(),
            eligible_rows: self.trend_source.eligible_rows(),
            used_fast_subset: self.trend_source.used_fast_subset(),
            full_data_reason: self
                .trend_source
                .full_data_reason()
                .map(|reason| reason.as_str()),
            transformed_rows: self.transformed.n_rows(),
            transformed_cols: self.transformed.n_cols(),
            trend_fit_rows: self.trend_fit.counts_summary.n_genes,
            trend_fit_cols: self.trend_fit.counts_summary.n_samples,
            trend_fit_type: self
                .trend_fit
                .dispersion_trend
                .as_ref()
                .map(|trend| trend.fit_type_label()),
            fast_subset_rows: self
                .fast_subset
                .as_ref()
                .map(|subset| subset.counts.n_genes()),
            fast_subset_indices: self
                .fast_subset
                .as_ref()
                .map(|subset| subset.row_indices.clone()),
        }
    }
}

impl RlogGlmMuOutput {
    /// Metadata view for wrappers, diagnostics, and benchmark logs.
    pub fn metadata(&self) -> RlogGlmMuMetadata {
        RlogGlmMuMetadata {
            rlog: self.rlog.metadata(),
            design_mode: self.design_mode.as_str(),
            fit_rows: self.fit.counts_summary.n_genes,
            fit_cols: self.fit.counts_summary.n_samples,
            trend_fit_type: self
                .fit
                .dispersion_trend
                .as_ref()
                .map(|trend| trend.fit_type_label()),
        }
    }
}

impl FrozenRlogGlmMuOutput {
    /// Metadata view for wrappers, diagnostics, and benchmark logs.
    pub fn metadata(&self) -> FrozenRlogGlmMuMetadata {
        FrozenRlogGlmMuMetadata {
            source_rlog: self.source_rlog.metadata(),
            frozen_rlog: self.frozen_rlog.metadata(),
            design_mode: self.design_mode.as_str(),
            fit_rows: self.fit.counts_summary.n_genes,
            fit_cols: self.fit.counts_summary.n_samples,
            trend_fit_type: self
                .fit
                .dispersion_trend
                .as_ref()
                .map(|trend| trend.fit_type_label()),
        }
    }
}

impl RlogDesignMode {
    /// Stable label for wrappers and benchmark logs.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::DesignAware => "designAware",
            Self::Blind => "blind",
        }
    }
}

impl VstTrendSource {
    /// Stable DESeq2-shaped label for the automatic VST trend source.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::FastSubset { .. } => "fastSubset",
            Self::FullData { .. } => "fullData",
        }
    }

    /// Requested fast-subset size considered by the automatic VST helper.
    pub fn nsub(&self) -> usize {
        match self {
            Self::FastSubset { nsub, .. } | Self::FullData { nsub, .. } => *nsub,
        }
    }

    /// Number of rows passing the `baseMean > 5` fast-VST eligibility rule.
    pub fn eligible_rows(&self) -> usize {
        match self {
            Self::FastSubset { eligible_rows, .. } | Self::FullData { eligible_rows, .. } => {
                *eligible_rows
            }
        }
    }

    /// Whether automatic VST fit the trend on the deterministic fast subset.
    pub fn used_fast_subset(&self) -> bool {
        matches!(self, Self::FastSubset { .. })
    }

    /// Reason the full-data trend path was selected, if applicable.
    pub fn full_data_reason(&self) -> Option<VstFullDataReason> {
        match self {
            Self::FastSubset { .. } => None,
            Self::FullData { reason, .. } => Some(*reason),
        }
    }
}

impl VstFullDataReason {
    /// Stable label for why automatic VST selected the full-data trend path.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::InsufficientEligibleRows => "insufficientEligibleRows",
        }
    }
}

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

    /// Run only the implemented initial normalization stages.
    pub fn fit_size_factors_and_base_means(
        &self,
        counts: &CountMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let stages = self.normalization_stages(counts)?;
        Ok(Self::base_fit(counts, None, stages.into_base_fit_input()))
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
        let stages = self.normalization_stages_for_design(counts, design)?;
        Ok(Self::base_fit(
            counts,
            Some(design.clone()),
            stages.into_base_fit_input(),
        ))
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
        self.ensure_no_observation_weights("native linear-mu dispersion estimation")?;
        design.validate_full_rank("linear-mu dispersion")?;
        let stages = self.normalization_stages(counts)?;
        let dispersion = estimate_gene_wise_dispersions_linear_mu(
            GeneWiseDispersionInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                normalized_counts: &stages.normalized,
                base_mean: &stages.base_mean,
                base_var: &stages.base_var,
                all_zero: &stages.all_zero,
                observation_weights: None,
            },
            self.gene_wise_dispersion_options,
        )?;
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
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
        design.validate_full_rank("GLM-mu dispersion")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let dispersion = estimate_gene_wise_dispersions_glm_mu(
            GeneWiseDispersionInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                normalized_counts: &stages.normalized,
                base_mean: &stages.base_mean,
                base_var: &stages.base_var,
                all_zero: &stages.all_zero,
                observation_weights: stages.observation_weights.as_ref(),
            },
            self.gene_wise_dispersion_options,
            self.irls_options.clone(),
        )?;
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.disp_gene_est = Some(dispersion.disp_gene_est.clone());
        fit.disp_gene_iter = Some(dispersion.disp_iter);
        fit.dispersion_converged = Some(dispersion.converged);
        fit.mu = Some(dispersion.mu);
        Ok(fit)
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
        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the mean trend.
    ///
    /// This mirrors the implemented subset of `estimateDispersionsGeneEst`
    /// followed by `estimateDispersionsFit(fitType="mean")`. It fills
    /// `dispGeneEst`, a constant `dispFit` for non-all-zero rows, and the
    /// linear fitted mean matrix, but it does not yet estimate prior variance
    /// or MAP dispersions for mean-trend fits.
    pub fn fit_mean_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
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
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_parametric_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the mean trend.
    pub fn fit_mean_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_mean_dispersion_trend(fit)
    }

    /// Run linear-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_linear_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    /// Run GLM-mu gene-wise dispersion estimation and fit the local trend.
    pub fn fit_local_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        let fit = self.fit_gene_wise_dispersions_glm_mu(counts, design)?;
        self.attach_local_dispersion_trend(fit)
    }

    fn attach_parametric_dispersion_trend(
        &self,
        mut fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_parametric_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            ParametricDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..ParametricDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Parametric(trend_fit));
        Ok(fit)
    }

    fn attach_mean_dispersion_trend(&self, mut fit: DeseqFit) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_mean_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            MeanDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..MeanDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Mean(trend_fit));
        Ok(fit)
    }

    fn attach_local_dispersion_trend(&self, mut fit: DeseqFit) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before trend fitting".to_string(),
                })?;
        let trend_fit = fit_local_dispersion_trend(
            &fit.base_mean,
            disp_gene_est,
            LocalDispersionTrendOptions {
                min_disp: self.gene_wise_dispersion_options.min_disp,
                ..LocalDispersionTrendOptions::default()
            },
        )?;
        fit.disp_fit = Some(trend_fit.disp_fit.clone());
        fit.dispersion_trend = Some(DispersionTrendFit::Local(trend_fit));
        Ok(fit)
    }

    /// Run the implemented linear-mu dispersion trend path selected by `fit_type`.
    ///
    /// `Parametric`, `Local`, and `Mean` are currently implemented.
    /// `GlmGamPoi` returns `UnsupportedFeature` until a parity implementation
    /// is added.
    pub fn fit_dispersion_trend_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        match self.fit_type {
            FitType::Parametric => self.fit_parametric_dispersion_trend_linear_mu(counts, design),
            FitType::Mean => self.fit_mean_dispersion_trend_linear_mu(counts, design),
            FitType::Local => self.fit_local_dispersion_trend_linear_mu(counts, design),
            FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
                feature: "linear-mu glmGamPoi dispersion trend fitting".to_string(),
            }),
        }
    }

    /// Run the implemented GLM-mu dispersion trend path selected by `fit_type`.
    ///
    /// `Parametric`, `Local`, and `Mean` are currently implemented.
    /// `GlmGamPoi` returns `UnsupportedFeature` until a parity implementation
    /// is added.
    pub fn fit_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        match self.fit_type {
            FitType::Parametric => self.fit_parametric_dispersion_trend_glm_mu(counts, design),
            FitType::Mean => self.fit_mean_dispersion_trend_glm_mu(counts, design),
            FitType::Local => self.fit_local_dispersion_trend_glm_mu(counts, design),
            FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
                feature: "GLM-mu glmGamPoi dispersion trend fitting".to_string(),
            }),
        }
    }

    /// Fit the selected GLM-mu dispersion trend on DESeq2's fast-VST subset.
    ///
    /// This is a building block for the high-level fast `vst()` workflow:
    /// size factors and normalization factors are derived from the full
    /// dataset, the deterministic fast-VST row subset is selected from the
    /// full-data `baseMean`, and the dispersion trend is fit on the subset
    /// count matrix. The returned subset keeps the original row indices and
    /// aligned normalized counts/factors for inspection.
    pub fn fit_fast_vst_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<(DeseqFit, FastVstSubset), DeseqError> {
        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "fast VST subset size must be positive".to_string(),
            });
        }
        design.validate_full_rank("fast VST GLM-mu dispersion trend")?;
        let stages = self.normalization_stages_for_design(counts, design)?;
        let subset = build_fast_vst_subset(
            counts,
            &stages.normalized,
            &stages.base_mean,
            nsub,
            stages.normalization_factors.as_ref(),
            stages.observation_weights.as_ref(),
        )?;
        let mut subset_builder = self.clone();
        subset_builder.size_factor_options.supplied_size_factors = Some(stages.size_factors);
        subset_builder.normalization_factors = subset.normalization_factors.clone();
        subset_builder.observation_weights = subset.observation_weights.clone();
        let fit = subset_builder.fit_dispersion_trend_glm_mu(&subset.counts, design)?;
        Ok((fit, subset))
    }

    /// Fit the selected GLM-mu dispersion trend on the default fast-VST subset.
    ///
    /// Uses DESeq2's default `nsub=1000`.
    pub fn fit_default_fast_vst_dispersion_trend_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, FastVstSubset), DeseqError> {
        self.fit_fast_vst_dispersion_trend_glm_mu(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply a fast-VST transform using a GLM-mu trend fit on the fast subset.
    ///
    /// This mirrors the implemented part of DESeq2's fast `vst()` workflow:
    /// the dispersion trend is estimated on the deterministic subset, then the
    /// fitted trend is applied to the full normalized count matrix. The subset
    /// fit and row-aligned subset bundle are returned for diagnostics.
    pub fn fast_vst_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<FastVstGlmMuOutput, DeseqError> {
        let (subset_fit, subset) =
            self.fit_fast_vst_dispersion_trend_glm_mu(counts, design, nsub)?;
        let trend_fit =
            subset_fit
                .dispersion_trend
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "a fitted fast-VST dispersion trend is required".to_string(),
                })?;
        let full_fit = Self::base_fit(
            counts,
            Some(design.clone()),
            BaseFitInput {
                size_factors: subset_fit.size_factors.clone(),
                normalization_factors: self.normalization_factors.clone(),
                observation_weights: None,
                weights_fail: None,
                weights_design_rank: None,
                base_mean: vec![f64::NAN; counts.n_genes()],
                base_var: vec![f64::NAN; counts.n_genes()],
                all_zero: counts.all_zero_flags(),
            },
        );
        let transformed = full_fit.variance_stabilizing_transform_with_trend(counts, trend_fit)?;
        Ok(FastVstGlmMuOutput {
            transformed,
            subset_fit,
            subset,
        })
    }

    /// Apply fast VST using DESeq2's default `nsub=1000`.
    pub fn default_fast_vst_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<FastVstGlmMuOutput, DeseqError> {
        self.fast_vst_glm_mu(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply the implemented GLM-mu VST path with DESeq2-like fast-path selection.
    ///
    /// When at least `nsub` rows have `baseMean > 5`, the dispersion trend is
    /// fit on the deterministic fast-VST subset and applied to the full count
    /// matrix. Otherwise, the selected GLM-mu dispersion trend is fit on the
    /// full count matrix before transforming the full matrix.
    pub fn vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        nsub: usize,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        if nsub == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "automatic VST subset size must be positive".to_string(),
            });
        }
        let base_fit = self.fit_size_factors_and_base_means_with_design(counts, design)?;
        let eligible_rows = base_fit.fast_vst_eligible_count()?;
        if eligible_rows >= nsub {
            let output = self.fast_vst_glm_mu(counts, design, nsub)?;
            return Ok(VstGlmMuOutput {
                transformed: output.transformed,
                trend_fit: output.subset_fit,
                trend_source: VstTrendSource::FastSubset {
                    nsub,
                    eligible_rows,
                },
                fast_subset: Some(output.subset),
            });
        }

        let trend_fit = self.fit_dispersion_trend_glm_mu(counts, design)?;
        let transformed = trend_fit.variance_stabilizing_transform(counts)?;
        Ok(VstGlmMuOutput {
            transformed,
            trend_fit,
            trend_source: VstTrendSource::FullData {
                nsub,
                eligible_rows,
                reason: VstFullDataReason::InsufficientEligibleRows,
            },
            fast_subset: None,
        })
    }

    /// Apply automatic GLM-mu VST using DESeq2's default `nsub=1000`.
    pub fn default_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        self.vst_glm_mu_auto(counts, design, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply automatic GLM-mu VST with an intercept-only design.
    ///
    /// This mirrors the implemented part of DESeq2's `blind=TRUE` VST shape:
    /// the transform ignores sample groups by fitting the selected dispersion
    /// trend with a one-column all-ones design, then uses the same automatic
    /// fast-subset/full-data decision as [`Self::vst_glm_mu_auto`].
    pub fn blind_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
        nsub: usize,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        self.vst_glm_mu_auto(counts, &design, nsub)
    }

    /// Apply blind automatic GLM-mu VST using DESeq2's default `nsub=1000`.
    pub fn default_blind_vst_glm_mu_auto(
        &self,
        counts: &CountMatrix,
    ) -> Result<VstGlmMuOutput, DeseqError> {
        self.blind_vst_glm_mu_auto(counts, DEFAULT_FAST_VST_NSUB)
    }

    /// Apply the implemented GLM-mu rlog workflow in one builder call.
    ///
    /// This fits GLM-mu gene-wise, trend, prior, and MAP dispersion stages
    /// using the builder's current options, then applies the fit-state rlog
    /// transform with default IRLS options.
    pub fn rlog_glm_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogOutput, DeseqError> {
        self.rlog_glm_mu_with_fit(counts, design)
            .map(|output| output.rlog)
    }

    /// Apply the implemented GLM-mu rlog workflow with explicit rlog IRLS options.
    pub fn rlog_glm_mu_with_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply the implemented GLM-mu rlog workflow and retain the dispersion fit state.
    pub fn rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        self.rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Apply the implemented GLM-mu rlog workflow with explicit options and retained fit state.
    pub fn rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        let rlog = fit.regularized_log_transform_with_options(counts, rlog_irls_options)?;
        Ok(RlogGlmMuOutput {
            rlog,
            fit,
            design_mode: RlogDesignMode::DesignAware,
        })
    }

    /// Learn rlog intercepts and immediately run a frozen-intercept rlog reuse pass.
    pub fn frozen_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        self.frozen_rlog_glm_mu_with_fit_and_options(counts, design, IrlsOptions::default())
    }

    /// Learn rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        let source = self.rlog_glm_mu_with_fit_and_options(counts, design, rlog_irls_options)?;
        let frozen_rlog = source.fit.frozen_rlog(
            counts,
            &source.rlog.intercept,
            source.rlog.sample_prior_variance,
        )?;
        Ok(FrozenRlogGlmMuOutput {
            source_rlog: source.rlog,
            frozen_rlog,
            fit: source.fit,
            design_mode: source.design_mode,
        })
    }

    /// Apply the implemented GLM-mu rlog workflow with an intercept-only design.
    ///
    /// This mirrors the implemented part of DESeq2's `blind=TRUE` rlog shape:
    /// dispersion fitting and rlog prior estimation ignore sample groups by
    /// using a one-column all-ones design.
    pub fn blind_rlog_glm_mu(&self, counts: &CountMatrix) -> Result<RlogOutput, DeseqError> {
        self.blind_rlog_glm_mu_with_fit(counts)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog with explicit rlog IRLS options.
    pub fn blind_rlog_glm_mu_with_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogOutput, DeseqError> {
        self.blind_rlog_glm_mu_with_fit_and_options(counts, rlog_irls_options)
            .map(|output| output.rlog)
    }

    /// Apply blind GLM-mu rlog and retain the intercept-only dispersion fit state.
    pub fn blind_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        self.blind_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Apply blind GLM-mu rlog with explicit options and retained fit state.
    pub fn blind_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<RlogGlmMuOutput, DeseqError> {
        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        let mut output =
            self.rlog_glm_mu_with_fit_and_options(counts, &design, rlog_irls_options)?;
        output.design_mode = RlogDesignMode::Blind;
        Ok(output)
    }

    /// Learn blind rlog intercepts and run a frozen-intercept rlog reuse pass.
    pub fn blind_frozen_rlog_glm_mu_with_fit(
        &self,
        counts: &CountMatrix,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        self.blind_frozen_rlog_glm_mu_with_fit_and_options(counts, IrlsOptions::default())
    }

    /// Learn blind rlog intercepts and run frozen rlog with explicit IRLS options.
    pub fn blind_frozen_rlog_glm_mu_with_fit_and_options(
        &self,
        counts: &CountMatrix,
        rlog_irls_options: IrlsOptions,
    ) -> Result<FrozenRlogGlmMuOutput, DeseqError> {
        let design = DesignMatrix::intercept_only(counts.n_samples())?;
        let mut output =
            self.frozen_rlog_glm_mu_with_fit_and_options(counts, &design, rlog_irls_options)?;
        output.design_mode = RlogDesignMode::Blind;
        Ok(output)
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
        let fit = self.fit_dispersion_trend_linear_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
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
        let fit = self.fit_parametric_dispersion_trend_linear_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
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
        let fit = self.fit_dispersion_trend_glm_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
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
        let fit = self.fit_parametric_dispersion_trend_glm_mu(counts, design)?;
        self.attach_map_dispersions(counts, design, fit)
    }

    fn attach_map_dispersions(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        mut fit: DeseqFit,
    ) -> Result<DeseqFit, DeseqError> {
        let disp_gene_est =
            fit.disp_gene_est
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidDispersion {
                    reason: "gene-wise dispersions are required before MAP fitting".to_string(),
                })?;
        let disp_fit = fit
            .disp_fit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted dispersion trend is required before MAP fitting".to_string(),
            })?;
        let mu = fit
            .mu
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "fitted means are required before MAP fitting".to_string(),
            })?;
        let prior_variance = estimate_dispersion_prior_variance(
            disp_gene_est,
            disp_fit,
            self.gene_wise_dispersion_options.min_disp,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        let map = estimate_map_dispersions(
            MapDispersionInput {
                counts,
                design,
                mu,
                disp_gene_est,
                disp_fit,
                all_zero: &fit.all_zero,
                observation_weights: fit.observation_weights.as_ref(),
                disp_prior_var: prior_variance.disp_prior_var,
                var_log_disp_estimates: prior_variance.var_log_disp_estimates,
            },
            MapDispersionOptions::from(self.gene_wise_dispersion_options),
        )?;
        fit.disp_prior_var = Some(prior_variance.disp_prior_var);
        fit.disp_map = Some(map.disp_map);
        fit.disp_iter = Some(map.disp_iter);
        fit.disp_outlier = Some(map.disp_outlier);
        fit.dispersion_converged = Some(map.converged);
        fit.dispersion = Some(map.dispersion);
        Ok(fit)
    }

    /// Run the implemented native dispersion path and then a Wald test.
    ///
    /// This is an early, explicitly scoped analogue of `DESeq(..., test="Wald")`:
    /// size factors, base means, linear-mu gene-wise dispersions, selected
    /// trend, deterministic prior variance, MAP dispersions, fixed-dispersion
    /// GLM fitting, Cook's distances, and selected-coefficient Wald results.
    /// Parametric, local, and mean trends are currently implemented. It does not yet
    /// implement DESeq2's general mean/dispersion iteration, beta priors,
    /// contrasts, exact locfit smoothing, or observation weights.
    pub fn fit_wald_linear_mu(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the linear-mu native Wald path for a named primitive contrast specification.
    pub fn fit_wald_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_linear_mu_contrast(counts, design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the linear-mu native Wald path for a caller-supplied factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_linear_mu(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the parametric linear-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_linear_mu_contrast_parametric(counts, design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the parametric linear-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
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
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_glm_mu_contrast(counts, design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_glm_mu(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the linear-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_linear_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_linear_mu_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the linear-mu native LRT path for a factor-level full-model contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_linear_mu(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the parametric linear-mu native LRT path and report a named contrast.
    pub fn fit_lrt_linear_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_linear_mu_contrast_parametric(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the parametric linear-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_linear_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_linear_mu_parametric(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_glm_mu(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_glm_mu(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the GLM-mu native LRT path and report a named full-model contrast.
    pub fn fit_lrt_glm_mu_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) =
            self.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_glm_mu(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        self.attach_native_lrt(counts, full_design, reduced_design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        self.attach_native_lrt_contrast(counts, full_design, reduced_design, contrast, None, fit)
    }

    /// Run the parametric GLM-mu native LRT path and report a named contrast.
    pub fn fit_lrt_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_lrt_glm_mu_contrast_parametric(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
        )?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the parametric GLM-mu native LRT path for a factor-level contrast.
    pub fn fit_lrt_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, full_design)?;
        let (fit, mut results) = self.attach_native_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
    }

    /// Run the current GLM-mu native Wald path with limited Cook's replacement refit.
    ///
    /// This is an explicitly scoped analogue of the `replaceOutliers` /
    /// `refitWithoutOutliers` part of DESeq2 for the currently implemented
    /// native GLM-mu Wald branch. It preserves the original size factors,
    /// builds replacement counts from original Cook's distances, reruns the
    /// implemented GLM-mu dispersion/MAP/Wald path on replacement counts, and
    /// merges only rows marked for refit. It does not yet implement beta
    /// priors or Bioconductor object slots.
    pub fn fit_wald_glm_mu_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_wald_glm_mu(counts, design, coefficient)?;
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        let normalized = match original_fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &original_fit.size_factors)?,
        };
        let refit_plan = prepare_cooks_replacement_refit(
            counts,
            &normalized,
            &original_fit.size_factors,
            original_fit.normalization_factors.as_ref(),
            original_cooks,
            design,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_wald_glm_mu(
                &refit_plan.replacement.replaced_counts,
                design,
                coefficient,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_wald_glm_mu_contrast(counts, design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_wald_glm_mu_contrast(
                &refit_plan.replacement.replaced_counts,
                design,
                contrast,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
    }

    /// Run native GLM-mu Wald replacement refit for a named primitive contrast specification.
    pub fn fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let mut output = self.fit_wald_glm_mu_contrast_with_cooks_replacement(
            counts,
            design,
            &numeric_contrast,
            replacement_options,
        )?;
        apply_contrast_metadata_to_replacement_output(
            &mut output,
            contrast.result_name(),
            contrast.comparison(),
        );
        Ok(output)
    }

    /// Run native GLM-mu Wald replacement refit for a caller-supplied factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_wald_glm_mu_factor_level_contrast(counts, design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_wald_glm_mu_factor_level_contrast(
                &refit_plan.replacement.replaced_counts,
                design,
                contrast,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementWaldOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)?;
        let original_cooks =
            original_fit
                .cooks
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "Cook's distances are required before replacement refit".to_string(),
                })?;
        let normalized = match original_fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &original_fit.size_factors)?,
        };
        let refit_plan = prepare_cooks_replacement_refit(
            counts,
            &normalized,
            &original_fit.size_factors,
            original_fit.normalization_factors.as_ref(),
            original_cooks,
            full_design,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_lrt_glm_mu(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
                coefficient,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) =
            raw_builder.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, contrast)?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            full_design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_lrt_glm_mu_contrast(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
                contrast,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
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
        let raw_builder = self
            .clone()
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        let (original_fit, original_results) = raw_builder.fit_lrt_glm_mu_factor_level_contrast(
            counts,
            full_design,
            reduced_design,
            contrast,
        )?;
        let refit_plan = replacement_refit_plan_from_original(
            counts,
            full_design,
            &original_fit,
            replacement_options,
        )?;

        let (refit_fit, refit_results) = if refit_plan.should_refit {
            let mut refit_builder = raw_builder.clone();
            refit_builder.size_factor_options.supplied_size_factors =
                Some(original_fit.size_factors.clone());
            let (fit, results) = refit_builder.fit_lrt_glm_mu_factor_level_contrast(
                &refit_plan.replacement.replaced_counts,
                full_design,
                reduced_design,
                contrast,
            )?;
            (Some(fit), Some(results))
        } else {
            (None, None)
        };

        let mut results = merge_replacement_refit_results(
            &original_results,
            refit_results.as_ref(),
            &refit_plan,
        )?;
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(CooksReplacementLrtOutput {
            original_fit,
            original_results,
            refit_plan,
            refit_fit,
            refit_results,
            results,
        })
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
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let mut output = self.fit_lrt_glm_mu_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            &numeric_contrast,
            replacement_options,
        )?;
        apply_lrt_contrast_metadata_to_replacement_output(
            &mut output,
            contrast.result_name(),
            contrast.comparison(),
        );
        Ok(output)
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
        validate_pipeline_wald_coefficient(design, coefficient)?;
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        self.attach_native_wald(counts, design, coefficient, fit)
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        self.attach_native_wald_contrast(counts, design, contrast, None, fit)
    }

    /// Run the parametric GLM-mu native Wald path for a named primitive contrast.
    pub fn fit_wald_glm_mu_contrast_spec_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) =
            self.fit_wald_glm_mu_contrast_parametric(counts, design, &numeric_contrast)?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run the parametric GLM-mu native Wald path for a factor-level contrast.
    pub fn fit_wald_glm_mu_factor_level_contrast_parametric(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let fit = self.fit_map_dispersions_glm_mu_parametric(counts, design)?;
        let (fit, mut results) = self.attach_native_wald_contrast(
            counts,
            design,
            &numeric_contrast,
            Some(&contrast_all_zero),
            fit,
        )?;
        results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, results))
    }

    fn attach_native_wald(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient: usize,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before Wald fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let wald_output = self.fixed_dispersion_wald_components(WaldPipelineInput {
            counts,
            design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient,
        })?;

        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald);
        Ok((fit, wald_output.results))
    }

    fn attach_native_wald_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before Wald fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let wald_output = self.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &fit.size_factors,
                normalization_factors: fit.normalization_factors.as_ref(),
                observation_weights: fit.observation_weights.as_ref(),
                all_zero: &fit.all_zero,
                dispersions,
            },
            &normalized,
            &fit.base_mean,
            contrast,
            contrast_all_zero_override,
        )?;

        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        Ok((fit, wald_output.results))
    }

    fn attach_native_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient: usize,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before LRT fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient,
        })?;

        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    fn attach_native_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
        mut fit: DeseqFit,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let dispersions = fit
            .dispersion
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "MAP dispersions are required before LRT fitting".to_string(),
            })?;
        let normalized = match fit.normalization_factors.as_ref() {
            Some(normalization_factors) => {
                normalized_counts_with_factors(counts, normalization_factors)?
            }
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            normalized: &normalized,
            base_mean: &fit.base_mean,
            all_zero: &fit.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &fit.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        let contrast_all_zero = match contrast_all_zero_override {
            Some(flags) => {
                if flags.len() != counts.n_genes() {
                    return Err(invalid_dimensions(
                        "contrastAllZero rows",
                        counts.n_genes(),
                        flags.len(),
                    ));
                }
                flags.to_vec()
            }
            None => contrast_all_zero_numeric(counts, full_design, contrast)?,
        };
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &fit.all_zero,
        )?;
        for (gene, all_zero) in fit.all_zero.iter().copied().enumerate() {
            lrt_output.results.rows[gene].max_cooks = lrt_output.cooks.max_cooks[gene];
            if all_zero {
                lrt_output.results.rows[gene].converged = None;
                lrt_output.results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut lrt_output.results, cooks_cutoff)?;
        apply_independent_filtering(&mut lrt_output.results, &self.independent_filtering_options)?;

        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    /// Run a supplied-dispersion Wald pipeline for one coefficient.
    ///
    /// This is the current closest analogue to the core `nbinomWaldTest` path,
    /// but it requires caller-supplied dispersions and does not yet implement
    /// contrasts or beta priors.
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
        let stages = self.normalization_stages_for_design(counts, design)?;
        let wald_output = self.fixed_dispersion_wald_components(WaldPipelineInput {
            counts,
            design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient,
        })?;
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald);
        Ok((fit, wald_output.results))
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
        let stages = self.normalization_stages_for_design(counts, design)?;
        let wald_output = self.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                observation_weights: stages.observation_weights.as_ref(),
                all_zero: &stages.all_zero,
                dispersions,
            },
            &stages.normalized,
            &stages.base_mean,
            contrast,
            None,
        )?;
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        Ok((fit, wald_output.results))
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
        let numeric_contrast = resolve_contrast(design, contrast)?;
        let (fit, mut results) = self.fit_fixed_dispersion_wald_contrast(
            counts,
            design,
            dispersions,
            &numeric_contrast,
        )?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
    }

    /// Run a supplied-dispersion Wald pipeline for a factor-level contrast.
    ///
    /// This resolves DESeq2-shaped coefficient names from the design matrix and
    /// applies DESeq2-style character contrast all-zero handling using
    /// caller-supplied sample levels. R formula parsing and colData ownership
    /// remain outside the Rust core.
    pub fn fit_fixed_dispersion_wald_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        dispersions: &[f64],
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
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
        let stages = self.normalization_stages_for_design(counts, design)?;
        let mut wald_output = self.fixed_dispersion_wald_contrast_components(
            FixedDispersionGlmInput {
                counts,
                design,
                size_factors: &stages.size_factors,
                normalization_factors: stages.normalization_factors.as_ref(),
                observation_weights: stages.observation_weights.as_ref(),
                all_zero: &stages.all_zero,
                dispersions,
            },
            &stages.normalized,
            &stages.base_mean,
            &numeric_contrast,
            Some(&contrast_all_zero),
        )?;
        let mut fit = Self::base_fit(counts, Some(design.clone()), stages.into_base_fit_input());
        fit.dispersion = Some(wald_output.expanded_dispersions);
        fit.cooks = Some(wald_output.cooks.cooks);
        fit.max_cooks = Some(wald_output.cooks.max_cooks);
        attach_glm_fit(&mut fit, wald_output.glm_fit);
        fit.wald = Some(wald_output.wald_contrast.wald);
        wald_output.results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        wald_output.results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        Ok((fit, wald_output.results))
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
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient,
        })?;

        let mut fit = Self::base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
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
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &stages.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, full_design, contrast)?;
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &stages.all_zero,
        )?;
        for (gene, all_zero) in stages.all_zero.iter().copied().enumerate() {
            lrt_output.results.rows[gene].max_cooks = lrt_output.cooks.max_cooks[gene];
            if all_zero {
                lrt_output.results.rows[gene].converged = None;
                lrt_output.results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut lrt_output.results, cooks_cutoff)?;
        apply_independent_filtering(&mut lrt_output.results, &self.independent_filtering_options)?;

        let mut fit = Self::base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
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
        let numeric_contrast = resolve_contrast(full_design, contrast)?;
        let (fit, mut results) = self.fit_fixed_dispersion_lrt_contrast(
            counts,
            full_design,
            reduced_design,
            dispersions,
            &numeric_contrast,
        )?;
        results.metadata.result_name = Some(contrast.result_name());
        results.metadata.comparison = Some(contrast.comparison());
        Ok((fit, results))
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
        let stages = self.normalization_stages_for_design(counts, full_design)?;
        let mut lrt_output = self.fixed_dispersion_lrt_components(LrtPipelineInput {
            counts,
            full_design,
            reduced_design,
            size_factors: &stages.size_factors,
            normalization_factors: stages.normalization_factors.as_ref(),
            observation_weights: stages.observation_weights.as_ref(),
            normalized: &stages.normalized,
            base_mean: &stages.base_mean,
            all_zero: &stages.all_zero,
            dispersions,
            coefficient: default_results_coefficient(full_design)?,
        })?;
        let contrast_output = wald_test_contrast_with_options(
            &lrt_output.full_fit,
            &numeric_contrast,
            &self.wald_test_options,
        )?;
        lrt_output.results = build_lrt_contrast_results(
            &stages.base_mean,
            &lrt_output.full_fit,
            &lrt_output.lrt,
            &contrast_output,
            counts.gene_names(),
            Some(&lrt_output.expanded_dispersions),
        )?;
        apply_contrast_all_zero_to_lrt_results(
            &mut lrt_output.results,
            &contrast_all_zero,
            &stages.all_zero,
        )?;
        for (gene, all_zero) in stages.all_zero.iter().copied().enumerate() {
            lrt_output.results.rows[gene].max_cooks = lrt_output.cooks.max_cooks[gene];
            if all_zero {
                lrt_output.results.rows[gene].converged = None;
                lrt_output.results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            full_design.n_samples(),
            full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut lrt_output.results, cooks_cutoff)?;
        apply_independent_filtering(&mut lrt_output.results, &self.independent_filtering_options)?;
        lrt_output.results.metadata.result_name = Some(format!(
            "{}_{}_vs_{}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));
        lrt_output.results.metadata.comparison = Some(format!(
            "factor-level contrast: {} {} vs {}",
            contrast.factor, contrast.numerator, contrast.denominator
        ));

        let mut fit = Self::base_fit(
            counts,
            Some(full_design.clone()),
            stages.into_base_fit_input(),
        );
        fit.reduced_design = Some(reduced_design.clone());
        fit.dispersion = Some(lrt_output.expanded_dispersions);
        fit.cooks = Some(lrt_output.cooks.cooks);
        fit.max_cooks = Some(lrt_output.cooks.max_cooks);
        fit.reduced_log_like = Some(lrt_output.reduced_fit.log_like.clone());
        fit.reduced_beta_converged = Some(lrt_output.reduced_fit.beta_converged.clone());
        fit.reduced_beta_iter = Some(lrt_output.reduced_fit.beta_iter.clone());
        fit.reduced_mu = Some(lrt_output.reduced_fit.mu.clone());
        fit.reduced_hat_diagonal = Some(lrt_output.reduced_fit.hat_diagonal.clone());
        attach_glm_fit(&mut fit, lrt_output.full_fit);
        fit.lrt = Some(lrt_output.lrt);
        Ok((fit, lrt_output.results))
    }

    fn normalization_stages(
        &self,
        counts: &CountMatrix,
    ) -> Result<NormalizationStages, DeseqError> {
        self.normalization_stages_inner(counts, None)
    }

    fn normalization_stages_for_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<NormalizationStages, DeseqError> {
        self.normalization_stages_inner(counts, Some(design))
    }

    fn normalization_stages_inner(
        &self,
        counts: &CountMatrix,
        design: Option<&DesignMatrix>,
    ) -> Result<NormalizationStages, DeseqError> {
        let (size_factors, normalization_factors, normalized) = match &self.normalization_factors {
            Some(factors) => {
                validate_normalization_factors(counts, factors)?;
                let size_factors = match &self.size_factor_options.supplied_size_factors {
                    Some(size_factors) => {
                        normalized_counts(counts, size_factors)?;
                        size_factors.clone()
                    }
                    None => vec![1.0; counts.n_samples()],
                };
                let factors = factors.clone();
                let normalized = normalized_counts_with_factors(counts, &factors)?;
                (size_factors, Some(factors), normalized)
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
                let normalized = normalized_counts(counts, &size_factors)?;
                (size_factors, None, normalized)
            }
        };
        let weighted_metadata = self.weighted_base_metadata(counts, design, &normalized)?;
        let raw_all_zero = counts.all_zero_flags();
        let all_zero = match &weighted_metadata.weights_fail {
            Some(weights_fail) => raw_all_zero
                .iter()
                .copied()
                .zip(weights_fail.iter().copied())
                .map(|(all_zero, weights_fail)| all_zero || weights_fail)
                .collect(),
            None => raw_all_zero,
        };
        Ok(NormalizationStages {
            size_factors,
            base_mean: weighted_metadata.base_mean,
            base_var: weighted_metadata.base_var,
            all_zero,
            normalized,
            normalization_factors,
            observation_weights: weighted_metadata.observation_weights,
            weights_fail: weighted_metadata.weights_fail,
            weights_design_rank: weighted_metadata.weights_design_rank,
        })
    }

    fn weighted_base_metadata(
        &self,
        counts: &CountMatrix,
        design: Option<&DesignMatrix>,
        normalized: &RowMajorMatrix<f64>,
    ) -> Result<WeightedBaseMetadata, DeseqError> {
        let Some(weights) = &self.observation_weights else {
            return Ok(WeightedBaseMetadata {
                base_mean: base_mean(normalized)?,
                base_var: base_variance(normalized)?,
                observation_weights: None,
                weights_fail: None,
                weights_design_rank: None,
            });
        };
        validate_observation_weights_for_counts(counts, weights)?;
        match design {
            Some(design) => {
                let checked = preprocess_observation_weights_with_options(
                    weights,
                    design,
                    self.observation_weight_options,
                )?;
                let base_mean = base_mean_with_weights(normalized, weights)?;
                let base_var = base_variance_with_weights(normalized, weights)?;
                Ok(WeightedBaseMetadata {
                    base_mean,
                    base_var,
                    observation_weights: Some(checked.weights),
                    weights_fail: Some(checked.weights_fail),
                    weights_design_rank: Some(checked.design_rank),
                })
            }
            None => Ok(WeightedBaseMetadata {
                base_mean: base_mean_with_weights(normalized, weights)?,
                base_var: base_variance_with_weights(normalized, weights)?,
                observation_weights: Some(weights.clone()),
                weights_fail: None,
                weights_design_rank: None,
            }),
        }
    }

    fn ensure_no_observation_weights(&self, feature: &str) -> Result<(), DeseqError> {
        if self.observation_weights.is_some() {
            return Err(DeseqError::UnsupportedFeature {
                feature: format!("{feature} with observation weights"),
            });
        }
        Ok(())
    }

    fn fixed_dispersion_wald_components(
        &self,
        input: WaldPipelineInput<'_>,
    ) -> Result<WaldPipelineOutput, DeseqError> {
        if input.coefficient >= input.design.n_coefficients() {
            return Err(DeseqError::InvalidDimensions {
                context: "pipeline Wald coefficient index".to_string(),
                expected: input.design.n_coefficients().saturating_sub(1),
                actual: input.coefficient,
            });
        }
        let FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        } = self.fixed_dispersion_glm_components(FixedDispersionGlmInput {
            counts: input.counts,
            design: input.design,
            size_factors: input.size_factors,
            normalization_factors: input.normalization_factors,
            observation_weights: input.observation_weights,
            all_zero: input.all_zero,
            dispersions: input.dispersions,
        })?;
        let mut wald = wald_test_coefficient_with_options(
            &glm_fit,
            input.coefficient,
            &self.wald_test_options,
        )?;
        mask_wald_degrees_of_freedom_for_all_zero_rows(&mut wald, input.all_zero)?;
        let cooks = calculate_cooks_distance(
            input.counts,
            input.normalized,
            &glm_fit.mu,
            &glm_fit.hat_diagonal,
            input.design,
        )?;
        let mut results = build_wald_results_from_wald(
            input.base_mean,
            &glm_fit,
            input.coefficient,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
            &wald,
        )?;
        results.apply_wald_test_options(&self.wald_test_options);
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.design.n_samples(),
            input.design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(WaldPipelineOutput {
            glm_fit,
            wald,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_wald_contrast_components(
        &self,
        input: FixedDispersionGlmInput<'_>,
        normalized: &RowMajorMatrix<f64>,
        base_mean: &[f64],
        contrast: &[f64],
        contrast_all_zero_override: Option<&[bool]>,
    ) -> Result<WaldContrastPipelineOutput, DeseqError> {
        let FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        } = self.fixed_dispersion_glm_components(input)?;
        let mut wald_contrast =
            wald_test_contrast_with_options(&glm_fit, contrast, &self.wald_test_options)?;
        mask_wald_degrees_of_freedom_for_all_zero_rows(&mut wald_contrast.wald, input.all_zero)?;
        let contrast_all_zero = match contrast_all_zero_override {
            Some(flags) => {
                if flags.len() != input.counts.n_genes() {
                    return Err(invalid_dimensions(
                        "contrastAllZero rows",
                        input.counts.n_genes(),
                        flags.len(),
                    ));
                }
                flags.to_vec()
            }
            None => contrast_all_zero_numeric(input.counts, input.design, contrast)?,
        };
        apply_contrast_all_zero_to_wald_contrast(
            &mut wald_contrast,
            &contrast_all_zero,
            input.all_zero,
        )?;
        let cooks = calculate_cooks_distance(
            input.counts,
            normalized,
            &glm_fit.mu,
            &glm_fit.hat_diagonal,
            input.design,
        )?;
        let mut results = build_wald_contrast_results(
            base_mean,
            &glm_fit,
            &wald_contrast,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
        )?;
        results.apply_wald_test_options(&self.wald_test_options);
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.design.n_samples(),
            input.design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(WaldContrastPipelineOutput {
            glm_fit,
            wald_contrast,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_lrt_components(
        &self,
        input: LrtPipelineInput<'_>,
    ) -> Result<LrtPipelineOutput, DeseqError> {
        validate_lrt_pipeline_input(&input)?;
        let nonzero_gene_indices = nonzero_gene_indices(input.all_zero);
        let (full_fit, reduced_fit) = if nonzero_gene_indices.is_empty() {
            (
                all_zero_glm_fit(input.counts, input.full_design)?,
                all_zero_glm_fit(input.counts, input.reduced_design)?,
            )
        } else {
            let compact_counts = compact_counts(input.counts, &nonzero_gene_indices)?;
            let compact_dispersions = nonzero_gene_indices
                .iter()
                .map(|gene| input.dispersions[*gene])
                .collect::<Vec<_>>();
            let compact_normalization_factors = input
                .normalization_factors
                .map(|factors| compact_matrix_rows(factors, &nonzero_gene_indices))
                .transpose()?;
            let compact_weights = input
                .observation_weights
                .map(|weights| compact_matrix_rows(weights, &nonzero_gene_indices))
                .transpose()?;
            let compact_full_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.full_design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            let compact_reduced_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.reduced_design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            (
                expand_glm_fit(compact_full_fit, input.all_zero)?,
                expand_glm_fit(compact_reduced_fit, input.all_zero)?,
            )
        };

        let expanded_dispersions =
            mask_all_zero_values_with_nan_rows(input.dispersions, input.all_zero)?;
        let lrt = lrt_test(&full_fit, &reduced_fit)?;
        let cooks = calculate_cooks_distance(
            input.counts,
            input.normalized,
            &full_fit.mu,
            &full_fit.hat_diagonal,
            input.full_design,
        )?;
        let mut results = build_lrt_results(
            input.base_mean,
            &full_fit,
            &lrt,
            input.coefficient,
            input.counts.gene_names(),
            Some(&expanded_dispersions),
        )?;
        for (gene, all_zero) in input.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            input.full_design.n_samples(),
            input.full_design.n_coefficients(),
        )?;
        apply_cooks_cutoff(&mut results, cooks_cutoff)?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;

        Ok(LrtPipelineOutput {
            full_fit,
            reduced_fit,
            lrt,
            cooks,
            results,
            expanded_dispersions,
        })
    }

    fn fixed_dispersion_glm_components(
        &self,
        input: FixedDispersionGlmInput<'_>,
    ) -> Result<FixedDispersionGlmOutput, DeseqError> {
        input.design.validate_full_rank("GLM")?;
        if input.dispersions.len() != input.counts.n_genes() {
            return Err(invalid_dimensions(
                "pipeline dispersions",
                input.counts.n_genes(),
                input.dispersions.len(),
            ));
        }

        let nonzero_gene_indices = nonzero_gene_indices(input.all_zero);
        let glm_fit = if nonzero_gene_indices.is_empty() {
            all_zero_glm_fit(input.counts, input.design)?
        } else {
            let compact_counts = compact_counts(input.counts, &nonzero_gene_indices)?;
            let compact_dispersions = nonzero_gene_indices
                .iter()
                .map(|gene| input.dispersions[*gene])
                .collect::<Vec<_>>();
            let compact_normalization_factors = input
                .normalization_factors
                .map(|factors| compact_matrix_rows(factors, &nonzero_gene_indices))
                .transpose()?;
            let compact_weights = input
                .observation_weights
                .map(|weights| compact_matrix_rows(weights, &nonzero_gene_indices))
                .transpose()?;
            let compact_fit = fit_fixed_dispersion_model(
                &compact_counts,
                input.design,
                input.size_factors,
                compact_normalization_factors.as_ref(),
                compact_weights.as_ref(),
                &compact_dispersions,
                self.irls_options.clone(),
            )?;
            expand_glm_fit(compact_fit, input.all_zero)?
        };
        let expanded_dispersions =
            mask_all_zero_values_with_nan_rows(input.dispersions, input.all_zero)?;
        Ok(FixedDispersionGlmOutput {
            glm_fit,
            expanded_dispersions,
        })
    }

    fn base_fit(
        counts: &CountMatrix,
        design: Option<DesignMatrix>,
        input: BaseFitInput,
    ) -> DeseqFit {
        DeseqFit {
            counts_summary: counts.summary(),
            design,
            reduced_design: None,
            size_factors: input.size_factors,
            normalization_factors: input.normalization_factors,
            observation_weights: input.observation_weights,
            weights_fail: input.weights_fail,
            weights_design_rank: input.weights_design_rank,
            base_mean: input.base_mean,
            base_var: input.base_var,
            all_zero: input.all_zero,
            disp_gene_est: None,
            disp_gene_iter: None,
            disp_fit: None,
            dispersion_trend: None,
            disp_map: None,
            dispersion: None,
            disp_iter: None,
            disp_outlier: None,
            disp_prior_var: None,
            dispersion_converged: None,
            beta: None,
            beta_se: None,
            beta_covariance: None,
            beta_converged: None,
            beta_iter: None,
            log_like: None,
            full_deviance: None,
            reduced_log_like: None,
            reduced_beta_converged: None,
            reduced_beta_iter: None,
            reduced_mu: None,
            reduced_hat_diagonal: None,
            mu: None,
            cooks: None,
            max_cooks: None,
            hat_diagonal: None,
            wald: None,
            lrt: None,
        }
    }

    /// Run the currently implemented DESeq-like workflow.
    ///
    /// For `test=Wald`, this follows the implemented GLM-mu native path and
    /// reports the last design coefficient, matching DESeq2's default
    /// coefficient selection shape. For `test=Lrt`, callers must first store a
    /// reduced design with [`DeseqBuilder::reduced_design`].
    pub fn fit(&self, counts: &CountMatrix, design: &DesignMatrix) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results(counts, design)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented DESeq-like workflow and return result rows.
    ///
    /// This is the result-table-producing companion to [`DeseqBuilder::fit`].
    pub fn fit_with_results(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results(counts, design, reduced_design)
            }
        }
    }

    /// Run the top-level Wald workflow and report a named design coefficient.
    pub fn fit_with_results_name(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = design.coefficient_index(coefficient_name)?;
                self.fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_name(counts, design, reduced_design, coefficient_name)
            }
        }
    }

    /// Run the top-level Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.fit_wald_glm_mu_with_cooks_replacement(
                    counts,
                    design,
                    coefficient,
                    replacement_options,
                )
            }
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            }),
        }
    }

    /// Run the top-level named Wald workflow with limited Cook's replacement refit.
    pub fn fit_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => {
                let coefficient = design.coefficient_index(coefficient_name)?;
                self.fit_wald_glm_mu_with_cooks_replacement(
                    counts,
                    design,
                    coefficient,
                    replacement_options,
                )
            }
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit without a reduced design".to_string(),
            }),
        }
    }

    /// Run the top-level workflow with limited Cook's replacement refit.
    ///
    /// Unlike [`Self::fit_with_results_with_cooks_replacement`], this method
    /// returns an enum so `test=Lrt` can route through a stored reduced design.
    pub fn fit_with_test_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_with_cooks_replacement(counts, design, replacement_options)
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the top-level named workflow with limited Cook's replacement refit.
    ///
    /// The returned enum keeps Wald and LRT replacement output types explicit
    /// while allowing `test=Lrt` to use the builder's stored reduced design.
    pub fn fit_with_test_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_name_with_cooks_replacement(
                    counts,
                    design,
                    coefficient_name,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_name_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    coefficient_name,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the currently implemented top-level Wald workflow for a numeric contrast.
    ///
    /// This is the primitive contrast companion to [`Self::fit_with_results`].
    /// It follows the implemented GLM-mu native Wald path when `test=Wald`.
    /// Top-level LRT remains explicit because a contrast is not an LRT reduced
    /// model.
    pub fn fit_with_results_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast(counts, design, reduced_design, contrast)
            }
        }
    }

    /// Run the currently implemented top-level Wald workflow for a numeric contrast.
    pub fn fit_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast".to_string(),
            }),
        }
    }

    /// Run the top-level numeric-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_spec(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_spec(counts, design, reduced_design, contrast)
            }
        }
    }

    /// Run the top-level Wald workflow for a named primitive contrast specification.
    pub fn fit_contrast_spec(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_contrast_spec(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level named Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald contrast specification"
                    .to_string(),
            }),
        }
    }

    /// Run the top-level named-contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_contrast_spec_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_factor_level_contrast(counts, design, contrast),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_factor_level_contrast(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                )
            }
        }
    }

    /// Run the top-level Wald workflow for a caller-supplied factor-level contrast.
    pub fn fit_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_with_results_factor_level_contrast(counts, design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the top-level factor-level Wald contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementWaldOutput, DeseqError> {
        match self.test {
            TestType::Wald => self.fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
                counts,
                design,
                contrast,
                replacement_options,
            ),
            TestType::Lrt => Err(DeseqError::UnsupportedFeature {
                feature: "top-level LRT replacement refit for a Wald factor-level contrast"
                    .to_string(),
            }),
        }
    }

    /// Run the top-level factor-level contrast workflow with limited Cook's replacement refit.
    pub fn fit_with_test_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementTestOutput, DeseqError> {
        match self.test {
            TestType::Wald => self
                .fit_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Wald),
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                    counts,
                    design,
                    reduced_design,
                    contrast,
                    replacement_options,
                )
                .map(CooksReplacementTestOutput::Lrt)
            }
        }
    }

    /// Run the currently implemented top-level LRT workflow with a reduced design.
    ///
    /// This follows the implemented GLM-mu native LRT path and reports the last
    /// full-design coefficient by default.
    pub fn fit_lrt(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results(counts, full_design, reduced_design)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and return result rows.
    pub fn fit_lrt_with_results(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let coefficient = default_results_coefficient(full_design)?;
        self.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    pub fn fit_lrt_with_results_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let coefficient = full_design.coefficient_index(coefficient_name)?;
        self.fit_lrt_glm_mu(counts, full_design, reduced_design, coefficient)
    }

    /// Run the currently implemented top-level LRT workflow and report a named full-design coefficient.
    pub fn fit_lrt_name(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_name(counts, full_design, reduced_design, coefficient_name)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    pub fn fit_lrt_with_results_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_contrast(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow and report a numeric contrast.
    pub fn fit_lrt_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_contrast(counts, full_design, reduced_design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    pub fn fit_lrt_with_results_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_contrast_spec(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    pub fn fit_lrt_with_results_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        self.fit_lrt_glm_mu_factor_level_contrast(counts, full_design, reduced_design, contrast)
    }

    /// Run the currently implemented top-level LRT workflow and report a named contrast.
    pub fn fit_lrt_contrast_spec(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_contrast_spec(counts, full_design, reduced_design, contrast)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow for a caller-supplied factor-level contrast.
    pub fn fit_lrt_factor_level_contrast(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
    ) -> Result<DeseqFit, DeseqError> {
        self.fit_lrt_with_results_factor_level_contrast(
            counts,
            full_design,
            reduced_design,
            contrast,
        )
        .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented top-level LRT workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let coefficient = default_results_coefficient(full_design)?;
        self.fit_lrt_glm_mu_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            coefficient,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level LRT replacement-refit workflow and report a named full-design coefficient.
    pub fn fit_lrt_with_results_name_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        coefficient_name: &str,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        let coefficient = full_design.coefficient_index(coefficient_name)?;
        self.fit_lrt_glm_mu_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            coefficient,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &[f64],
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level named LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_contrast_spec_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: &ContrastSpec,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    /// Run the currently implemented top-level factor-level LRT contrast workflow with limited Cook's replacement refit.
    pub fn fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
        &self,
        counts: &CountMatrix,
        full_design: &DesignMatrix,
        reduced_design: &DesignMatrix,
        contrast: FactorLevelContrast<'_>,
        replacement_options: &CooksReplacementOptions,
    ) -> Result<CooksReplacementLrtOutput, DeseqError> {
        self.fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
            counts,
            full_design,
            reduced_design,
            contrast,
            replacement_options,
        )
    }

    fn reduced_design_for_top_level_lrt(&self) -> Result<&DesignMatrix, DeseqError> {
        self.reduced_design
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "top-level LRT fit without a reduced design".to_string(),
            })
    }
}

fn default_results_coefficient(design: &DesignMatrix) -> Result<usize, DeseqError> {
    design
        .n_coefficients()
        .checked_sub(1)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "default results coefficient".to_string(),
            expected: 1,
            actual: 0,
        })
}

fn validate_observation_weights_for_counts(
    counts: &CountMatrix,
    weights: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
        return Err(DeseqError::InvalidDimensions {
            context: "observation weights".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: weights.len(),
        });
    }
    Ok(())
}

fn validate_lrt_pipeline_input(input: &LrtPipelineInput<'_>) -> Result<(), DeseqError> {
    if input.dispersions.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "pipeline dispersions",
            input.counts.n_genes(),
            input.dispersions.len(),
        ));
    }
    if input.full_design.n_samples() != input.reduced_design.n_samples() {
        return Err(invalid_dimensions(
            "LRT reduced design samples",
            input.full_design.n_samples(),
            input.reduced_design.n_samples(),
        ));
    }
    if input.full_design.n_coefficients() <= input.reduced_design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT full/reduced coefficients".to_string(),
            expected: input.reduced_design.n_coefficients() + 1,
            actual: input.full_design.n_coefficients(),
        });
    }
    if input.coefficient >= input.full_design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT result coefficient index".to_string(),
            expected: input.full_design.n_coefficients().saturating_sub(1),
            actual: input.coefficient,
        });
    }
    if input.all_zero.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "LRT all-zero rows",
            input.counts.n_genes(),
            input.all_zero.len(),
        ));
    }
    if input.base_mean.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "LRT baseMean rows",
            input.counts.n_genes(),
            input.base_mean.len(),
        ));
    }
    if input.normalized.n_rows() != input.counts.n_genes()
        || input.normalized.n_cols() != input.counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT normalized counts".to_string(),
            expected: input.counts.n_genes() * input.counts.n_samples(),
            actual: input.normalized.len(),
        });
    }
    input.full_design.validate_full_rank("LRT full")?;
    input.reduced_design.validate_full_rank("LRT reduced")?;
    Ok(())
}

fn is_intercept_only_design(design: &DesignMatrix) -> bool {
    design.n_coefficients() == 1
        && design
            .matrix()
            .as_slice()
            .iter()
            .all(|value| (*value - 1.0).abs() <= f64::EPSILON)
}

fn fit_fixed_dispersion_model(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    normalization_factors: Option<&RowMajorMatrix<f64>>,
    weights: Option<&RowMajorMatrix<f64>>,
    dispersions: &[f64],
    irls_options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    if is_intercept_only_design(design)
        && irls_options.uses_intercept_shortcut_for_coefficients(design.n_coefficients())?
    {
        match normalization_factors {
            Some(factors) => fit_intercept_only_fixed_dispersion_with_normalization_factors(
                counts,
                factors,
                dispersions,
                weights,
            ),
            None => fit_intercept_only_fixed_dispersion_with_weights(
                counts,
                size_factors,
                dispersions,
                weights,
            ),
        }
    } else {
        match normalization_factors {
            Some(factors) => fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
                counts,
                design,
                factors,
                dispersions,
                weights,
                irls_options,
            ),
            None => fit_fixed_dispersion_irls_with_weights(
                counts,
                design,
                size_factors,
                dispersions,
                weights,
                irls_options,
            ),
        }
    }
}

fn merge_replacement_refit_results(
    original_results: &DeseqResults,
    refit_results: Option<&DeseqResults>,
    refit_plan: &CooksRefitPlan,
) -> Result<DeseqResults, DeseqError> {
    if original_results.rows.len() != refit_plan.replacement.replace.len() {
        return Err(invalid_dimensions(
            "replacement-refit result rows",
            refit_plan.replacement.replace.len(),
            original_results.rows.len(),
        ));
    }
    if let Some(refit_results) = refit_results {
        if refit_results.rows.len() != original_results.rows.len() {
            return Err(invalid_dimensions(
                "replacement-refit refit result rows",
                original_results.rows.len(),
                refit_results.rows.len(),
            ));
        }
    }
    if refit_plan.replaced_base_mean.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "replacement-refit baseMean rows",
            original_results.rows.len(),
            refit_plan.replaced_base_mean.len(),
        ));
    }
    if refit_plan.post_refit_max_cooks.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "replacement-refit maxCooks rows",
            original_results.rows.len(),
            refit_plan.post_refit_max_cooks.len(),
        ));
    }

    let mut merged = original_results.clone();
    for (gene, row) in merged.rows.iter_mut().enumerate() {
        row.base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            row.max_cooks = refit_plan.post_refit_max_cooks[gene];
            row.cooks_outlier = None;
            row.filtered = None;
        }
    }

    if let Some(refit_results) = refit_results {
        for gene in refit_plan.refit_rows.iter().copied() {
            merged.rows[gene] = refit_results.rows[gene].clone();
            merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
            merged.rows[gene].cooks_outlier = None;
            merged.rows[gene].filtered = None;
        }
    }

    for gene in refit_plan.new_all_zero_rows.iter().copied() {
        clear_replacement_all_zero_result(&mut merged.rows[gene]);
        merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
        }
    }

    merged.independent_filtering = None;
    Ok(merged)
}

fn replacement_refit_plan_from_original(
    counts: &CountMatrix,
    design: &DesignMatrix,
    original_fit: &DeseqFit,
    replacement_options: &CooksReplacementOptions,
) -> Result<CooksRefitPlan, DeseqError> {
    let original_cooks = original_fit
        .cooks
        .as_ref()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "Cook's distances are required before replacement refit".to_string(),
        })?;
    let normalized = match original_fit.normalization_factors.as_ref() {
        Some(normalization_factors) => {
            normalized_counts_with_factors(counts, normalization_factors)?
        }
        None => normalized_counts(counts, &original_fit.size_factors)?,
    };
    prepare_cooks_replacement_refit(
        counts,
        &normalized,
        &original_fit.size_factors,
        original_fit.normalization_factors.as_ref(),
        original_cooks,
        design,
        replacement_options,
    )
}

fn apply_contrast_metadata_to_replacement_output(
    output: &mut CooksReplacementWaldOutput,
    result_name: String,
    comparison: String,
) {
    output.original_results.metadata.result_name = Some(result_name.clone());
    output.original_results.metadata.comparison = Some(comparison.clone());
    if let Some(refit_results) = &mut output.refit_results {
        refit_results.metadata.result_name = Some(result_name.clone());
        refit_results.metadata.comparison = Some(comparison.clone());
    }
    output.results.metadata.result_name = Some(result_name);
    output.results.metadata.comparison = Some(comparison);
}

fn apply_lrt_contrast_metadata_to_replacement_output(
    output: &mut CooksReplacementLrtOutput,
    result_name: String,
    comparison: String,
) {
    output.original_results.metadata.result_name = Some(result_name.clone());
    output.original_results.metadata.comparison = Some(comparison.clone());
    if let Some(refit_results) = &mut output.refit_results {
        refit_results.metadata.result_name = Some(result_name.clone());
        refit_results.metadata.comparison = Some(comparison.clone());
    }
    output.results.metadata.result_name = Some(result_name);
    output.results.metadata.comparison = Some(comparison);
}

fn clear_replacement_all_zero_result(row: &mut DeseqResultRow) {
    row.log2_fold_change = Some(0.0);
    row.lfc_se = Some(0.0);
    row.stat = Some(0.0);
    row.pvalue = Some(1.0);
    row.padj = None;
    row.dispersion = None;
    row.converged = None;
    row.cooks_outlier = None;
    row.filtered = None;
}

fn compact_counts(counts: &CountMatrix, gene_indices: &[usize]) -> Result<CountMatrix, DeseqError> {
    counts.select_rows(gene_indices)
}

fn compact_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    row_indices: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * matrix.n_cols());
    for row in row_indices {
        values.extend_from_slice(matrix.row(*row)?);
    }
    RowMajorMatrix::from_row_major(row_indices.len(), matrix.n_cols(), values)
}

fn compact_f64_values(values: &[f64], row_indices: &[usize]) -> Result<Vec<f64>, DeseqError> {
    let mut compact = Vec::with_capacity(row_indices.len());
    for row in row_indices {
        let Some(value) = values.get(*row) else {
            return Err(invalid_dimensions(
                "compact vector rows",
                row + 1,
                values.len(),
            ));
        };
        compact.push(*value);
    }
    Ok(compact)
}

fn expand_rlog_output_with_all_zero_rows(
    compact_output: RlogOutput,
    all_zero: &[bool],
    n_samples: usize,
) -> Result<RlogOutput, DeseqError> {
    if compact_output.transformed.n_cols() != n_samples {
        return Err(invalid_dimensions(
            "expanded rlog columns",
            n_samples,
            compact_output.transformed.n_cols(),
        ));
    }
    let mut values = vec![0.0; all_zero.len() * n_samples];
    let mut compact_row = 0_usize;
    for (gene, is_zero) in all_zero.iter().copied().enumerate() {
        if is_zero {
            continue;
        }
        let src = compact_output.transformed.row(compact_row)?;
        let start = gene * n_samples;
        values[start..start + n_samples].copy_from_slice(src);
        compact_row += 1;
    }
    if compact_row != compact_output.transformed.n_rows() {
        return Err(invalid_dimensions(
            "expanded rlog non-zero rows",
            compact_row,
            compact_output.transformed.n_rows(),
        ));
    }
    Ok(RlogOutput {
        transformed: RowMajorMatrix::from_row_major(all_zero.len(), n_samples, values)?,
        intercept: expand_rlog_intercepts_with_all_zero_rows(&compact_output, all_zero)?,
        sample_prior_variance: compact_output.sample_prior_variance,
        offset_mode: compact_output.offset_mode,
    })
}

fn expand_rlog_intercepts_with_all_zero_rows(
    compact_output: &RlogOutput,
    all_zero: &[bool],
) -> Result<Vec<f64>, DeseqError> {
    let expected_nonzero = all_zero.iter().filter(|is_zero| !**is_zero).count();
    if compact_output.intercept.len() != expected_nonzero {
        return Err(invalid_dimensions(
            "expanded rlog intercepts",
            expected_nonzero,
            compact_output.intercept.len(),
        ));
    }
    let mut compact_row = 0usize;
    let mut intercept = Vec::with_capacity(all_zero.len());
    for is_zero in all_zero {
        if *is_zero {
            intercept.push(0.0);
        } else {
            intercept.push(compact_output.intercept[compact_row]);
            compact_row += 1;
        }
    }
    Ok(intercept)
}

fn expand_frozen_rlog_output_with_all_zero_rows(
    compact_output: RlogOutput,
    all_zero: &[bool],
    n_samples: usize,
    frozen_intercept: &[f64],
) -> Result<RlogOutput, DeseqError> {
    if frozen_intercept.len() != all_zero.len() {
        return Err(invalid_dimensions(
            "expanded frozen rlog intercepts",
            all_zero.len(),
            frozen_intercept.len(),
        ));
    }
    if compact_output.transformed.n_cols() != n_samples {
        return Err(invalid_dimensions(
            "expanded frozen rlog columns",
            n_samples,
            compact_output.transformed.n_cols(),
        ));
    }
    let mut values = Vec::with_capacity(all_zero.len() * n_samples);
    let mut compact_row = 0usize;
    for (gene, is_zero) in all_zero.iter().enumerate() {
        if *is_zero {
            values.extend(std::iter::repeat(frozen_intercept[gene]).take(n_samples));
        } else {
            let src = compact_output.transformed.row(compact_row)?;
            values.extend_from_slice(src);
            compact_row += 1;
        }
    }
    if compact_row != compact_output.transformed.n_rows() {
        return Err(invalid_dimensions(
            "expanded frozen rlog non-zero rows",
            compact_row,
            compact_output.transformed.n_rows(),
        ));
    }
    Ok(RlogOutput {
        transformed: RowMajorMatrix::from_row_major(all_zero.len(), n_samples, values)?,
        intercept: frozen_intercept.to_vec(),
        sample_prior_variance: compact_output.sample_prior_variance,
        offset_mode: compact_output.offset_mode,
    })
}

fn all_zero_glm_fit(
    counts: &CountMatrix,
    design: &DesignMatrix,
) -> Result<NbinomGlmFit, DeseqError> {
    Ok(NbinomGlmFit {
        log_like: vec![f64::NAN; counts.n_genes()],
        beta_converged: vec![false; counts.n_genes()],
        beta: RowMajorMatrix::from_elem(counts.n_genes(), design.n_coefficients(), f64::NAN)?,
        beta_se: RowMajorMatrix::from_elem(counts.n_genes(), design.n_coefficients(), f64::NAN)?,
        beta_covariance: Some(RowMajorMatrix::from_elem(
            counts.n_genes(),
            design.n_coefficients() * design.n_coefficients(),
            f64::NAN,
        )?),
        mu: RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), f64::NAN)?,
        beta_iter: vec![0; counts.n_genes()],
        model_matrix: design.clone(),
        n_terms: design.n_coefficients(),
        hat_diagonal: RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), f64::NAN)?,
    })
}

fn expand_glm_fit(
    compact_fit: NbinomGlmFit,
    all_zero: &[bool],
) -> Result<NbinomGlmFit, DeseqError> {
    Ok(NbinomGlmFit {
        log_like: expand_gene_values_with_nan_rows(&compact_fit.log_like, all_zero)?,
        beta_converged: expand_gene_values_with_fill_rows(
            &compact_fit.beta_converged,
            all_zero,
            false,
        )?,
        beta: expand_matrix_with_nan_rows(&compact_fit.beta, all_zero)?,
        beta_se: expand_matrix_with_nan_rows(&compact_fit.beta_se, all_zero)?,
        beta_covariance: compact_fit
            .beta_covariance
            .as_ref()
            .map(|matrix| expand_matrix_with_nan_rows(matrix, all_zero))
            .transpose()?,
        mu: expand_matrix_with_nan_rows(&compact_fit.mu, all_zero)?,
        beta_iter: expand_gene_values_with_fill_rows(&compact_fit.beta_iter, all_zero, 0)?,
        model_matrix: compact_fit.model_matrix,
        n_terms: compact_fit.n_terms,
        hat_diagonal: expand_matrix_with_nan_rows(&compact_fit.hat_diagonal, all_zero)?,
    })
}

fn mask_wald_degrees_of_freedom_for_all_zero_rows(
    wald: &mut WaldOutput,
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let Some(degrees_of_freedom) = &mut wald.degrees_of_freedom else {
        return Ok(());
    };
    if degrees_of_freedom.len() != all_zero.len() {
        return Err(invalid_dimensions(
            "Wald degrees of freedom all-zero mask",
            all_zero.len(),
            degrees_of_freedom.len(),
        ));
    }
    for (df, is_all_zero) in degrees_of_freedom.iter_mut().zip(all_zero.iter().copied()) {
        if is_all_zero {
            *df = None;
        }
    }
    Ok(())
}

fn apply_contrast_all_zero_to_wald_contrast(
    contrast: &mut WaldContrastOutput,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let n_genes = contrast.log2_fold_change.len();
    if contrast_all_zero.len() != n_genes {
        return Err(invalid_dimensions(
            "contrastAllZero rows",
            n_genes,
            contrast_all_zero.len(),
        ));
    }
    if all_zero.len() != n_genes {
        return Err(invalid_dimensions("allZero rows", n_genes, all_zero.len()));
    }
    for gene in 0..n_genes {
        if contrast_all_zero[gene] && !all_zero[gene] {
            contrast.log2_fold_change[gene] = Some(0.0);
            contrast.wald.stat[gene] = Some(0.0);
            contrast.wald.pvalue[gene] = Some(1.0);
        }
    }
    Ok(())
}

fn apply_contrast_all_zero_to_lrt_results(
    results: &mut DeseqResults,
    contrast_all_zero: &[bool],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    let n_genes = results.rows.len();
    if contrast_all_zero.len() != n_genes {
        return Err(invalid_dimensions(
            "contrastAllZero rows",
            n_genes,
            contrast_all_zero.len(),
        ));
    }
    if all_zero.len() != n_genes {
        return Err(invalid_dimensions("allZero rows", n_genes, all_zero.len()));
    }
    for gene in 0..n_genes {
        if contrast_all_zero[gene] && !all_zero[gene] {
            results.rows[gene].log2_fold_change = Some(0.0);
        }
    }
    Ok(())
}

fn attach_glm_fit(fit: &mut DeseqFit, glm_fit: NbinomGlmFit) {
    let full_deviance = glm_fit
        .log_like
        .iter()
        .map(|log_like| full_deviance_from_log_like(*log_like))
        .collect();
    fit.beta = Some(glm_fit.beta);
    fit.beta_se = Some(glm_fit.beta_se);
    fit.beta_covariance = glm_fit.beta_covariance;
    fit.beta_converged = Some(glm_fit.beta_converged);
    fit.beta_iter = Some(glm_fit.beta_iter);
    fit.log_like = Some(glm_fit.log_like);
    fit.full_deviance = Some(full_deviance);
    fit.mu = Some(glm_fit.mu);
    fit.hat_diagonal = Some(glm_fit.hat_diagonal);
}

fn full_deviance_from_log_like(log_like: f64) -> f64 {
    checked_product2(-2.0, log_like).unwrap_or(f64::NAN)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let deviance = left * right;
    if left.is_finite() && right.is_finite() && deviance.is_finite() {
        Some(deviance)
    } else {
        None
    }
}

fn validate_pipeline_wald_coefficient(
    design: &DesignMatrix,
    coefficient: usize,
) -> Result<(), DeseqError> {
    if coefficient >= design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "pipeline Wald coefficient index".to_string(),
            expected: design.n_coefficients().saturating_sub(1),
            actual: coefficient,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{full_deviance_from_log_like, CountMatrix};

    #[test]
    fn count_matrix_rejects_bad_length() {
        let err = CountMatrix::from_row_major_u32(2, 3, vec![1, 2]).unwrap_err();
        assert!(err.to_string().contains("invalid dimensions"));
    }

    #[test]
    fn count_matrix_row_access() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 2, 3, 4, 5, 6]).unwrap();
        assert_eq!(counts.row(1).unwrap(), &[4, 5, 6]);
    }

    #[test]
    fn count_matrix_accepts_u64_when_values_fit() {
        let counts = CountMatrix::from_row_major_u64(1, 3, vec![1, 2, 3]).unwrap();
        assert_eq!(counts.as_slice(), &[1, 2, 3]);
    }

    #[test]
    fn all_zero_gene_detection() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![0, 0, 0, 4, 5, 6]).unwrap();
        assert!(counts.is_all_zero_gene(0).unwrap());
        assert!(!counts.is_all_zero_gene(1).unwrap());
    }

    #[test]
    fn full_deviance_from_log_like_masks_nonfinite_values() {
        assert_eq!(full_deviance_from_log_like(-2.0), 4.0);
        assert!(full_deviance_from_log_like(f64::NAN).is_nan());
        assert!(full_deviance_from_log_like(f64::MAX).is_nan());
        assert!(full_deviance_from_log_like(-f64::MAX).is_nan());
    }
}
