use crate::all_zero::{
    expand_gene_values_with_fill_rows, expand_gene_values_with_nan_rows,
    expand_matrix_with_nan_rows, mask_all_zero_values_with_nan_rows, nonzero_gene_indices,
};
use crate::contrasts::{
    contrast_all_zero_factor_levels, contrast_all_zero_numeric, resolve_coefficient_index,
    resolve_contrast, resolve_results_contrast, factor_level_contrast_from_model_frame,
    ContrastSpec, FactorLevelContrast, ResultsContrast,
};
use crate::cooks::{
    calculate_cooks_distance, prepare_cooks_replacement_refit, CooksOutput, CooksRefitPlan,
    CooksReplacementOptions,
};
use crate::design::{r_like_name_candidates as design_name_candidates, DesignMatrix, FormulaModelFrame};
use crate::design::{
    expanded_formula_design_from_model_frame, expanded_formula_design_with_offsets_from_model_frame,
    formula_has_offset_terms, ExpandedAdditiveFactorDesign, ExpandedFormulaDesignWithOffsets,
};
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
use crate::matrix::{normalize_index_range, RowMajorMatrix};
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
    apply_cooks_cutoff, apply_cooks_cutoff_with_low_count_heuristic,
    build_lrt_contrast_results, build_lrt_results, build_wald_contrast_results,
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
        if let Some(names) = &gene_names
            && names.len() != n_genes {
                return Err(invalid_dimensions("gene names", n_genes, names.len()));
            }
        if let Some(names) = &sample_names
            && names.len() != n_samples {
                return Err(invalid_dimensions("sample names", n_samples, names.len()));
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

    /// Reusable gene-index span.
    pub fn gene_indices(&self) -> core::range::Range<usize> {
        core::range::Range {
            start: 0,
            end: self.n_genes,
        }
    }

    /// Reusable sample-index span.
    pub fn sample_indices(&self) -> core::range::Range<usize> {
        core::range::Range {
            start: 0,
            end: self.n_samples,
        }
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

    /// Return a contiguous block of gene rows in row-major order.
    ///
    /// The range accepts both legacy range syntax (`1..3`) and the newer
    /// `core::range` types. The returned slice contains
    /// `n_genes_in_range * n_samples` counts.
    pub fn gene_rows<R: core::ops::RangeBounds<usize>>(
        &self,
        genes: R,
    ) -> Result<&[u32], DeseqError> {
        let (start_gene, end_gene) = normalize_index_range(genes, self.n_genes, "gene range")?;
        let start = start_gene.checked_mul(self.n_samples).ok_or_else(|| {
            DeseqError::InvalidDimensions {
                context: "gene range start overflow".to_string(),
                expected: self.counts.len(),
                actual: usize::MAX,
            }
        })?;
        let end =
            end_gene
                .checked_mul(self.n_samples)
                .ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "gene range end overflow".to_string(),
                    expected: self.counts.len(),
                    actual: usize::MAX,
                })?;
        Ok(&self.counts[start..end])
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
        self.gene_indices()
            .into_iter()
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
        let all_zero_genes = self
            .gene_indices()
            .into_iter()
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
