//! Validated genes x samples count storage.

use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::normalize_index_range;

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
            && names.len() != n_genes
        {
            return Err(invalid_dimensions("gene names", n_genes, names.len()));
        }
        if let Some(names) = &sample_names
            && names.len() != n_samples
        {
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
