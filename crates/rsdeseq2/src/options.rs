use serde::{Deserialize, Serialize};

use crate::errors::{invalid_dimensions, DeseqError};

/// Size-factor estimation method.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum SizeFactorMethod {
    /// DESeq2 median-ratio method.
    #[default]
    Ratio,
    /// DESeq2 `poscounts` method, using the nth root of the product of positive counts.
    PosCounts,
}

/// Rust representation of DESeq2's `controlGenes` argument.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ControlGenes {
    /// Zero-based row indices.
    Indices(Vec<usize>),
    /// Logical mask with one value per gene.
    Mask(Vec<bool>),
}

impl ControlGenes {
    /// Convert control genes to zero-based indices after validating dimensions.
    pub fn to_indices(&self, n_genes: usize) -> Result<Vec<usize>, DeseqError> {
        match self {
            Self::Indices(indices) => {
                if indices.is_empty() {
                    return Err(DeseqError::NoUsableGenesForSizeFactors);
                }
                for index in indices {
                    if *index >= n_genes {
                        return Err(DeseqError::InvalidDimensions {
                            context: "control gene index".to_string(),
                            expected: n_genes.saturating_sub(1),
                            actual: *index,
                        });
                    }
                }
                Ok(indices.clone())
            }
            Self::Mask(mask) => {
                if mask.len() != n_genes {
                    return Err(invalid_dimensions("control gene mask", n_genes, mask.len()));
                }
                let indices = mask
                    .iter()
                    .copied()
                    .enumerate()
                    .filter_map(|(idx, keep)| keep.then_some(idx))
                    .collect::<Vec<_>>();
                if indices.is_empty() {
                    return Err(DeseqError::NoUsableGenesForSizeFactors);
                }
                Ok(indices)
            }
        }
    }
}

/// Options for DESeq2-like size-factor estimation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SizeFactorOptions {
    /// Size-factor estimator.
    pub method: SizeFactorMethod,
    /// Optional caller-supplied size factors.
    pub supplied_size_factors: Option<Vec<f64>>,
    /// Optional supplied per-gene geometric means.
    pub geo_means: Option<Vec<f64>>,
    /// Optional control-gene subset.
    pub control_genes: Option<ControlGenes>,
}

/// Cook's distance p-value filtering option for result assembly.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum CooksCutoff {
    /// Use DESeq2's default `qf(.99, p, m - p)` cutoff.
    #[default]
    Default,
    /// Disable Cook's p-value filtering.
    Disabled,
    /// Use an explicit Cook's distance threshold.
    Threshold(f64),
}

impl SizeFactorOptions {
    /// Create options for a method.
    pub fn new(method: SizeFactorMethod) -> Self {
        Self {
            method,
            supplied_size_factors: None,
            geo_means: None,
            control_genes: None,
        }
    }
}

/// Dispersion trend fit type.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum FitType {
    /// Parametric dispersion trend.
    #[default]
    Parametric,
    /// Local dispersion trend.
    Local,
    /// Mean dispersion fit.
    Mean,
    /// Future glmGamPoi-like mode.
    GlmGamPoi,
}

/// Test type for the future full workflow.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum TestType {
    /// Wald test.
    #[default]
    Wald,
    /// Likelihood-ratio test.
    Lrt,
}

/// Execution policy for deterministic parity, speed, or diagnostics.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum ExecutionMode {
    /// Deterministic, parity-first behavior.
    #[default]
    Strict,
    /// Future aggressive parallel and approximate paths.
    Fast,
    /// Expose extra intermediate state and convergence diagnostics.
    Diagnostic,
}
