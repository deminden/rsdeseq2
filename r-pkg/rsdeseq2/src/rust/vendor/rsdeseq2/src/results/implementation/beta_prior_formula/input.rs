//! Data-only adapters; formula parsing and normalization policy stay separate.
use super::*;

pub(super) struct FormulaPriorInput<'a> {
    pub(super) counts: &'a CountMatrix,
    pub(super) size_factors: &'a [f64],
    pub(super) weights: Option<&'a RowMajorMatrix<f64>>,
    pub(super) dispersions: &'a [f64],
    pub(super) base_mean: &'a [f64],
    pub(super) disp_fit: &'a [f64],
    pub(super) gene_names: Option<&'a [String]>,
    pub(super) options: BetaPriorRefitOptions,
}

impl<'a> From<ExpandedFormulaBetaPriorWaldResultsInput<'a>> for FormulaPriorInput<'a> {
    fn from(input: ExpandedFormulaBetaPriorWaldResultsInput<'a>) -> Self {
        Self {
            counts: input.counts,
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        }
    }
}

impl<'a> From<ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'a>> for FormulaPriorInput<'a> {
    fn from(input: ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'a>) -> Self {
        Self {
            counts: input.counts,
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        }
    }
}

pub(super) struct FormulaPriorNormalizedInput<'a> {
    pub(super) counts: &'a CountMatrix,
    pub(super) normalization_factors: &'a RowMajorMatrix<f64>,
    pub(super) weights: Option<&'a RowMajorMatrix<f64>>,
    pub(super) dispersions: &'a [f64],
    pub(super) base_mean: &'a [f64],
    pub(super) disp_fit: &'a [f64],
    pub(super) gene_names: Option<&'a [String]>,
    pub(super) options: BetaPriorRefitOptions,
}

impl<'a> From<ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'a>>
    for FormulaPriorNormalizedInput<'a>
{
    fn from(input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'a>) -> Self {
        Self {
            counts: input.counts,
            normalization_factors: input.normalization_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        }
    }
}

impl<'a> From<ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'a>>
    for FormulaPriorNormalizedInput<'a>
{
    fn from(input: ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'a>) -> Self {
        Self {
            counts: input.counts,
            normalization_factors: input.normalization_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        }
    }
}
