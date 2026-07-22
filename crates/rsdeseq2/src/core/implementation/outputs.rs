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
    /// Test type selected for this replacement-refit output.
    pub fn test_type(&self) -> TestType {
        match self {
            Self::Wald(_) => TestType::Wald,
            Self::Lrt(_) => TestType::Lrt,
        }
    }

    /// Final merged result rows after replacement/refit.
    pub fn results(&self) -> &DeseqResults {
        match self {
            Self::Wald(output) => &output.results,
            Self::Lrt(output) => &output.results,
        }
    }

    /// Original result rows before replacement/refit.
    pub fn original_results(&self) -> &DeseqResults {
        match self {
            Self::Wald(output) => &output.original_results,
            Self::Lrt(output) => &output.original_results,
        }
    }

    /// Original fit before replacement/refit.
    pub fn original_fit(&self) -> &DeseqFit {
        match self {
            Self::Wald(output) => &output.original_fit,
            Self::Lrt(output) => &output.original_fit,
        }
    }

    /// Refit on replacement counts, when any non-all-zero replacement row exists.
    pub fn refit_fit(&self) -> Option<&DeseqFit> {
        match self {
            Self::Wald(output) => output.refit_fit.as_ref(),
            Self::Lrt(output) => output.refit_fit.as_ref(),
        }
    }

    /// Result rows from the replacement-count refit, before merge.
    pub fn refit_results(&self) -> Option<&DeseqResults> {
        match self {
            Self::Wald(output) => output.refit_results.as_ref(),
            Self::Lrt(output) => output.refit_results.as_ref(),
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

/// Output from the fast-VST GLM-mu helper.
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
    /// The trend was fit on the deterministic fast-VST subset.
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
    /// Metadata for explicit fast-VST diagnostics.
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
    /// Metadata for automatic VST diagnostics.
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
    /// Stable label for the automatic VST trend source.
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
