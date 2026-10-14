\name{rsdeseq2-primitives}
\alias{rsdeseq2-primitives}
\alias{DESeq}
\alias{fitRust}
\alias{results}
\alias{estimateSizeFactorsRust}
\alias{normalizedCountsRust}
\alias{baseMeanRust}
\alias{baseMetadataRust}
\alias{applyCooksCutoffRust}
\alias{resultsTableRust}
\alias{waldFitRust}
\alias{lrtFitRust}
\alias{resultsNames}
\alias{resultsNamesRust}
\alias{resolveResultsContrastRust}
\alias{applyIndependentFilteringRust}
\alias{deseq2McolsDiagnosticsRust}
\alias{rsdeseq2DiagnosticSchemaRust}
\alias{estimateDispersionsRust}
\alias{nbinomWaldTestRust}
\alias{nbinomLRTRust}
\alias{vst}
\alias{varianceStabilizingTransformation}
\alias{rlog}
\alias{lfcShrink}
\title{rsdeseq2 R wrapper helpers}
\description{
R-facing entry points to the Rust rsdeseq2 core for matrix stages, Wald/LRT
fitting, and DESeq2-compatible result-table semantics. The package builds a
static library from the vendored Rust source using Cargo.
}
\usage{
DESeq(object, fitType = c("parametric", "local", "mean"),
  test = c("Wald", "LRT"), engine = "rust", parityMode = "strict",
  nproc = 1L, design = NULL, reduced = NULL, sizeFactors = NULL,
  normalizationFactors = NULL, weights = NULL, ...)

fitRust(object, fitType = c("parametric", "local", "mean"),
  test = c("Wald", "LRT"), engine = "rust", parityMode = "strict",
  nproc = 1L, design = NULL, reduced = NULL, sizeFactors = NULL,
  normalizationFactors = NULL, weights = NULL, ...)

estimateSizeFactorsRust(counts, method = c("ratio", "poscounts"),
  geoMeans = NULL, controlGenes = NULL, native = TRUE)

normalizedCountsRust(counts, sizeFactors = NULL,
  normalizationFactors = NULL, native = TRUE)

baseMeanRust(counts, sizeFactors = NULL, normalizationFactors = NULL,
  native = TRUE)

baseMetadataRust(counts, sizeFactors = NULL, normalizationFactors = NULL,
  weights = NULL, native = TRUE)

applyCooksCutoffRust(pvalue, maxCooks, cooksCutoff, counts = NULL,
  cooks = NULL, lowCountHeuristic = FALSE, native = TRUE)

resultsTableRust(baseMean, log2FoldChange = NULL, lfcSE = NULL,
  stat = NULL, pvalue = NULL, padj = NULL, dispersion = NULL,
  converged = NULL, contrast = NULL, rowNames = NULL)

waldFitRust(baseMean, beta, betaCovariance = NULL, betaSE = NULL,
  resultsNames = NULL, counts = NULL, sampleLevels = NULL,
  factorLevels = NULL, factorReferences = NULL, modelMatrix = NULL,
  colData = NULL, rowNames = NULL)

lrtFitRust(baseMean, beta, betaCovariance = NULL, betaSE = NULL,
  lrtStat, lrtPvalue = NULL, lrtDf = NULL, resultsNames = NULL,
  counts = NULL, sampleLevels = NULL, factorLevels = NULL,
  factorReferences = NULL, modelMatrix = NULL, colData = NULL,
  rowNames = NULL)

resultsNames(object, ...)

resultsNamesRust(object)

resolveResultsContrastRust(contrast, resultsNames, listValues = c(1, -1),
  reference = NULL)

applyIndependentFilteringRust(results, alpha = 0.1, theta = NULL,
  enabled = TRUE)

deseq2McolsDiagnosticsRust(nGenes, test = c("none", "Wald", "LRT"),
  rowNames = NULL, betaConv = NULL, fullBetaConv = NULL,
  reducedBetaConv = NULL, betaIter = NULL, reducedBetaIter = NULL,
  deviance = NULL, maxCooks = NULL)

rsdeseq2DiagnosticSchemaRust(native = FALSE)

DESeq(object, fitType = c("parametric", "local", "mean"),
  test = c("Wald", "LRT"), engine = "rust",
  parityMode = "strict", nproc = 1, design = NULL, reduced = NULL,
  sizeFactors = NULL, normalizationFactors = NULL, weights = NULL, ...)

results(object, contrast = NULL, name = NULL, listValues = c(1, -1),
  reference = NULL, lfcThreshold = 0,
  altHypothesis = c("greaterAbs", "lessAbs", "greater", "less",
  "greaterAbs2014", "greaterAbsUPSHOT"), useT = FALSE,
  degreesOfFreedom = NULL, alpha = 0.1,
  independentFiltering = TRUE, ...)

estimateDispersionsRust(counts, design,
  fitType = c("parametric", "local", "mean"), nproc = 1L)

nbinomWaldTestRust(baseMean = NULL, beta = NULL, betaCovariance = NULL,
  betaSE = NULL, resultsNames = NULL, counts = NULL,
  sampleLevels = NULL, factorLevels = NULL, factorReferences = NULL,
  modelMatrix = NULL, colData = NULL, rowNames = NULL, object = NULL, ...)

nbinomLRTRust(baseMean = NULL, beta = NULL, betaCovariance = NULL,
  betaSE = NULL, lrtStat = NULL, lrtPvalue = NULL, lrtDf = NULL,
  resultsNames = NULL, counts = NULL, sampleLevels = NULL, factorLevels = NULL,
  factorReferences = NULL, modelMatrix = NULL, colData = NULL,
  rowNames = NULL, object = NULL, ...)

vst(object, blind = TRUE, nsub = 1000L,
  fitType = c("parametric", "local", "mean"), nproc = 1L)

varianceStabilizingTransformation(object, blind = TRUE,
  fitType = c("parametric", "local", "mean"), nproc = 1L)

rlog(object, blind = TRUE, intercept, betaPriorVar,
  fitType = c("parametric", "local", "mean"), nproc = 1L)

lfcShrink(dds, coef = NULL, type = "normal", res = NULL, contrast = NULL,
  lfcThreshold = 0, svalue = FALSE, independentFiltering = TRUE,
  alpha = 0.1, apeAdapt = TRUE, apeMethod = "nbinomCR",
  quiet = FALSE, ...)
}
\arguments{
\item{counts}{Numeric or integer matrix with genes in rows and samples in columns.}
\item{method}{Size-factor method, either \code{"ratio"} or \code{"poscounts"}.}
\item{geoMeans}{Optional supplied geometric means, one per gene.}
\item{controlGenes}{Optional control-gene subset, as logical, numeric, or row-name values.}
\item{native}{Whether to use the registered native \code{.Call} bridge when
available.}
\item{sizeFactors}{Positive numeric size factors, one per sample. In high-level
fitting, these are used when normalization factors are absent.}
\item{normalizationFactors}{Positive numeric matrix with the same dimensions as
\code{counts}; these preempt size factors in high-level fitting.}
\item{weights}{Optional finite non-negative observation-weight matrix for base
metadata or high-level fitting.}
\item{pvalue}{Numeric p-value vector with values in \code{[0, 1]} or \code{NA}.}
\item{maxCooks}{Numeric vector of maximum Cook's distances, one per p-value.}
\item{cooksCutoff}{Numeric cutoff, \code{FALSE}, or \code{NULL}.}
\item{cooks}{Optional Cook's distance matrix for the low-count heuristic.}
\item{lowCountHeuristic}{Whether to apply the explicit two-group low-count Cook's heuristic.}
\item{baseMean}{Numeric base mean vector, one per gene.}
\item{log2FoldChange, lfcSE, stat}{Optional numeric result vectors, one per gene.}
\item{padj}{Optional BH-adjusted p-values. If omitted and \code{pvalue} is supplied,
it is computed with \code{stats::p.adjust(..., method = "BH")}.}
\item{dispersion}{Optional numeric final dispersion vector, one per gene.}
\item{converged}{Optional logical convergence vector, one per gene.}
\item{beta}{Numeric matrix of log2-scale beta estimates, genes by coefficients.}
\item{betaCovariance}{Optional numeric array of log2-scale beta covariances,
genes by coefficients by coefficients. Required for multi-coefficient
contrast standard errors.}
\item{betaSE}{Optional numeric matrix of coefficient standard errors with the
same shape as \code{beta}.}
\item{lrtStat}{Non-negative LRT statistic vector, one per gene.}
\item{lrtPvalue}{Optional LRT p-value vector, one per gene. If omitted,
\code{lrtDf} is required and p-values are computed from chi-square tails.}
\item{lrtDf}{Optional scalar or per-gene positive LRT degrees of freedom.}
\item{resultsNames}{Coefficient names corresponding to \code{beta} columns.}
\item{sampleLevels}{Optional character vector of factor levels, one per sample,
used for character contrast all-zero handling.}
\item{factorLevels}{Optional named list of factor levels. For character
contrasts, the first level is used as the reference when \code{reference} is
not supplied, matching DESeq2 fitted-object behavior.}
\item{factorReferences}{Optional named character vector or list of factor
reference levels used for character contrast resolution.}
\item{modelMatrix}{Optional numeric model matrix, samples by coefficients, used
for numeric/list contrast all-zero handling.}
\item{colData}{Optional sample metadata used to derive factor levels and
character contrast sample-level vectors.}
\item{contrast}{DESeq2-style contrast specification: character triplet, list, or
numeric vector. For normal \code{lfcShrink}, numeric vectors refer to expanded
design coefficients and require a supplied result table.}
\item{name}{Optional coefficient name for \code{results()}.}
\item{listValues}{Two numeric weights for list contrasts.}
\item{reference}{Optional reference level for factor-level contrast resolution.}
\item{lfcThreshold}{Non-negative log2 fold-change threshold for thresholded
Wald tests.}
\item{altHypothesis}{DESeq2-style threshold alternative.}
\item{useT}{Whether to use Student t tails instead of Normal tails for
primitive Wald p-values.}
\item{degreesOfFreedom}{Scalar or per-gene t degrees of freedom when
\code{useT = TRUE}.}
\item{independentFiltering}{Whether \code{results()} applies the implemented
independent-filtering helper to primitive result objects.}
\item{results}{Primitive result data frame containing at least
\code{baseMean}, \code{pvalue}, and \code{padj}.}
\item{alpha}{Target FDR threshold for independent-filter rejection counting.}
\item{theta}{Optional independent-filter quantile grid. When omitted, the
DESeq2-style default grid starts at the fraction of zero \code{baseMean} values
and ends at \code{0.95}, or \code{1.0} if that fraction is at least
\code{0.95}.}
\item{enabled}{Whether to apply independent filtering. If false, regular BH
adjustment is recomputed and filtering metadata records the disabled state.}
\item{nGenes}{Number of genes for diagnostic metadata helpers.}
\item{test}{Diagnostic schema mode.}
\item{rowNames}{Optional row names for diagnostic metadata.}
\item{betaConv, fullBetaConv, reducedBetaConv}{Optional logical convergence vectors.}
\item{betaIter, reducedBetaIter}{Optional non-negative iteration-count vectors.}
\item{deviance}{Optional numeric deviance diagnostic vector.}
\item{object}{A primitive \code{waldFitRust()} or \code{lrtFitRust()} object
for implemented \code{results()} behavior, a list-like primitive Wald/LRT input,
or an already-fitted DESeq2-shaped list/S4-like object with compatible \code{mcols}
fields for direct \code{results()} dispatch, \code{nbinomWaldTestRust()}, or
\code{nbinomLRTRust()}.}
\item{fitType, engine, parityMode, nproc}{Trend fit, Rust engine, execution
mode (\code{strict}, \code{fast}, or \code{diagnostic}), and thread count.}
\item{blind}{Whether a transformation should ignore the dataset design and use
an intercept-only dispersion fit.}
\item{nsub}{Positive integer number of genes used to fit the fast VST trend.
Genes are selected by ordinary normalized-count means above 5, including when
observation weights are present. The transform errors if fewer than \code{nsub}
genes qualify.}
\item{intercept, betaPriorVar}{Optional frozen rlog gene intercepts on the
log2 scale and their shared sample-effect prior variance. Supply both together.}
\item{dds}{A DESeqDataSet fitted with \code{rsdeseq2::DESeq()}.}
\item{coef}{A one-based coefficient index or coefficient name from
\code{resultsNames(dds)}.}
\item{type}{Shrinkage method. Dataset workflows support \code{"normal"} and
optional \code{"ashr"} or \code{"apeglm"} methods.}
\item{res}{Optional result table containing \code{log2FoldChange} and
\code{lfcSE}.}
\item{svalue}{Whether to request s-values. Normal shrinkage requires
\code{FALSE}; \code{ashr} and \code{apeglm} support s-values.}
\item{apeAdapt}{Whether apeglm adapts its prior using the supplied result
estimate and standard error.}
\item{apeMethod}{Apeglm optimization method, such as \code{"nbinomCR"},
\code{"nbinomR"}, or \code{"general"}.}
\item{quiet}{Whether to suppress method messages.}
\item{design, reduced}{Full and reduced numeric model matrices, or formulas for
\code{DESeqDataSet} input.}
\item{...}{Additional options passed to the selected optional shrinkage
package. Normal shrinkage rejects unsupported arguments.}
}
\details{
The implemented helpers operate on primitive R vectors and matrices. When
\code{normalizationFactors} are supplied, they preempt \code{sizeFactors},
matching DESeq2 normalized-count behavior. \code{resultsTableRust} assembles
already-computed vectors into a DESeq2-shaped primitive result table but does
not run Wald, LRT, outlier replacement, or independent filtering itself.
\code{waldFitRust} packages already-computed Wald beta/covariance outputs.
\code{nbinomWaldTestRust} accepts those primitive fields, or a list-like
primitive object containing them, and returns the same primitive Wald fit
object.
\code{lrtFitRust} and \code{nbinomLRTRust} package already-computed primitive
LRT statistic/p-value diagnostics together with full-model beta/covariance
outputs. For primitive LRT fits, \code{results} reports the selected
coefficient or contrast LFC/SE while preserving the LRT statistic and p-value.
\code{results} can consume primitive Wald and LRT objects for selected
coefficients and DESeq2-style \code{contrast = ...} requests, including
thresholded
\code{lfcThreshold}/\code{altHypothesis} Wald p-values with Normal or supplied
t tails. It can also consume list-like objects containing the same primitive
Wald/LRT fields directly, including already-fitted DESeq2-shaped list/S4 metadata
with coefficient columns and matching \code{SE_} columns in \code{mcols}.
\code{resolveResultsContrastRust} resolves character triplet,
list/listValues, and numeric contrast forms while preserving the matching
character or numeric all-zero rule. Named coefficient results do not apply
contrast all-zero handling. \code{resultsNames} and
\code{resultsNamesRust} return coefficient names from primitive Wald/LRT fit
objects, list-like primitive objects, fitted-object metadata, beta matrices, or
supplied character vectors for contrast-resolution workflows.
\code{applyIndependentFilteringRust} applies baseMean-driven filtered BH
adjustment to primitive result tables and stores filtering metadata in the
\code{"independentFiltering"} attribute. \code{DESeq} fits matrix and list
inputs into an \code{rsdeseq2RustFit}; for a \code{DESeqDataSet}, it returns the
dataset with the Rust fit in metadata and fitted size or normalization factors.
The wrapper does not populate every DESeq2 diagnostic assay or implement every
DESeq2 fitting option.
\code{vst}, \code{varianceStabilizingTransformation}, and \code{rlog} accept
DESeqDataSet inputs and return DESeqTransform objects while Rust performs the
statistical transformation. \code{rlog} accepts frozen gene intercepts and a
shared sample-effect prior variance. \code{lfcShrink} implements normal
beta-prior shrinkage for a named or indexed Wald coefficient and expanded-model
factor-level, numeric, and list contrasts. Numeric contrasts require a supplied
result table; list contrasts average selected expanded coefficients and need
supplied results when their names are absent from the unshrunk design. The
contrast route uses a Rust expanded-model refit and returns a DESeqResults
object with prior information. Optional \code{ashr} shrinkage returns
posterior means/standard deviations and s-values. Optional \code{apeglm}
shrinkage supports named/indexed Wald coefficients and returns posterior
standard deviations and s-values. Formula interactions and transformed-factor
normal contrasts are not implemented.
}
\value{
Depending on the helper, a numeric vector, numeric matrix, data frame, or an
explicit unsupported-feature error.
}
