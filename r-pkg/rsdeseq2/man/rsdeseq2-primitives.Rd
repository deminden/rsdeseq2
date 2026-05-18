\name{rsdeseq2-primitives}
\alias{rsdeseq2-primitives}
\alias{DESeq}
\alias{results}
\alias{estimateSizeFactorsRust}
\alias{normalizedCountsRust}
\alias{baseMeanRust}
\alias{baseMetadataRust}
\alias{applyCooksCutoffRust}
\alias{resultsTableRust}
\alias{applyIndependentFilteringRust}
\alias{deseq2McolsDiagnosticsRust}
\alias{rsdeseq2DiagnosticSchemaRust}
\alias{estimateDispersionsRust}
\alias{nbinomWaldTestRust}
\alias{nbinomLRTRust}
\title{Primitive rsdeseq2 R wrapper helpers}
\description{
Early R-facing helpers for implemented rsdeseq2 primitive matrix stages.
These functions are not a DESeq2 drop-in replacement and do not implement
DESeqDataSet integration yet.
}
\usage{
estimateSizeFactorsRust(counts, method = c("ratio", "poscounts"),
  geoMeans = NULL, controlGenes = NULL, native = FALSE)

normalizedCountsRust(counts, sizeFactors = NULL,
  normalizationFactors = NULL, native = FALSE)

baseMeanRust(counts, sizeFactors = NULL, normalizationFactors = NULL,
  native = FALSE)

baseMetadataRust(counts, sizeFactors = NULL, normalizationFactors = NULL,
  weights = NULL, native = FALSE)

applyCooksCutoffRust(pvalue, maxCooks, cooksCutoff, counts = NULL,
  cooks = NULL, lowCountHeuristic = FALSE, native = FALSE)

resultsTableRust(baseMean, log2FoldChange = NULL, lfcSE = NULL,
  stat = NULL, pvalue = NULL, padj = NULL, dispersion = NULL,
  converged = NULL, rowNames = NULL)

applyIndependentFilteringRust(results, alpha = 0.1, theta = NULL,
  enabled = TRUE)

deseq2McolsDiagnosticsRust(nGenes, test = c("none", "Wald", "LRT"),
  rowNames = NULL, betaConv = NULL, fullBetaConv = NULL,
  reducedBetaConv = NULL, betaIter = NULL, reducedBetaIter = NULL,
  deviance = NULL, maxCooks = NULL)

rsdeseq2DiagnosticSchemaRust(native = FALSE)

DESeq(object, fitType = "parametric", test = "Wald", engine = "rust",
  parityMode = "strict", nproc = 1, ...)

results(object, ...)

estimateDispersionsRust(...)

nbinomWaldTestRust(...)

nbinomLRTRust(...)
}
\arguments{
\item{counts}{Numeric or integer matrix with genes in rows and samples in columns.}
\item{method}{Size-factor method, either \code{"ratio"} or \code{"poscounts"}.}
\item{geoMeans}{Optional supplied geometric means, one per gene.}
\item{controlGenes}{Optional control-gene subset, as logical, numeric, or row-name values.}
\item{native}{Whether to try the registered native \code{.Call} bridge before falling back to R.}
\item{sizeFactors}{Positive numeric size factors, one per sample.}
\item{normalizationFactors}{Positive numeric matrix with the same dimensions as \code{counts}.}
\item{weights}{Optional finite non-negative observation-weight matrix.}
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
\item{object}{Future DESeqDataSet-like input. Currently unsupported.}
\item{fitType, engine, parityMode, nproc}{Future high-level DESeq options.}
\item{...}{Additional arguments reserved for future compatibility.}
}
\details{
The implemented helpers operate on primitive R vectors and matrices. When
\code{normalizationFactors} are supplied, they preempt \code{sizeFactors},
matching DESeq2 normalized-count behavior. \code{resultsTableRust} assembles
already-computed vectors into a DESeq2-shaped primitive result table but does
not run Wald, LRT, contrast handling, outlier replacement, or independent
filtering itself. \code{applyIndependentFilteringRust} applies baseMean-driven
filtered BH adjustment to such primitive result tables and stores filtering
metadata in the \code{"independentFiltering"} attribute. Unsupported high-level
entry points raise explicit errors.
}
\value{
Depending on the helper, a numeric vector, numeric matrix, data frame, or an
explicit unsupported-feature error.
}
