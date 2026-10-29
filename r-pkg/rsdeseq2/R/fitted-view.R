# Validate mutable dataset input before using the retained native fit.
.rsdeseq2_fitted_view <- function(object) {
  fit <- S4Vectors::metadata(object)$rsdeseq2RustFit
  if (is.null(fit) || !inherits(fit, "rsdeseq2RustFit")) {
    stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before using fitted state", call. = FALSE)
  }
  invalid <- function(field) stop(sprintf("Fitted dataset %s changed; refit with rsdeseq2::DESeq()", field), call. = FALSE)
  counts <- SummarizedExperiment::assay(object, "counts")
  if (!identical(colnames(counts), colnames(fit$counts)) || ncol(counts) != ncol(fit$counts)) invalid("samples")
  if (is.null(rownames(counts)) || is.null(rownames(fit$counts)) || anyDuplicated(rownames(counts))) invalid("gene identities")
  rows <- match(rownames(counts), rownames(fit$counts))
  if (anyNA(rows) || !length(rows)) invalid("genes")
  equal <- function(a, b) {
    identical(dim(a), dim(b)) && length(a) == length(b) &&
      isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = 0, check.attributes = FALSE))
  }
  if (!equal(counts, fit$counts[rows, , drop = FALSE])) invalid("counts")
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  frame <- stats::model.frame(methods::slot(object, "design"), data = cd)
  design <- stats::model.matrix(methods::slot(object, "design"), data = frame)
  if (!identical(dimnames(design), dimnames(fit$design)) || !equal(design, fit$design)) invalid("design")
  state <- S4Vectors::metadata(object)$rsdeseq2Normalization
  offsets <- stats::model.offset(frame)
  if (is.null(state) && !is.null(offsets)) invalid("offset provenance")
  if (!is.null(state) && !equal(offsets, state$offsets)) invalid("offsets")
  nf <- DESeq2::normalizationFactors(object)
  storedNF <- if (is.null(fit$normalizationFactors)) NULL else fit$normalizationFactors[rows, , drop = FALSE]
  if (!equal(nf, storedNF)) invalid("normalization factors")
  if (!equal(DESeq2::sizeFactors(object), fit$sizeFactors)) invalid("size factors")
  weights <- if ("weights" %in% SummarizedExperiment::assayNames(object)) SummarizedExperiment::assay(object, "weights") else NULL
  storedWeights <- if (is.null(fit$weights)) NULL else fit$weights[rows, , drop = FALSE]
  if (!equal(weights, storedWeights)) invalid("weights")
  if (!equal(S4Vectors::mcols(object)$dispersion, fit$dispersions[rows])) invalid("dispersions")
  for (name in fit$resultsNames) {
    j <- match(name, fit$resultsNames)
    if (!equal(S4Vectors::mcols(object)[[paste0("MLE_", name)]], fit$beta[rows, j])) invalid("coefficients")
    if (!equal(S4Vectors::mcols(object)[[paste0("SE_", name)]], fit$betaSE[rows, j])) invalid("standard errors")
  }
  for (name in c("mu", "cooks")) {
    if (!equal(SummarizedExperiment::assay(object, name), fit[[name]][rows, , drop = FALSE])) invalid(name)
  }
  for (name in c("dispFit", "baseMean", "baseVar", "allZero")) {
    if (!equal(S4Vectors::mcols(object)[[name]], fit[[name]][rows])) invalid(name)
  }
  for (name in c("betaConv", "maxCooks")) {
    stored <- if (name == "betaConv") fit$results$converged else fit$results$maxCooks
    if (!equal(S4Vectors::mcols(object)[[name]], stored[rows])) invalid(name)
  }
  eligible <- any(fit$replaceable)
  if (!equal(S4Vectors::mcols(object)$replace,
             if (eligible) fit$replace[rows] else NULL)) invalid("replacement flags")
  if (!equal(SummarizedExperiment::colData(object)$replaceable,
             if (eligible) fit$replaceable else NULL)) invalid("replacement eligibility")
  for (name in c("replaceCounts", "replaceCooks")) {
    current <- if (name %in% SummarizedExperiment::assayNames(object)) SummarizedExperiment::assay(object, name) else NULL
    stored <- if (fit$nRefit > 0L) {
      if (name == "replaceCounts") fit$replaceCounts[rows, , drop = FALSE] else fit$cooks[rows, , drop = FALSE]
    } else NULL
    if (!equal(current, stored)) invalid(name)
  }
  vectors <- c("baseMean", "baseVar", "allZero", "replace", "refitRows", "newAllZeroRows",
               "dispersions", "dispFit", "rawPvalue", "rawStat", "lrtStat", "lrtPvalue")
  matrices <- c("beta", "betaSE", "mu", "cooks", "replaceCounts", "normalizationFactors", "weights", "counts")
  for (name in vectors) if (!is.null(fit[[name]])) fit[[name]] <- fit[[name]][rows]
  for (name in matrices) if (!is.null(fit[[name]])) fit[[name]] <- fit[[name]][rows, , drop = FALSE]
  fit$betaCovariance <- fit$betaCovariance[rows, , , drop = FALSE]
  fit$results <- fit$results[rows, , drop = FALSE]
  fit$nRefit <- sum(fit$refitRows | fit$newAllZeroRows, na.rm = TRUE)
  fit$colData <- cd
  fit
}
