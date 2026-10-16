# Primitive Wald/LRT fits and result-name entry points.

waldFitRust <- function(baseMean,
                        beta,
                        betaCovariance = NULL,
                        betaSE = NULL,
                        resultsNames = NULL,
                        counts = NULL,
                        sampleLevels = NULL,
                        factorLevels = NULL,
                        factorReferences = NULL,
                        modelMatrix = NULL,
                        colData = NULL,
                        rowNames = NULL) {
  beta <- .rsdeseq2_validate_beta_matrix(beta, resultsNames)
  resultsNames <- colnames(beta)
  nGenes <- nrow(beta)
  geneNames <- rownames(beta)
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  baseMean <- .rsdeseq2_align_named_gene_vector(baseMean, geneNames, "baseMean")
  if (length(baseMean) != nGenes) {
    stop("baseMean must have one value per beta row", call. = FALSE)
  }
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% geneNames %||% names(baseMean), nGenes)
  betaCovariance <- .rsdeseq2_validate_beta_covariance(betaCovariance, nGenes, resultsNames)
  betaCovariance <- .rsdeseq2_align_named_gene_array(betaCovariance, geneNames, "betaCovariance")
  betaSE <- .rsdeseq2_validate_beta_se_matrix(betaSE, nGenes, resultsNames)
  betaSE <- .rsdeseq2_align_named_gene_matrix(betaSE, geneNames, "betaSE")
  counts <- .rsdeseq2_validate_optional_contrast_counts(counts, nGenes)
  counts <- .rsdeseq2_align_named_gene_matrix(counts, geneNames, "counts")
  sampleLevels <- .rsdeseq2_validate_optional_sample_levels(sampleLevels, counts)
  colData <- .rsdeseq2_validate_optional_col_data(colData, counts)
  factorLevels <- factorLevels %||% .rsdeseq2_factor_levels_from_col_data(colData)
  factorLevels <- .rsdeseq2_validate_optional_factor_levels(factorLevels)
  factorReferences <- .rsdeseq2_validate_optional_factor_references(factorReferences, factorLevels)
  modelMatrix <- .rsdeseq2_validate_optional_model_matrix(modelMatrix, counts, resultsNames)

  out <- list(
    baseMean = baseMean,
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    resultsNames = resultsNames,
    counts = counts,
    sampleLevels = sampleLevels,
    colData = colData,
    factorLevels = factorLevels,
    factorReferences = factorReferences,
    modelMatrix = modelMatrix,
    rowNames = rowNames
  )
  class(out) <- c("rsdeseq2PrimitiveWaldFit", "list")
  out
}

lrtFitRust <- function(baseMean,
                       beta,
                       betaCovariance = NULL,
                       betaSE = NULL,
                       lrtStat,
                       lrtPvalue = NULL,
                       lrtDf = NULL,
                       resultsNames = NULL,
                       counts = NULL,
                       sampleLevels = NULL,
                       factorLevels = NULL,
                       factorReferences = NULL,
                       modelMatrix = NULL,
                       colData = NULL,
                       rowNames = NULL) {
  beta <- .rsdeseq2_validate_beta_matrix(beta, resultsNames)
  resultsNames <- colnames(beta)
  nGenes <- nrow(beta)
  geneNames <- rownames(beta)
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  baseMean <- .rsdeseq2_align_named_gene_vector(baseMean, geneNames, "baseMean")
  if (length(baseMean) != nGenes) {
    stop("baseMean must have one value per beta row", call. = FALSE)
  }
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% geneNames %||% names(baseMean), nGenes)
  betaCovariance <- .rsdeseq2_validate_beta_covariance(betaCovariance, nGenes, resultsNames)
  betaCovariance <- .rsdeseq2_align_named_gene_array(betaCovariance, geneNames, "betaCovariance")
  betaSE <- .rsdeseq2_validate_beta_se_matrix(betaSE, nGenes, resultsNames)
  betaSE <- .rsdeseq2_align_named_gene_matrix(betaSE, geneNames, "betaSE")
  lrtStat <- .rsdeseq2_validate_lrt_stat(lrtStat, nGenes)
  lrtStat <- .rsdeseq2_align_named_gene_vector(lrtStat, geneNames, "lrtStat")
  lrtDf <- .rsdeseq2_align_named_gene_vector(lrtDf, geneNames, "lrtDf", allowScalar = TRUE)
  lrtDf <- .rsdeseq2_validate_lrt_df(lrtDf, nGenes)
  if (is.null(lrtPvalue)) {
    if (is.null(lrtDf)) {
      stop("lrtDf is required when lrtPvalue is NULL", call. = FALSE)
    }
    lrtPvalue <- stats::pchisq(lrtStat, df = lrtDf, lower.tail = FALSE)
    lrtPvalue[is.na(lrtStat) | is.na(lrtDf)] <- NA_real_
  } else {
    lrtPvalue <- .rsdeseq2_align_named_gene_vector(lrtPvalue, geneNames, "lrtPvalue")
    lrtPvalue <- .rsdeseq2_validate_result_pvalues(lrtPvalue, nGenes)
  }
  counts <- .rsdeseq2_validate_optional_contrast_counts(counts, nGenes)
  counts <- .rsdeseq2_align_named_gene_matrix(counts, geneNames, "counts")
  sampleLevels <- .rsdeseq2_validate_optional_sample_levels(sampleLevels, counts)
  colData <- .rsdeseq2_validate_optional_col_data(colData, counts)
  factorLevels <- factorLevels %||% .rsdeseq2_factor_levels_from_col_data(colData)
  factorLevels <- .rsdeseq2_validate_optional_factor_levels(factorLevels)
  factorReferences <- .rsdeseq2_validate_optional_factor_references(factorReferences, factorLevels)
  modelMatrix <- .rsdeseq2_validate_optional_model_matrix(modelMatrix, counts, resultsNames)

  out <- list(
    baseMean = baseMean,
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    lrtStat = lrtStat,
    lrtPvalue = lrtPvalue,
    lrtDf = lrtDf,
    resultsNames = resultsNames,
    counts = counts,
    sampleLevels = sampleLevels,
    colData = colData,
    factorLevels = factorLevels,
    factorReferences = factorReferences,
    modelMatrix = modelMatrix,
    rowNames = rowNames
  )
  class(out) <- c("rsdeseq2PrimitiveLrtFit", "list")
  out
}

resolveResultsContrastRust <- function(contrast,
                                       resultsNames,
                                       listValues = c(1, -1),
                                       reference = NULL) {
  resultsNames <- .rsdeseq2_validate_results_names(resultsNames)
  parsed <- .rsdeseq2_parse_results_contrast(
    contrast = contrast,
    resultsNames = resultsNames,
    listValues = listValues,
    reference = reference
  )
  class(parsed) <- c("rsdeseq2ResultsContrast", "list")
  parsed
}

resultsNamesRust <- function(object) {
  if (inherits(object, "DESeqDataSet")) {
    fit <- S4Vectors::metadata(object)$rsdeseq2RustFit
    if (!is.null(fit)) return(.rsdeseq2_validate_results_names(fit$resultsNames))
  }
  if (inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit"))) {
    return(.rsdeseq2_validate_results_names(object$resultsNames))
  }
  if (is.character(object)) {
    return(.rsdeseq2_validate_results_names(object))
  }
  if (is.matrix(object)) {
    beta <- .rsdeseq2_validate_beta_matrix(object, colnames(object))
    return(colnames(beta))
  }
  if (is.list(object)) {
    resultNames <- .rsdeseq2_object_field(object, "resultsNames")
    if (!is.null(resultNames)) {
      beta <- .rsdeseq2_object_field(object, "beta")
      if (!is.null(beta)) {
        beta <- .rsdeseq2_validate_beta_matrix(beta, resultNames)
        return(colnames(beta))
      }
      return(.rsdeseq2_validate_results_names(resultNames))
    }
    beta <- .rsdeseq2_object_field(object, "beta")
    if (!is.null(beta)) {
      return(resultsNamesRust(beta))
    }
  }
  fields <- .rsdeseq2_extract_wald_fit_fields(object)
  if (!is.null(fields$resultsNames)) {
    return(.rsdeseq2_validate_results_names(fields$resultsNames))
  }
  .rsdeseq2_not_available("resultsNamesRust object integration")
}

resultsNames <- function(object, ...) {
  resultsNamesRust(object)
}
