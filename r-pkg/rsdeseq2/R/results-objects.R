# Adapters for fitted objects, metadata columns and primitive fit inputs.

.rsdeseq2_is_list_like_lrt_input <- function(object) {
  !inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit")) &&
    {
      fields <- .rsdeseq2_extract_lrt_fit_fields(object)
      !is.null(fields$baseMean) &&
        !is.null(fields$beta) &&
        !is.null(fields$lrtStat)
    }
}

.rsdeseq2_is_list_like_wald_input <- function(object) {
  !inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit")) &&
    {
      fields <- .rsdeseq2_extract_wald_fit_fields(object)
      !is.null(fields$baseMean) &&
        !is.null(fields$beta)
    }
}

.rsdeseq2_extract_wald_fit_fields <- function(object) {
  mcols <- .rsdeseq2_coerce_mcols_frame(.rsdeseq2_object_field(object, "mcols"))
  assays <- .rsdeseq2_object_field(object, "assays")
  colData <- .rsdeseq2_object_field(object, "colData")
  explicitResultsNames <- .rsdeseq2_object_field(object, "resultsNames")
  beta <- .rsdeseq2_object_field(object, "beta") %||%
    .rsdeseq2_object_field(object, "betaMatrix") %||%
    .rsdeseq2_mcols_beta_matrix(mcols, explicitResultsNames)
  resultsNames <- explicitResultsNames %||% colnames(beta)
  list(
    baseMean = .rsdeseq2_object_field(object, "baseMean") %||%
      .rsdeseq2_mcols_column(mcols, "baseMean"),
    beta = beta,
    betaCovariance = .rsdeseq2_object_field(object, "betaCovariance"),
    betaSE = .rsdeseq2_object_field(object, "betaSE") %||%
      .rsdeseq2_mcols_prefixed_matrix(mcols, "SE_", resultsNames),
    resultsNames = resultsNames,
    counts = .rsdeseq2_object_field(object, "counts") %||%
      .rsdeseq2_object_field(assays, "counts"),
    sampleLevels = .rsdeseq2_object_field(object, "sampleLevels"),
    colData = colData,
    factorLevels = .rsdeseq2_object_field(object, "factorLevels") %||%
      .rsdeseq2_factor_levels_from_col_data(colData),
    factorReferences = .rsdeseq2_object_field(object, "factorReferences"),
    modelMatrix = .rsdeseq2_object_field(object, "modelMatrix") %||%
      .rsdeseq2_object_field(object, "design"),
    rowNames = .rsdeseq2_object_field(object, "rowNames") %||%
      rownames(beta) %||%
      rownames(mcols)
  )
}

.rsdeseq2_extract_lrt_fit_fields <- function(object) {
  fields <- .rsdeseq2_extract_wald_fit_fields(object)
  mcols <- .rsdeseq2_coerce_mcols_frame(.rsdeseq2_object_field(object, "mcols"))
  fields$lrtStat <- .rsdeseq2_object_field(object, "lrtStat") %||%
    .rsdeseq2_object_field(object, "stat") %||%
    .rsdeseq2_mcols_column(mcols, "stat")
  fields$lrtPvalue <- .rsdeseq2_object_field(object, "lrtPvalue") %||%
    .rsdeseq2_object_field(object, "pvalue") %||%
    .rsdeseq2_mcols_column(mcols, "pvalue")
  fields$lrtDf <- .rsdeseq2_object_field(object, "lrtDf")
  fields
}

.rsdeseq2_object_field <- function(object, name) {
  if (is.null(object)) {
    return(NULL)
  }
  if (is.list(object) && !is.null(object[[name]])) {
    return(object[[name]])
  }
  for (alias in .rsdeseq2_object_field_aliases(name)) {
    if (is.list(object) && !is.null(object[[alias]])) {
      return(object[[alias]])
    }
  }
  accessed <- .rsdeseq2_optional_object_accessor(object, name)
  if (!is.null(accessed)) {
    return(accessed)
  }
  attributed <- .rsdeseq2_object_attribute_field(object, name)
  if (!is.null(attributed)) {
    return(attributed)
  }
  if (isS4(object)) {
    slots <- methods::slotNames(object)
    if (name %in% slots) {
      return(methods::slot(object, name))
    }
    for (alias in .rsdeseq2_object_field_aliases(name)) {
      if (alias %in% slots) {
        return(methods::slot(object, alias))
      }
    }
  }
  bracketAccessed <- .rsdeseq2_optional_bracket_field(object, name)
  if (!is.null(bracketAccessed)) {
    return(bracketAccessed)
  }
  NULL
}

.rsdeseq2_object_attribute_field <- function(object, name) {
  if (is.null(object) || is.null(name)) {
    return(NULL)
  }
  value <- attr(object, name, exact = TRUE)
  if (!is.null(value)) {
    return(value)
  }
  for (alias in .rsdeseq2_object_field_aliases(name)) {
    value <- attr(object, alias, exact = TRUE)
    if (!is.null(value)) {
      return(value)
    }
  }
  NULL
}

.rsdeseq2_object_field_aliases <- function(name) {
  switch(
    name,
    beta = c("coef", "coefficients"),
    betaMatrix = c("coefMatrix", "coefficientMatrix"),
    betaCovariance = c("betaCov", "coefCovariance", "coefficientCovariance"),
    factorLevels = c("factorLevelNames", "factorLevelsList", "levelNames"),
    factorReferences = c("factorReference", "factorReferenceLevels", "referenceLevels"),
    mcols = c("elementMetadata", "rowData"),
    modelMatrix = c("model.matrix", "designMatrix", "fullModelMatrix"),
    resultsNames = c("resultNames", "coefNames", "coefficientNames"),
    rowNames = c("NAMES"),
    character()
  )
}

.rsdeseq2_optional_bracket_field <- function(object, name) {
  if (is.null(object) || is.null(name)) {
    return(NULL)
  }
  tryCatch(
    object[[name]],
    error = function(e) NULL
  )
}

.rsdeseq2_optional_object_accessor <- function(object, name) {
  accessors <- switch(
    name,
    mcols = list(c("BiocGenerics", "mcols"), c("S4Vectors", "mcols")),
    assays = list(c("SummarizedExperiment", "assays")),
    colData = list(c("SummarizedExperiment", "colData")),
    counts = list(c("DESeq2", "counts")),
    sizeFactors = list(c("DESeq2", "sizeFactors")),
    rowNames = list(c("BiocGenerics", "rownames")),
    list()
  )
  for (accessor in accessors) {
    namespace <- accessor[[1L]]
    symbol <- accessor[[2L]]
    if (requireNamespace(namespace, quietly = TRUE) && exists(symbol, envir = asNamespace(namespace), mode = "function")) {
      value <- tryCatch(
        get(symbol, envir = asNamespace(namespace), mode = "function")(object),
        error = function(e) NULL
      )
      if (!is.null(value)) {
        return(value)
      }
    }
  }
  NULL
}

.rsdeseq2_mcols_column <- function(mcols, name) {
  mcols <- .rsdeseq2_coerce_mcols_frame(mcols)
  column <- .rsdeseq2_mcols_resolve_column_name(name, mcols)
  if (is.null(column) || is.na(column)) {
    return(NULL)
  }
  out <- mcols[[column]]
  .rsdeseq2_apply_mcols_row_names(out, rownames(mcols))
}

.rsdeseq2_apply_mcols_row_names <- function(value, rowNames) {
  if (is.null(value) || is.null(rowNames)) {
    return(value)
  }
  if (is.vector(value) && is.null(names(value)) &&
    .rsdeseq2_has_explicit_row_names(rowNames, length(value))) {
    names(value) <- rowNames
    return(value)
  }
  if (is.matrix(value) && is.null(rownames(value)) &&
    .rsdeseq2_has_explicit_row_names(rowNames, nrow(value))) {
    rownames(value) <- rowNames
    return(value)
  }
  if (is.array(value) && length(dim(value)) >= 1L &&
    .rsdeseq2_has_explicit_row_names(rowNames, dim(value)[[1L]])) {
    names <- dimnames(value)
    if (is.null(names)) {
      names <- vector("list", length(dim(value)))
    }
    if (is.null(names[[1L]])) {
      names[[1L]] <- rowNames
      dimnames(value) <- names
    }
  }
  value
}

.rsdeseq2_coerce_mcols_frame <- function(mcols) {
  if (is.null(mcols)) {
    return(NULL)
  }
  if (is.data.frame(mcols)) {
    return(mcols)
  }
  coerced <- tryCatch(
    as.data.frame(mcols),
    error = function(e) NULL
  )
  if (is.null(coerced) || !is.data.frame(coerced)) {
    return(mcols)
  }
  coerced
}

.rsdeseq2_factor_levels_from_col_data <- function(colData) {
  if (is.null(colData)) {
    return(NULL)
  }
  colData <- .rsdeseq2_coerce_col_data_frame(colData)
  if (is.null(colData)) {
    return(NULL)
  }
  names <- colnames(colData) %||% names(colData) %||% character()
  out <- list()
  for (name in names) {
    column <- colData[[name]]
    if (is.factor(column)) {
      out[[name]] <- levels(column)
    } else if (is.character(column) && length(column) > 0L && !anyNA(column) && all(nzchar(column))) {
      out[[name]] <- sort(unique(column))
    }
  }
  if (length(out) == 0L) {
    return(NULL)
  }
  out
}

.rsdeseq2_mcols_beta_matrix <- function(mcols, resultsNames = NULL) {
  mcols <- .rsdeseq2_coerce_mcols_frame(mcols)
  if (is.null(mcols)) {
    return(NULL)
  }
  betaMatrix <- .rsdeseq2_mcols_column(mcols, "betaMatrix")
  if (!is.null(betaMatrix)) {
    return(betaMatrix)
  }
  if (is.null(colnames(mcols))) {
    return(NULL)
  }
  if (!is.null(resultsNames)) {
    betaColumns <- .rsdeseq2_mcols_resolve_column_names(mcols, resultsNames)
    if (is.null(betaColumns)) {
      return(NULL)
    }
    out <- as.matrix(mcols[, betaColumns, drop = FALSE])
    colnames(out) <- resultsNames
    return(out)
  }
  seColumns <- grep("^SE_", colnames(mcols), value = TRUE)
  betaColumns <- sub("^SE_", "", seColumns)
  betaColumns <- betaColumns[betaColumns %in% colnames(mcols)]
  if (length(betaColumns) == 0L) {
    return(NULL)
  }
  as.matrix(mcols[, betaColumns, drop = FALSE])
}

.rsdeseq2_mcols_prefixed_matrix <- function(mcols, prefix, resultsNames) {
  mcols <- .rsdeseq2_coerce_mcols_frame(mcols)
  if (is.null(mcols) || is.null(resultsNames)) {
    return(NULL)
  }
  columns <- .rsdeseq2_mcols_resolve_column_names(mcols, paste0(prefix, resultsNames))
  if (is.null(columns)) {
    return(NULL)
  }
  out <- as.matrix(mcols[, columns, drop = FALSE])
  colnames(out) <- resultsNames
  out
}

.rsdeseq2_mcols_resolve_column_names <- function(mcols, names) {
  mcols <- .rsdeseq2_coerce_mcols_frame(mcols)
  columns <- vapply(
    names,
    .rsdeseq2_mcols_resolve_column_name,
    character(1L),
    mcols = mcols,
    USE.NAMES = FALSE
  )
  if (anyNA(columns)) {
    return(NULL)
  }
  columns
}

.rsdeseq2_mcols_resolve_column_name <- function(name, mcols) {
  mcols <- .rsdeseq2_coerce_mcols_frame(mcols)
  columns <- colnames(mcols)
  if (is.null(columns) || is.null(name)) {
    return(NA_character_)
  }
  exact <- which(columns == name)
  if (length(exact) == 1L) {
    return(columns[[exact]])
  }
  if (length(exact) > 1L) {
    stop(sprintf("mcols column %s is duplicated", name), call. = FALSE)
  }
  alias <- make.names(name)
  aliases <- make.names(columns)
  matches <- which(aliases == alias)
  if (length(matches) == 1L) {
    return(columns[[matches]])
  }
  if (length(matches) > 1L) {
    stop(sprintf("mcols column alias %s is ambiguous", name), call. = FALSE)
  }
  NA_character_
}

.rsdeseq2_results_primitive_wald_fit <- function(object,
                                                 contrast,
                                                 name,
                                                 listValues,
                                                 reference,
                                                 lfcThreshold,
                                                 altHypothesis,
                                                 useT,
                                                 degreesOfFreedom,
                                                 alpha,
                                                 independentFiltering) {
  resultsNames <- object$resultsNames
  reference <- reference %||% .rsdeseq2_reference_for_character_contrast(object, contrast)
  contrast <- .rsdeseq2_canonicalize_character_contrast(object, contrast)
  reference <- .rsdeseq2_canonicalize_character_contrast_reference(object, contrast, reference)
  if (is.null(contrast)) {
    if (is.null(name)) {
      if (length(resultsNames) < 2L) {
        stop("default results require at least two result names", call. = FALSE)
      }
      name <- resultsNames[[length(resultsNames)]]
    }
    coefficient <- .rsdeseq2_resolve_results_name_index(name, resultsNames)
    numericContrast <- numeric(length(resultsNames))
    numericContrast[[coefficient]] <- 1
    names(numericContrast) <- resultsNames
    resultName <- resultsNames[[coefficient]]
    comparison <- sprintf("coefficient %s", resultName)
    allZero <- list(type = "none")
  } else {
    resolved <- resolveResultsContrastRust(
      contrast = contrast,
      resultsNames = resultsNames,
      listValues = listValues,
      reference = reference
    )
    numericContrast <- resolved$numeric
    resultName <- resolved$resultName
    comparison <- resolved$comparison
    allZero <- resolved$allZero
  }

  waldOptions <- .rsdeseq2_validate_wald_result_options(
    lfcThreshold = lfcThreshold,
    altHypothesis = altHypothesis,
    useT = useT,
    degreesOfFreedom = degreesOfFreedom,
    nGenes = nrow(object$beta)
  )
  wald <- .rsdeseq2_wald_from_primitive_fit(object, numericContrast, waldOptions)
  wald <- .rsdeseq2_apply_results_contrast_all_zero(
    wald = wald,
    fit = object,
    numericContrast = numericContrast,
    allZero = allZero
  )
  out <- resultsTableRust(
    baseMean = object$baseMean,
    log2FoldChange = wald$log2FoldChange,
    lfcSE = wald$lfcSE,
    stat = wald$stat,
    pvalue = wald$pvalue,
    contrast = numericContrast,
    rowNames = object$rowNames
  )
  attr(out, "resultName") <- resultName
  attr(out, "comparison") <- comparison
  attr(out, "contrastAllZero") <- allZero$type
  attr(out, "lfcThreshold") <- waldOptions$lfcThreshold
  attr(out, "altHypothesis") <- waldOptions$altHypothesis
  if (!is.null(waldOptions$degreesOfFreedom)) {
    attr(out, "degreesOfFreedom") <- waldOptions$degreesOfFreedom
  }
  applyIndependentFilteringRust(
    out,
    alpha = alpha,
    enabled = independentFiltering
  )
}

.rsdeseq2_results_primitive_lrt_fit <- function(object,
                                                contrast,
                                                name,
                                                listValues,
                                                reference,
                                                lfcThreshold,
                                                altHypothesis,
                                                useT,
                                                degreesOfFreedom,
                                                alpha,
                                                independentFiltering) {
  if (!is.numeric(lfcThreshold) || length(lfcThreshold) != 1L || !is.finite(lfcThreshold) || lfcThreshold != 0 ||
    !identical(altHypothesis, "greaterAbs") || isTRUE(useT) || !is.null(degreesOfFreedom)) {
    stop("lfcThreshold, altHypothesis, useT, and degreesOfFreedom are only supported for Wald results", call. = FALSE)
  }
  resultsNames <- object$resultsNames
  reference <- reference %||% .rsdeseq2_reference_for_character_contrast(object, contrast)
  contrast <- .rsdeseq2_canonicalize_character_contrast(object, contrast)
  reference <- .rsdeseq2_canonicalize_character_contrast_reference(object, contrast, reference)
  if (is.null(contrast)) {
    if (is.null(name)) {
      if (length(resultsNames) < 2L) {
        stop("default results require at least two result names", call. = FALSE)
      }
      name <- resultsNames[[length(resultsNames)]]
    }
    coefficient <- .rsdeseq2_resolve_results_name_index(name, resultsNames)
    numericContrast <- numeric(length(resultsNames))
    numericContrast[[coefficient]] <- 1
    names(numericContrast) <- resultsNames
    resultName <- resultsNames[[coefficient]]
    comparison <- sprintf("coefficient %s; LRT full vs reduced", resultName)
    allZero <- list(type = "none")
  } else {
    resolved <- resolveResultsContrastRust(
      contrast = contrast,
      resultsNames = resultsNames,
      listValues = listValues,
      reference = reference
    )
    numericContrast <- resolved$numeric
    resultName <- resolved$resultName
    comparison <- sprintf("%s; LRT full vs reduced", resolved$comparison)
    allZero <- resolved$allZero
  }

  lrt <- .rsdeseq2_lrt_from_primitive_fit(object, numericContrast)
  lrt <- .rsdeseq2_apply_lrt_results_contrast_all_zero(
    lrt = lrt,
    fit = object,
    numericContrast = numericContrast,
    allZero = allZero
  )
  out <- resultsTableRust(
    baseMean = object$baseMean,
    log2FoldChange = lrt$log2FoldChange,
    lfcSE = lrt$lfcSE,
    stat = object$lrtStat,
    pvalue = object$lrtPvalue,
    contrast = numericContrast,
    rowNames = object$rowNames
  )
  attr(out, "resultName") <- resultName
  attr(out, "comparison") <- comparison
  attr(out, "contrastAllZero") <- allZero$type
  if (!is.null(object$lrtDf)) {
    attr(out, "lrtDf") <- object$lrtDf
  }
  applyIndependentFilteringRust(
    out,
    alpha = alpha,
    enabled = independentFiltering
  )
}
