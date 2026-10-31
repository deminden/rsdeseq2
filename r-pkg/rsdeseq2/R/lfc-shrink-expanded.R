.rsdeseq2_expanded_contrast_design <- function(dds, standardDesign, contrast) {
  formula <- methods::slot(dds, "design")
  if (!inherits(formula, "formula")) {
    stop("normal contrast shrinkage requires a formula design", call. = FALSE)
  }
  if (any(attr(stats::terms(formula), "order") > 1L)) {
    stop("normal contrast shrinkage does not support formula interactions", call. = FALSE)
  }
  colData <- as.data.frame(SummarizedExperiment::colData(dds))
  modelFrame <- .rsdeseq2_model_frame(formula, data = colData)
  factorNames <- names(modelFrame)[vapply(modelFrame, is.factor, logical(1))]
  if (!all(factorNames %in% names(colData))) {
    stop("normal contrast shrinkage currently requires factors used directly in the formula", call. = FALSE)
  }
  if (is.character(contrast)) {
    if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
      stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
    }
    factorName <- contrast[[1L]]
    if (!factorName %in% factorNames) {
      stop(sprintf("contrast factor '%s' must name a factor in the formula", factorName), call. = FALSE)
    }
    selectedLevels <- base::levels(modelFrame[[factorName]])
    if (!contrast[[2L]] %in% selectedLevels || !contrast[[3L]] %in% selectedLevels) {
      stop("contrast numerator and denominator must be declared factor levels", call. = FALSE)
    }
    if (identical(contrast[[2L]], contrast[[3L]])) {
      stop("contrast numerator and denominator must differ", call. = FALSE)
    }
  } else if (is.numeric(contrast) && is.null(dim(contrast))) {
    factorName <- NULL
  } else if (is.list(contrast)) {
    if (length(contrast) == 1L) contrast <- list(contrast[[1L]], character())
    if (length(contrast) != 2L || !all(vapply(contrast, is.character, logical(1)))) {
      stop("list contrast must contain one or two character vectors", call. = FALSE)
    }
    if (anyNA(unlist(contrast, use.names = FALSE)) || any(!nzchar(unlist(contrast, use.names = FALSE)))) {
      stop("list contrast coefficient names must be non-empty", call. = FALSE)
    }
    contrastNames <- unlist(contrast, use.names = FALSE)
    if (anyDuplicated(contrastNames) || length(contrastNames) == 0L) {
      stop("list contrast must select distinct expanded-design coefficients", call. = FALSE)
    }
    factorName <- NULL
  } else {
    stop("contrast must be a factor-level triplet, list, or numeric expanded-design vector", call. = FALSE)
  }
  nullLevel <- "_null_level_"
  if (any(vapply(factorNames, function(name) nullLevel %in% levels(modelFrame[[name]]), logical(1)))) {
    stop("a formula factor already uses the reserved expanded-model level '_null_level_'", call. = FALSE)
  }
  expandedData <- rbind(colData, colData[nrow(colData), , drop = FALSE])
  if (length(factorNames)) {
    values <- .rsdeseq2_formula_request(list(operation="expanded_factors",
      values=unname(expandedData[factorNames]),null_level=nullLevel,
      collation=.rsdeseq2_formula_collation()))
    for (i in seq_along(factorNames)) expandedData[[factorNames[[i]]]] <- values[[i]]
  }
  expanded <- .rsdeseq2_model_matrix(formula, data = expandedData)
  nullColumns <- grepl(nullLevel, colnames(expanded), fixed = TRUE)
  expanded <- expanded[-nrow(expanded), !nullColumns, drop = FALSE]
  colnames(expanded)[colnames(expanded) == "(Intercept)"] <- "Intercept"
  colnames(expanded) <- make.names(colnames(expanded))
  rowOrder <- match(rownames(SummarizedExperiment::colData(dds)), rownames(expanded))
  if (anyNA(rowOrder)) stop("expanded design rows do not match dataset samples", call. = FALSE)
  expanded <- expanded[rowOrder, , drop = FALSE]

  expandedGroups <- integer(ncol(expanded))
  standardGroups <- integer(ncol(standardDesign))
  expandedPriorMap <- integer(ncol(expanded))
  standardNames <- make.names(ifelse(colnames(standardDesign) == "(Intercept)", "Intercept", colnames(standardDesign)))
  expandedFactorColumns <- vector("list", length(factorNames))
  names(expandedFactorColumns) <- factorNames
  for (group in seq_along(factorNames)) {
    name <- factorNames[[group]]
    factorLevels <- levels(modelFrame[[name]])
    expandedNames <- make.names(paste0(name, factorLevels))
    expandedColumns <- match(expandedNames, colnames(expanded))
    if (anyNA(expandedColumns)) {
      stop(sprintf("could not map expanded coefficients for factor '%s'", name), call. = FALSE)
    }
    expandedGroups[expandedColumns] <- group
    expandedFactorColumns[[name]] <- expandedColumns
    standardFactorLevels <- if (attr(stats::terms(formula), "intercept") == 1L) {
      factorLevels[-1L]
    } else {
      factorLevels
    }
    standardNamesForLevels <- make.names(paste0(name, standardFactorLevels))
    standardColumns <- match(standardNamesForLevels, standardNames)
    standardColumns <- standardColumns[!is.na(standardColumns)]
    if (length(standardColumns)) standardGroups[standardColumns] <- group
  }
  for (column in seq_len(ncol(expanded))) {
    if (expandedGroups[[column]] > 0L) next
    standardColumn <- match(colnames(expanded)[[column]], standardNames)
    if (is.na(standardColumn)) {
      stop(sprintf("could not map expanded coefficient '%s' to the fitted design", colnames(expanded)[[column]]), call. = FALSE)
    }
    expandedPriorMap[[column]] <- standardColumn
  }
  if (is.numeric(contrast)) {
    if (length(contrast) != ncol(expanded)) {
      stop("numeric contrast vector must have one value per expanded-design coefficient", call. = FALSE)
    }
    if (any(!is.finite(contrast)) || all(contrast == 0)) {
      stop("numeric contrast vector must be finite and contain a non-zero value", call. = FALSE)
    }
    contrastVector <- as.numeric(contrast)
  } else if (is.list(contrast)) {
    positive <- match(contrast[[1L]], colnames(expanded))
    negative <- match(contrast[[2L]], colnames(expanded))
    if (anyNA(positive) || anyNA(negative)) {
      stop("all list-contrast coefficients must name expanded-design columns", call. = FALSE)
    }
    if (length(intersect(positive, negative))) {
      stop("list contrast coefficients must not appear in both numerator and denominator", call. = FALSE)
    }
    contrastVector <- .Call("rsdeseq2_average_contrast", as.integer(ncol(expanded)),
      as.integer(positive), as.integer(negative), PACKAGE="rsdeseq2")
  } else {
    contrastVector <- numeric(ncol(expanded))
    selectedColumns <- expandedFactorColumns[[factorName]]
    contrastVector[[selectedColumns[match(contrast[[2L]], selectedLevels)]]] <- 1
    contrastVector[[selectedColumns[match(contrast[[3L]], selectedLevels)]]] <- -1
  }
  list(
    design = expanded,
    expandedGroups = as.integer(expandedGroups),
    standardGroups = as.integer(standardGroups),
    expandedPriorMap = as.integer(expandedPriorMap),
    contrast = contrastVector,
    interceptIndex = match("Intercept", colnames(expanded)),
    standardInterceptIndex = match("Intercept", standardNames)
  )
}
