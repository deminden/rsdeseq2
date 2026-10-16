# Result-table assembly, independent filtering and Cook's filtering.

resultsTableRust <- function(baseMean,
                             log2FoldChange = NULL,
                             lfcSE = NULL,
                             stat = NULL,
                             pvalue = NULL,
                             padj = NULL,
                             dispersion = NULL,
                             converged = NULL,
                             contrast = NULL,
                             rowNames = NULL) {
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  nGenes <- length(baseMean)
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% names(baseMean), nGenes)
  baseMean <- .rsdeseq2_align_named_gene_vector(baseMean, rowNames, "baseMean")

  if (is.null(pvalue)) {
    if (!is.null(padj)) {
      stop("padj cannot be supplied when pvalue is NULL", call. = FALSE)
    }
    padj <- NULL
  } else {
    pvalue <- .rsdeseq2_validate_result_pvalues(pvalue, nGenes)
    pvalue <- .rsdeseq2_align_named_gene_vector(pvalue, rowNames, "pvalue")
    if (is.null(padj)) {
      padj <- stats::p.adjust(pvalue, method = "BH")
    } else {
      padj <- .rsdeseq2_validate_result_pvalues(padj, nGenes, name = "padj")
      padj <- .rsdeseq2_align_named_gene_vector(padj, rowNames, "padj")
    }
  }

  log2FoldChange <- .rsdeseq2_validate_result_optional_numeric(log2FoldChange, nGenes, "log2FoldChange")
  log2FoldChange <- .rsdeseq2_align_named_gene_vector(log2FoldChange, rowNames, "log2FoldChange")
  lfcSE <- .rsdeseq2_validate_result_optional_numeric(lfcSE, nGenes, "lfcSE")
  lfcSE <- .rsdeseq2_align_named_gene_vector(lfcSE, rowNames, "lfcSE")
  stat <- .rsdeseq2_validate_result_optional_numeric(stat, nGenes, "stat")
  stat <- .rsdeseq2_align_named_gene_vector(stat, rowNames, "stat")

  out <- data.frame(
    baseMean = baseMean,
    log2FoldChange = log2FoldChange,
    lfcSE = lfcSE,
    stat = stat,
    pvalue = pvalue %||% rep(NA_real_, nGenes),
    padj = padj %||% rep(NA_real_, nGenes),
    check.names = FALSE
  )

  if (!is.null(dispersion)) {
    dispersion <- .rsdeseq2_validate_result_optional_numeric(dispersion, nGenes, "dispersion")
    out$dispersion <- .rsdeseq2_align_named_gene_vector(dispersion, rowNames, "dispersion")
  }
  if (!is.null(converged)) {
    converged <- .rsdeseq2_validate_result_optional_logical(converged, nGenes, "converged")
    out$converged <- .rsdeseq2_align_named_gene_vector(converged, rowNames, "converged")
  }
  if (!is.null(contrast)) {
    contrast <- .rsdeseq2_validate_result_table_contrast(contrast)
    attr(out, "contrast") <- contrast
  }
  if (!is.null(rowNames)) {
    rownames(out) <- rowNames
  }
  out
}

.rsdeseq2_validate_result_table_contrast <- function(contrast) {
  if (!is.numeric(contrast) || is.list(contrast)) {
    stop("contrast metadata must be a numeric vector", call. = FALSE)
  }
  if (length(contrast) == 0L) {
    stop("contrast metadata cannot be empty", call. = FALSE)
  }
  if (any(!is.finite(contrast))) {
    stop("contrast metadata must contain finite values", call. = FALSE)
  }
  out <- as.numeric(contrast)
  names(out) <- names(contrast)
  out
}

applyIndependentFilteringRust <- function(results,
                                          alpha = 0.1,
                                          theta = NULL,
                                          enabled = TRUE) {
  results <- .rsdeseq2_validate_result_table(results)
  alpha <- .rsdeseq2_validate_independent_filter_alpha(alpha)

  if (!isTRUE(enabled)) {
    results$padj <- stats::p.adjust(results$pvalue, method = "BH")
    results$filtered <- rep(NA, nrow(results))
    attr(results, "independentFiltering") <- list(
      enabled = FALSE,
      theta = numeric(),
      numRej = integer(),
      selected = NA_integer_,
      filterTheta = NA_real_,
      filterThreshold = NA_real_,
      lo.fit = NULL,
      alpha = alpha
    )
    return(results)
  }

  filter <- results$baseMean
  theta <- .rsdeseq2_validate_independent_filter_theta(theta, filter)
  cutoffs <- as.numeric(stats::quantile(filter, probs = theta, names = FALSE, type = 7))
  filteredPadj <- vector("list", length(cutoffs))
  numRej <- integer(length(cutoffs))

  for (idx in seq_along(cutoffs)) {
    selected <- filter >= cutoffs[[idx]]
    adjusted <- rep(NA_real_, length(filter))
    adjusted[selected] <- stats::p.adjust(results$pvalue[selected], method = "BH")
    filteredPadj[[idx]] <- adjusted
    numRej[[idx]] <- sum(adjusted < alpha, na.rm = TRUE)
  }

  loFit <- stats::lowess(theta, numRej, f = 1 / 5)
  selectedIndex <- .rsdeseq2_select_independent_filter_index(numRej, loFit$y)
  filterThreshold <- cutoffs[[selectedIndex]]

  results$padj <- filteredPadj[[selectedIndex]]
  results$filtered <- ifelse(is.na(results$pvalue), NA, filter < filterThreshold)
  attr(results, "independentFiltering") <- list(
    enabled = TRUE,
    theta = theta,
    numRej = numRej,
    selected = selectedIndex,
    filterTheta = theta[[selectedIndex]],
    filterThreshold = filterThreshold,
    lo.fit = loFit,
    alpha = alpha
  )
  results
}

applyCooksCutoffRust <- function(pvalue,
                                 maxCooks,
                                 cooksCutoff,
                                 counts = NULL,
                                 cooks = NULL,
                                 lowCountHeuristic = FALSE,
                                 native = TRUE) {
  pvalue <- .rsdeseq2_validate_pvalues(pvalue)
  maxCooks <- .rsdeseq2_validate_optional_numeric_vector(maxCooks, length(pvalue), "maxCooks")
  maxCooks <- .rsdeseq2_align_named_gene_vector(maxCooks, names(pvalue), "maxCooks")
  cooksCutoff <- .rsdeseq2_validate_cooks_cutoff(cooksCutoff)

  if (is.null(cooksCutoff)) {
    return(.rsdeseq2_cooks_cutoff_output(pvalue, maxCooks, rep(NA, length(pvalue))))
  }

  if (isTRUE(lowCountHeuristic)) {
    counts <- .rsdeseq2_validate_counts(counts)
    if (nrow(counts) != length(pvalue)) {
      stop("counts must have one row per p-value", call. = FALSE)
    }
    counts <- .rsdeseq2_align_named_gene_matrix(counts, names(pvalue), "counts")
    cooks <- .rsdeseq2_validate_cooks_matrix(cooks, nrow(counts), ncol(counts))
    cooks <- .rsdeseq2_align_cooks_matrix(cooks, counts)
  } else {
    counts <- NULL
    cooks <- NULL
  }

  nativeMasked <- .Call(
    "rsdeseq2_apply_cooks_cutoff",
    pvalue,
    maxCooks,
    cooksCutoff,
    counts,
    cooks,
    isTRUE(lowCountHeuristic),
    PACKAGE = "rsdeseq2"
  )
  .rsdeseq2_cooks_cutoff_output(
    nativeMasked$pvalue,
    maxCooks,
    nativeMasked$cooksOutlier
  )
}

.rsdeseq2_align_cooks_matrix <- function(cooks, counts) {
  cooksNames <- rownames(cooks)
  countNames <- rownames(counts)
  if (!is.null(countNames) &&
    !is.null(cooksNames) &&
    .rsdeseq2_has_explicit_row_names(cooksNames, nrow(cooks))) {
    if (anyNA(cooksNames) || any(!nzchar(cooksNames)) || anyDuplicated(cooksNames)) {
      stop("cooks row names must be unique non-empty gene names", call. = FALSE)
    }
    missing <- setdiff(countNames, cooksNames)
    if (length(missing) > 0L) {
      stop("cooks row names must contain all count row names", call. = FALSE)
    }
    cooks <- cooks[countNames, , drop = FALSE]
  }

  cooksSamples <- colnames(cooks)
  countSamples <- colnames(counts)
  if (!is.null(countSamples) &&
    !is.null(cooksSamples) &&
    .rsdeseq2_has_explicit_row_names(cooksSamples, ncol(cooks))) {
    if (anyNA(cooksSamples) || any(!nzchar(cooksSamples)) || anyDuplicated(cooksSamples)) {
      stop("cooks column names must be unique non-empty sample names", call. = FALSE)
    }
    missing <- setdiff(countSamples, cooksSamples)
    if (length(missing) > 0L) {
      stop("cooks column names must contain all count sample names", call. = FALSE)
    }
    cooks <- cooks[, countSamples, drop = FALSE]
  }
  cooks
}
