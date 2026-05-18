results <- function(object, ...) {
  .rsdeseq2_not_available("results")
}

resultsTableRust <- function(baseMean,
                             log2FoldChange = NULL,
                             lfcSE = NULL,
                             stat = NULL,
                             pvalue = NULL,
                             padj = NULL,
                             dispersion = NULL,
                             converged = NULL,
                             rowNames = NULL) {
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  nGenes <- length(baseMean)
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% names(baseMean), nGenes)

  if (is.null(pvalue)) {
    if (!is.null(padj)) {
      stop("padj cannot be supplied when pvalue is NULL", call. = FALSE)
    }
    padj <- NULL
  } else {
    pvalue <- .rsdeseq2_validate_result_pvalues(pvalue, nGenes)
    if (is.null(padj)) {
      padj <- stats::p.adjust(pvalue, method = "BH")
    } else {
      padj <- .rsdeseq2_validate_result_pvalues(padj, nGenes, name = "padj")
    }
  }

  out <- data.frame(
    baseMean = baseMean,
    log2FoldChange = .rsdeseq2_validate_result_optional_numeric(log2FoldChange, nGenes, "log2FoldChange"),
    lfcSE = .rsdeseq2_validate_result_optional_numeric(lfcSE, nGenes, "lfcSE"),
    stat = .rsdeseq2_validate_result_optional_numeric(stat, nGenes, "stat"),
    pvalue = pvalue %||% rep(NA_real_, nGenes),
    padj = padj %||% rep(NA_real_, nGenes),
    check.names = FALSE
  )

  if (!is.null(dispersion)) {
    out$dispersion <- .rsdeseq2_validate_result_optional_numeric(dispersion, nGenes, "dispersion")
  }
  if (!is.null(converged)) {
    out$converged <- .rsdeseq2_validate_result_optional_logical(converged, nGenes, "converged")
  }
  if (!is.null(rowNames)) {
    rownames(out) <- rowNames
  }
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
                                 native = FALSE) {
  pvalue <- .rsdeseq2_validate_pvalues(pvalue)
  maxCooks <- .rsdeseq2_validate_optional_numeric_vector(maxCooks, length(pvalue), "maxCooks")
  cooksCutoff <- .rsdeseq2_validate_cooks_cutoff(cooksCutoff)

  if (is.null(cooksCutoff)) {
    return(.rsdeseq2_cooks_cutoff_output(pvalue, maxCooks, rep(NA, length(pvalue))))
  }

  if (isTRUE(lowCountHeuristic)) {
    counts <- .rsdeseq2_validate_counts(counts)
    cooks <- .rsdeseq2_validate_cooks_matrix(cooks, nrow(counts), ncol(counts))
    if (nrow(counts) != length(pvalue)) {
      stop("counts must have one row per p-value", call. = FALSE)
    }
  } else {
    counts <- NULL
    cooks <- NULL
  }

  if (isTRUE(native)) {
    nativeMasked <- .rsdeseq2_try_native_cooks_cutoff(
      pvalue = pvalue,
      maxCooks = maxCooks,
      cooksCutoff = cooksCutoff,
      counts = counts,
      cooks = cooks,
      lowCountHeuristic = lowCountHeuristic
    )
    if (!is.null(nativeMasked)) {
      return(.rsdeseq2_cooks_cutoff_output(
        nativeMasked$pvalue,
        maxCooks,
        nativeMasked$cooksOutlier
      ))
    }
  }

  .rsdeseq2_apply_cooks_cutoff_fallback(
    pvalue = pvalue,
    maxCooks = maxCooks,
    cooksCutoff = cooksCutoff,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = lowCountHeuristic
  )
}

.rsdeseq2_apply_cooks_cutoff_fallback <- function(pvalue,
                                                  maxCooks,
                                                  cooksCutoff,
                                                  counts,
                                                  cooks,
                                                  lowCountHeuristic) {
  cooksOutlier <- ifelse(is.na(maxCooks), NA, maxCooks > cooksCutoff)
  if (isTRUE(lowCountHeuristic)) {
    for (gene in which(cooksOutlier %in% TRUE)) {
      if (.rsdeseq2_low_count_cooks_spares_row(counts[gene, ], cooks[gene, ])) {
        cooksOutlier[gene] <- FALSE
      }
    }
  }

  masked <- pvalue
  masked[which(cooksOutlier %in% TRUE)] <- NA_real_
  .rsdeseq2_cooks_cutoff_output(masked, maxCooks, cooksOutlier)
}

.rsdeseq2_try_native_cooks_cutoff <- function(pvalue,
                                              maxCooks,
                                              cooksCutoff,
                                              counts,
                                              cooks,
                                              lowCountHeuristic) {
  args <- list(
    pvalue,
    maxCooks,
    cooksCutoff,
    counts %||% NULL,
    cooks %||% NULL,
    as.logical(lowCountHeuristic)
  )
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_apply_cooks_cutoff", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_apply_cooks_cutoff"), args)),
    error = function(error) NULL
  )
}
