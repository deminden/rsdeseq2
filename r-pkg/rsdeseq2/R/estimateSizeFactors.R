estimateSizeFactorsRust <- function(counts,
                                    method = c("ratio", "poscounts"),
                                    geoMeans = NULL,
                                    controlGenes = NULL,
                                    native = FALSE) {
  method <- match.arg(method)
  counts <- .rsdeseq2_validate_counts(counts)
  controlIdx <- .rsdeseq2_control_gene_indices(controlGenes, counts)
  incomingGeoMeans <- !is.null(geoMeans)

  logGeoMeans <- .rsdeseq2_validate_geo_means(geoMeans, nrow(counts))
  if (is.null(logGeoMeans)) {
    if (method == "ratio") {
      logGeoMeans <- apply(counts, 1L, function(row) {
        if (any(row == 0)) {
          -Inf
        } else {
          mean(log(row))
        }
      })
    } else {
      nSamples <- ncol(counts)
      logGeoMeans <- apply(counts, 1L, function(row) {
        if (all(row == 0)) {
          -Inf
        } else {
          sum(log(row[row > 0])) / nSamples
        }
      })
    }
  }

  if (all(is.infinite(logGeoMeans))) {
    stop("no usable genes for size-factor estimation", call. = FALSE)
  }

  if (isTRUE(native)) {
    sizeFactors <- .rsdeseq2_try_native_size_factors(
      counts = counts,
      logGeoMeans = logGeoMeans,
      controlIdx = controlIdx,
      stabilize = incomingGeoMeans
    )
    if (!is.null(sizeFactors)) {
      return(sizeFactors)
    }
  }
  .rsdeseq2_size_factors_fallback(
    counts = counts,
    logGeoMeans = logGeoMeans,
    controlIdx = controlIdx,
    stabilize = incomingGeoMeans
  )
}

.rsdeseq2_size_factors_fallback <- function(counts,
                                            logGeoMeans,
                                            controlIdx,
                                            stabilize) {
  sizeFactors <- numeric(ncol(counts))
  for (sample in seq_len(ncol(counts))) {
    usable <- is.finite(logGeoMeans[controlIdx]) & counts[controlIdx, sample] > 0
    if (!any(usable)) {
      stop(sprintf("sample %d has no usable positive count ratios", sample), call. = FALSE)
    }
    logRatios <- log(counts[controlIdx[usable], sample]) - logGeoMeans[controlIdx[usable]]
    sizeFactors[sample] <- exp(stats::median(logRatios))
  }
  if (any(!is.finite(sizeFactors)) || any(sizeFactors <= 0)) {
    stop("estimated size factors must be finite and positive", call. = FALSE)
  }
  if (stabilize) {
    sizeFactors <- .rsdeseq2_stabilize_size_factors(sizeFactors)
  }
  names(sizeFactors) <- colnames(counts)
  sizeFactors
}

.rsdeseq2_try_native_size_factors <- function(counts,
                                              logGeoMeans,
                                              controlIdx,
                                              stabilize) {
  args <- list(counts, logGeoMeans, as.integer(controlIdx), as.logical(stabilize))
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_estimate_size_factors", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_estimate_size_factors"), args)),
    error = function(error) NULL
  )
}
