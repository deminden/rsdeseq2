normalizedCountsRust <- function(counts,
                                 sizeFactors = NULL,
                                 normalizationFactors = NULL,
                                 native = FALSE) {
  counts <- .rsdeseq2_validate_counts(counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (!is.null(normalizationFactors)) {
    sizeFactors <- NULL
  } else {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  }
  if (isTRUE(native)) {
    normalized <- .rsdeseq2_try_native_normalized_counts(
      counts = counts,
      sizeFactors = sizeFactors,
      normalizationFactors = normalizationFactors
    )
    if (!is.null(normalized)) {
      return(normalized)
    }
  }
  .rsdeseq2_normalized_counts_fallback(
    counts = counts,
    sizeFactors = sizeFactors,
    normalizationFactors = normalizationFactors
  )
}

.rsdeseq2_normalized_counts_fallback <- function(counts,
                                                 sizeFactors,
                                                 normalizationFactors) {
  if (!is.null(normalizationFactors)) {
    normalized <- counts / normalizationFactors
  } else {
    normalized <- sweep(counts, 2L, sizeFactors, "/")
  }
  dimnames(normalized) <- dimnames(counts)
  normalized
}

baseMeanRust <- function(counts,
                         sizeFactors = NULL,
                         normalizationFactors = NULL,
                         native = FALSE) {
  counts <- .rsdeseq2_validate_counts(counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (is.null(normalizationFactors)) {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  } else {
    sizeFactors <- NULL
  }
  if (isTRUE(native)) {
    means <- .rsdeseq2_try_native_base_mean(
      counts = counts,
      sizeFactors = sizeFactors,
      normalizationFactors = normalizationFactors
    )
    if (!is.null(means)) {
      return(means)
    }
  }
  .rsdeseq2_base_mean_fallback(
    counts,
    sizeFactors = sizeFactors,
    normalizationFactors = normalizationFactors
  )
}

.rsdeseq2_base_mean_fallback <- function(counts,
                                         sizeFactors,
                                         normalizationFactors) {
  normalized <- normalizedCountsRust(
    counts,
    sizeFactors = sizeFactors,
    normalizationFactors = normalizationFactors
  )
  means <- rowMeans(normalized)
  names(means) <- rownames(normalized)
  means
}

baseMetadataRust <- function(counts,
                             sizeFactors = NULL,
                             normalizationFactors = NULL,
                             weights = NULL,
                             native = FALSE) {
  counts <- .rsdeseq2_validate_counts(counts)
  weights <- .rsdeseq2_validate_observation_weights(weights, counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (is.null(normalizationFactors)) {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  } else {
    sizeFactors <- NULL
  }
  if (isTRUE(native)) {
    nativeMetadata <- .rsdeseq2_try_native_base_metadata(
      counts = counts,
      sizeFactors = sizeFactors,
      normalizationFactors = normalizationFactors,
      weights = weights
    )
    if (!is.null(nativeMetadata)) {
      return(nativeMetadata)
    }
  }
  .rsdeseq2_base_metadata_fallback(
    counts = counts,
    sizeFactors = sizeFactors,
    normalizationFactors = normalizationFactors,
    weights = weights
  )
}

.rsdeseq2_base_metadata_fallback <- function(counts,
                                             sizeFactors,
                                             normalizationFactors,
                                             weights) {
  normalized <- normalizedCountsRust(
    counts,
    sizeFactors = sizeFactors,
    normalizationFactors = normalizationFactors
  )
  metadataCounts <- if (is.null(weights)) {
    normalized
  } else {
    normalized * weights
  }
  baseVar <- apply(metadataCounts, 1L, stats::var)
  out <- data.frame(
    baseMean = rowMeans(metadataCounts),
    baseVar = as.numeric(baseVar),
    allZero = rowSums(counts) == 0,
    check.names = FALSE
  )
  if (!is.null(rownames(counts))) {
    rownames(out) <- rownames(counts)
  }
  out
}

.rsdeseq2_try_native_base_metadata <- function(counts,
                                               sizeFactors,
                                               normalizationFactors,
                                               weights) {
  args <- list(
    counts,
    sizeFactors %||% NULL,
    normalizationFactors %||% NULL,
    weights %||% NULL
  )
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_base_metadata", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_base_metadata"), args)),
    error = function(error) NULL
  )
}

.rsdeseq2_try_native_normalized_counts <- function(counts,
                                                   sizeFactors,
                                                   normalizationFactors) {
  args <- list(
    counts,
    sizeFactors %||% NULL,
    normalizationFactors %||% NULL
  )
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_normalized_counts", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_normalized_counts"), args)),
    error = function(error) NULL
  )
}

.rsdeseq2_try_native_base_mean <- function(counts,
                                           sizeFactors,
                                           normalizationFactors) {
  args <- list(
    counts,
    sizeFactors %||% NULL,
    normalizationFactors %||% NULL
  )
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_base_mean", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_base_mean"), args)),
    error = function(error) NULL
  )
}
