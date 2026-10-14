normalizedCountsRust <- function(counts,
                                 sizeFactors = NULL,
                                 normalizationFactors = NULL,
                                 native = TRUE) {
  counts <- .rsdeseq2_validate_counts(counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (!is.null(normalizationFactors)) {
    sizeFactors <- NULL
  } else {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  }
  .Call(
    "rsdeseq2_normalized_counts", counts, sizeFactors, normalizationFactors,
    PACKAGE = "rsdeseq2"
  )
}

baseMeanRust <- function(counts,
                         sizeFactors = NULL,
                         normalizationFactors = NULL,
                         native = TRUE) {
  counts <- .rsdeseq2_validate_counts(counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (is.null(normalizationFactors)) {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  } else {
    sizeFactors <- NULL
  }
  .Call(
    "rsdeseq2_base_mean", counts, sizeFactors, normalizationFactors,
    PACKAGE = "rsdeseq2"
  )
}

baseMetadataRust <- function(counts,
                             sizeFactors = NULL,
                             normalizationFactors = NULL,
                             weights = NULL,
                             native = TRUE) {
  counts <- .rsdeseq2_validate_counts(counts)
  weights <- .rsdeseq2_validate_observation_weights(weights, counts)
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  if (is.null(normalizationFactors)) {
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  } else {
    sizeFactors <- NULL
  }
  .Call(
    "rsdeseq2_base_metadata", counts, sizeFactors, normalizationFactors, weights,
    PACKAGE = "rsdeseq2"
  )
}
