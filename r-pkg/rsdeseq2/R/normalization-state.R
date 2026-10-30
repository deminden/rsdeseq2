# Base factors and formula offsets are distinct from the factors written to a fit.
.rsdeseq2_prepare_normalization <- function(object, counts, sizeFactors,
                                           normalizationFactors, offsets,
                                           explicitSize = FALSE, explicitFactors = FALSE) {
  state <- S4Vectors::metadata(object)$rsdeseq2Normalization
  if (!is.null(normalizationFactors)) {
    normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
    normalizationFactors <- .rsdeseq2_align_fit_observations(normalizationFactors, counts, "normalizationFactors")
  }
  if (!is.null(state)) {
    rows <- match(rownames(counts), state$geneNames)
    samples <- match(colnames(counts), state$sampleNames)
    if (anyNA(rows) || anyNA(samples)) {
      stop("Normalization provenance does not match current genes/samples; use a fresh dataset with explicit base factors", call. = FALSE)
    }
    for (name in c("baseFactors", "effectiveFactors")) {
      if (!is.null(state[[name]])) state[[name]] <- state[[name]][rows, samples, drop = FALSE]
    }
    for (name in c("baseSizeFactors", "writtenSizeFactors")) {
      if (!is.null(state[[name]])) state[[name]] <- state[[name]][samples]
    }
    if (!explicitFactors && identical(normalizationFactors, state$effectiveFactors)) {
      normalizationFactors <- state$baseFactors
    }
    if (!explicitSize && identical(unname(sizeFactors), unname(state$writtenSizeFactors))) sizeFactors <- state$baseSizeFactors
  } else if (!is.null(offsets) && !is.null(S4Vectors::metadata(object)$rsdeseq2RustFit) &&
             !is.null(normalizationFactors) && !explicitFactors) {
    stop("This saved fit has offsets without base-factor provenance; supply explicit base normalizationFactors or use the original unfitted dataset", call. = FALSE)
  }
  if (!is.null(sizeFactors)) {
    if (!is.null(names(sizeFactors)) && !is.null(colnames(counts))) {
      if (anyDuplicated(names(sizeFactors))) stop("sizeFactors names must be unique", call. = FALSE)
      order <- match(colnames(counts), names(sizeFactors))
      if (anyNA(order)) stop("sizeFactors names must match count columns", call. = FALSE)
      sizeFactors <- sizeFactors[order]
    }
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  }
  baseFactors <- normalizationFactors
  if (!is.null(offsets)) {
    if (length(offsets) != ncol(counts) || any(!is.finite(offsets))) stop("formula offsets must be finite with one value per sample", call. = FALSE)
    if (is.null(baseFactors) && is.null(sizeFactors)) sizeFactors <- estimateSizeFactorsRust(counts)
    normalizationFactors <- .Call("rsdeseq2_offset_factors", counts,
      sizeFactors, baseFactors, as.numeric(offsets), PACKAGE="rsdeseq2")
    normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  }
  list(sizeFactors = sizeFactors, normalizationFactors = normalizationFactors,
       state = list(version = 1L, baseFactors = baseFactors, baseSizeFactors = sizeFactors,
                    offsets = offsets, effectiveFactors = normalizationFactors,
                    geneNames = rownames(counts), sampleNames = colnames(counts)))
}
