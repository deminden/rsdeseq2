estimateSizeFactorsRust <- function(counts,
                                    method = c("ratio", "poscounts"),
                                    geoMeans = NULL,
                                    controlGenes = NULL,
                                    native = TRUE) {
  method <- match.arg(method)
  counts <- .rsdeseq2_validate_counts(counts)
  controlIdx <- .rsdeseq2_control_gene_indices(controlGenes, counts)
  if (!is.null(geoMeans)) {
    .rsdeseq2_validate_geo_means(geoMeans, nrow(counts))
    geoMeans <- as.double(geoMeans)
  }
  if (is.null(controlIdx)) controlIdx <- seq_len(nrow(counts))
  sizeFactors <- .Call(
    "rsdeseq2_estimate_size_factors", counts, geoMeans, as.integer(controlIdx),
    as.integer(match(method, c("ratio", "poscounts")) - 1L), PACKAGE = "rsdeseq2"
  )
  if (!is.null(colnames(counts))) names(sizeFactors) <- colnames(counts)
  sizeFactors
}
