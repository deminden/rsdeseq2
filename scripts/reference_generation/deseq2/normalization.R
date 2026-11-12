dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = col_data,
  design = ~ condition
)

dds_ratio <- DESeq2::estimateSizeFactors(dds, type = "ratio")
dds_poscounts <- DESeq2::estimateSizeFactors(dds, type = "poscounts")

write_tsv(
  data.frame(
    sample = names(DESeq2::sizeFactors(dds_ratio)),
    size_factor = DESeq2::sizeFactors(dds_ratio)
  ),
  file.path(data_dir, "size_factors_ratio.tsv")
)

write_tsv(
  data.frame(
    sample = names(DESeq2::sizeFactors(dds_poscounts)),
    size_factor = DESeq2::sizeFactors(dds_poscounts)
  ),
  file.path(data_dir, "size_factors_poscounts.tsv")
)

normalized_ratio <- DESeq2::counts(dds_ratio, normalized = TRUE)
write_matrix_tsv(normalized_ratio, "gene", file.path(data_dir, "normalized_counts_ratio.tsv"))

base_mean <- rowMeans(normalized_ratio)
base_var <- apply(normalized_ratio, 1, stats::var)
all_zero <- rowSums(counts) == 0
write_tsv(
  data.frame(gene = names(base_mean), baseMean = base_mean),
  file.path(data_dir, "base_mean_ratio.tsv")
)

write_tsv(
  data.frame(
    gene = names(base_mean),
    baseMean = base_mean,
    baseVar = base_var,
    allZero = all_zero
  ),
  file.path(data_dir, "base_metadata_ratio.tsv")
)

normalization_factors <- matrix(
  c(
    1.0, 2.0, 1.0, 2.0,
    1.0, 1.5, 1.0, 1.5,
    2.0, 2.0, 2.0, 2.0,
    0.5, 1.0, 2.0, 4.0
  ),
  nrow = nrow(counts),
  byrow = TRUE,
  dimnames = dimnames(counts)
)
dds_norm_factors <- dds
DESeq2::normalizationFactors(dds_norm_factors) <- normalization_factors
normalized_nf <- DESeq2::counts(dds_norm_factors, normalized = TRUE)
base_mean_nf <- rowMeans(normalized_nf)
base_var_nf <- apply(normalized_nf, 1, stats::var)

write_matrix_tsv(
  normalization_factors,
  "gene",
  file.path(data_dir, "normalization_factors.tsv")
)
write_matrix_tsv(
  normalized_nf,
  "gene",
  file.path(data_dir, "normalized_counts_nf.tsv")
)
write_tsv(
  data.frame(
    gene = names(base_mean_nf),
    baseMean = base_mean_nf,
    baseVar = base_var_nf,
    allZero = all_zero
  ),
  file.path(data_dir, "base_metadata_nf.tsv")
)

observation_weights <- matrix(
  c(
    1.0, 2.0, 1.0, 2.0,
    1.0, 1.0, 0.0, 0.0,
    2.0, 1.0, 2.0, 1.0,
    1.0, 0.5, 1.0, 0.5
  ),
  nrow = nrow(counts),
  byrow = TRUE,
  dimnames = dimnames(counts)
)

dds_weighted <- dds_ratio
assays_weighted <- SummarizedExperiment::assays(dds_weighted, withDimnames = FALSE)
assays_weighted[["weights"]] <- observation_weights
dds_weighted <- SummarizedExperiment::`assays<-`(
  dds_weighted,
  withDimnames = FALSE,
  value = assays_weighted
)

weighted_metadata_reference <- tryCatch({
  dds_weighted_meta <- DESeq2:::getBaseMeansAndVariances(dds_weighted)
  wlist <- suppressWarnings(DESeq2:::getAndCheckWeights(dds_weighted_meta, full_design))
  weighted_meta <- S4Vectors::mcols(wlist$object)
  weights_fail <- if ("weightsFail" %in% names(weighted_meta)) {
    weighted_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  list(
    row_meta = weighted_meta,
    weights = wlist$weights,
    weights_fail = weights_fail
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "weighted_base_metadata_reference_error.txt"))
  NULL
})

write_matrix_tsv(
  observation_weights,
  "gene",
  file.path(data_dir, "observation_weights.tsv")
)

if (!is.null(weighted_metadata_reference)) {
  normalized_observation_weights <- weighted_metadata_reference$weights
  dimnames(normalized_observation_weights) <- dimnames(counts)
  write_matrix_tsv(
    normalized_observation_weights,
    "gene",
    file.path(data_dir, "observation_weights_normalized.tsv")
  )
  write_tsv(
    data.frame(
      gene = rownames(counts),
      baseMean = weighted_metadata_reference$row_meta$baseMean,
      baseVar = weighted_metadata_reference$row_meta$baseVar,
      allZero = weighted_metadata_reference$row_meta$allZero,
      weightsFail = weighted_metadata_reference$weights_fail,
      check.names = FALSE
    ),
    file.path(data_dir, "base_metadata_weighted.tsv")
  )
}
