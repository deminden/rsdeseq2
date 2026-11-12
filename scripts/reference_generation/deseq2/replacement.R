cooks_replacement_counts <- matrix(
  c(
    10L, 30L, 20L, 50L,
    0L, 0L, 0L, 0L,
    5L, 10L, 15L, 20L,
    40L, 80L, 60L, 100L
  ),
  nrow = 4,
  byrow = TRUE
)
rownames(cooks_replacement_counts) <- paste0("replace_gene", seq_len(nrow(cooks_replacement_counts)))
colnames(cooks_replacement_counts) <- paste0("sample", seq_len(ncol(cooks_replacement_counts)))
cooks_replacement_col_data <- S4Vectors::DataFrame(
  condition = factor(rep("A", ncol(cooks_replacement_counts))),
  row.names = colnames(cooks_replacement_counts)
)
cooks_replacement_design <- stats::model.matrix(~ 1, data = as.data.frame(cooks_replacement_col_data))
cooks_replacement_size_factors <- c(1, 2, 1, 2)
names(cooks_replacement_size_factors) <- colnames(cooks_replacement_counts)
cooks_replacement_cooks <- matrix(
  c(
    0, 9, 1, 0.5,
    8, 0, 0, 0,
    0, 0, 7, 0,
    0, 0, 0, 0
  ),
  nrow = nrow(cooks_replacement_counts),
  byrow = TRUE,
  dimnames = dimnames(cooks_replacement_counts)
)
cooks_replacement_trim <- 0.25
cooks_replacement_cutoff <- 5
cooks_replacement_min_replicates <- 3L
cooks_replacement_which_samples <- c(FALSE, TRUE, TRUE, FALSE)

write_matrix_tsv(
  cooks_replacement_counts,
  "gene",
  file.path(data_dir, "cooks_replacement_counts.tsv")
)
write_matrix_tsv(
  cooks_replacement_design,
  "sample",
  file.path(data_dir, "cooks_replacement_design.tsv")
)
write_matrix_tsv(
  cooks_replacement_cooks,
  "gene",
  file.path(data_dir, "cooks_replacement_cooks.tsv")
)
write_tsv(
  data.frame(
    sample = names(cooks_replacement_size_factors),
    size_factor = cooks_replacement_size_factors,
    replaceable = cooks_replacement_which_samples,
    check.names = FALSE
  ),
  file.path(data_dir, "cooks_replacement_size_factors.tsv")
)
write_tsv(
  data.frame(
    trim = cooks_replacement_trim,
    cooksCutoff = cooks_replacement_cutoff,
    minReplicates = cooks_replacement_min_replicates,
    check.names = FALSE
  ),
  file.path(data_dir, "cooks_replacement_options.tsv")
)

cooks_replacement_reference <- tryCatch({
  dds_replace <- DESeq2::DESeqDataSetFromMatrix(
    countData = cooks_replacement_counts,
    colData = cooks_replacement_col_data,
    design = ~ 1
  )
  DESeq2::sizeFactors(dds_replace) <- cooks_replacement_size_factors
  attr(dds_replace, "modelMatrix") <- cooks_replacement_design
  assays_replace <- SummarizedExperiment::assays(dds_replace, withDimnames = FALSE)
  assays_replace[["cooks"]] <- cooks_replacement_cooks
  dds_replace <- SummarizedExperiment::`assays<-`(
    dds_replace,
    withDimnames = FALSE,
    value = assays_replace
  )

  replaced <- DESeq2::replaceOutliers(
    dds_replace,
    trim = cooks_replacement_trim,
    cooksCutoff = cooks_replacement_cutoff,
    minReplicates = cooks_replacement_min_replicates,
    whichSamples = cooks_replacement_which_samples
  )
  replaced_meta <- S4Vectors::mcols(replaced)
  replaced_with_base <- DESeq2:::getBaseMeansAndVariances(replaced)
  replaced_base_meta <- S4Vectors::mcols(replaced_with_base)

  normalized <- DESeq2::counts(dds_replace, normalized = TRUE)
  trim_base_mean <- apply(normalized, 1, mean, trim = cooks_replacement_trim)
  candidate_counts <- as.integer(outer(trim_base_mean, cooks_replacement_size_factors, "*"))
  dim(candidate_counts) <- dim(cooks_replacement_counts)
  dimnames(candidate_counts) <- dimnames(cooks_replacement_counts)

  replace_cooks <- cooks_replacement_cooks
  replace_cooks[, cooks_replacement_which_samples] <- 0
  post_refit_max_cooks <- DESeq2:::recordMaxCooks(
    DESeq2::design(replaced),
    SummarizedExperiment::colData(replaced),
    cooks_replacement_design,
    replace_cooks,
    nrow(replaced)
  )

  replace_flags <- replaced_meta$replace
  new_all_zero <- replace_flags & replaced_base_meta$allZero
  refit_replace <- replace_flags & !replaced_base_meta$allZero

  list(
    replaced_counts = DESeq2::counts(replaced),
    original_counts = SummarizedExperiment::assays(replaced)[["originalCounts"]],
    candidate_counts = candidate_counts,
    replace_flags = replace_flags,
    replaceable = replaced$replaceable,
    base_meta = replaced_base_meta,
    new_all_zero = new_all_zero,
    refit_replace = refit_replace,
    post_refit_max_cooks = post_refit_max_cooks
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "cooks_replacement_reference_error.txt"))
  NULL
})

if (!is.null(cooks_replacement_reference)) {
  write_matrix_tsv(
    cooks_replacement_reference$candidate_counts,
    "gene",
    file.path(data_dir, "cooks_replacement_candidate_counts.tsv")
  )
  write_matrix_tsv(
    cooks_replacement_reference$replaced_counts,
    "gene",
    file.path(data_dir, "cooks_replacement_replaced_counts.tsv")
  )
  write_matrix_tsv(
    cooks_replacement_reference$original_counts,
    "gene",
    file.path(data_dir, "cooks_replacement_original_counts.tsv")
  )
  write_tsv(
    data.frame(
      gene = rownames(cooks_replacement_counts),
      replace = cooks_replacement_reference$replace_flags,
      allZero = cooks_replacement_reference$base_meta$allZero,
      baseMean = cooks_replacement_reference$base_meta$baseMean,
      baseVar = cooks_replacement_reference$base_meta$baseVar,
      newAllZero = cooks_replacement_reference$new_all_zero,
      refitReplace = cooks_replacement_reference$refit_replace,
      postRefitMaxCooks = cooks_replacement_reference$post_refit_max_cooks,
      check.names = FALSE
    ),
    file.path(data_dir, "cooks_replacement_rows.tsv")
  )
}
