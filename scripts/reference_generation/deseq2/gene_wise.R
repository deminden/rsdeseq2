native_nf_reference <- tryCatch({
  dds_nf_meta <- DESeq2:::getBaseMeansAndVariances(dds_norm_factors)
  row_meta <- S4Vectors::mcols(dds_nf_meta)
  nonzero_rows <- !row_meta$allZero
  dds_nf_nz <- dds_nf_meta[nonzero_rows, , drop = FALSE]
  rough_disp <- DESeq2:::roughDispEstimate(
    y = DESeq2::counts(dds_nf_nz, normalized = TRUE),
    x = full_design
  )
  moments_disp <- DESeq2:::momentsDispEstimate(dds_nf_nz)
  min_disp <- 1e-8
  max_disp <- max(10, ncol(dds_nf_meta))
  disp_init <- pmin(pmax(min_disp, pmin(rough_disp, moments_disp)), max_disp)
  dds_nf_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_norm_factors,
    modelMatrix = full_design,
    linearMu = TRUE,
    niter = 1,
    useCR = FALSE,
    quiet = TRUE
  )
  gene_meta <- S4Vectors::mcols(dds_nf_gene)
  gene_mu <- SummarizedExperiment::assays(dds_nf_gene)[["mu"]]

  list(
    row_meta = row_meta,
    nonzero_rows = nonzero_rows,
    rough_disp = rough_disp,
    moments_disp = moments_disp,
    disp_init = disp_init,
    mu = gene_mu,
    disp_gene_est = gene_meta$dispGeneEst,
    disp_gene_iter = gene_meta$dispGeneIter
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_nf_dispersion_reference_error.txt"))
  NULL
})

if (!is.null(native_nf_reference)) {
  nonzero_index <- 0L
  native_rows <- lapply(seq_len(nrow(counts)), function(idx) {
    if (native_nf_reference$nonzero_rows[[idx]]) {
      nonzero_index <<- nonzero_index + 1L
      rough <- native_nf_reference$rough_disp[[nonzero_index]]
      moments <- native_nf_reference$moments_disp[[nonzero_index]]
      disp_init <- native_nf_reference$disp_init[[nonzero_index]]
    } else {
      rough <- NA_real_
      moments <- NA_real_
      disp_init <- NA_real_
    }
    data.frame(
      gene = rownames(counts)[[idx]],
      baseMean = native_nf_reference$row_meta$baseMean[[idx]],
      baseVar = native_nf_reference$row_meta$baseVar[[idx]],
      allZero = native_nf_reference$row_meta$allZero[[idx]],
      roughDisp = rough,
      momentsDisp = moments,
      dispInit = disp_init,
      dispGeneEst = native_nf_reference$disp_gene_est[[idx]],
      dispGeneIter = native_nf_reference$disp_gene_iter[[idx]],
      check.names = FALSE
    )
  })
  write_tsv(
    do.call(rbind, native_rows),
    file.path(data_dir, "native_nf_dispersion_reference.tsv")
  )
  write_matrix_tsv(
    native_nf_reference$mu,
    "gene",
    file.path(data_dir, "native_nf_mu.tsv")
  )
}

native_glm_mu_cr_reference <- tryCatch({
  dds_cr_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_ratio,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = TRUE,
    quiet = TRUE
  )
  list(
    row_meta = S4Vectors::mcols(dds_cr_gene),
    dispersion_mu = SummarizedExperiment::assays(dds_cr_gene)[["mu"]]
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_glm_mu_cr_reference_error.txt"))
  NULL
})

if (!is.null(native_glm_mu_cr_reference)) {
  native_cr_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_glm_mu_cr_reference$row_meta$baseMean,
    baseVar = native_glm_mu_cr_reference$row_meta$baseVar,
    allZero = native_glm_mu_cr_reference$row_meta$allZero,
    dispGeneEst = native_glm_mu_cr_reference$row_meta$dispGeneEst,
    dispGeneIter = native_glm_mu_cr_reference$row_meta$dispGeneIter,
    check.names = FALSE
  )
  write_tsv(
    native_cr_rows,
    file.path(data_dir, "native_glm_mu_cr_reference.tsv")
  )
  write_matrix_tsv(
    native_glm_mu_cr_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_glm_mu_cr_dispersion_mu.tsv")
  )
}

native_weighted_glm_mu_cr_reference <- tryCatch({
  dds_weighted_cr_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_weighted,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = TRUE,
    quiet = TRUE
  )
  row_meta <- S4Vectors::mcols(dds_weighted_cr_gene)
  weights_fail <- if ("weightsFail" %in% names(row_meta)) {
    row_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  list(
    row_meta = row_meta,
    weights_fail = weights_fail,
    dispersion_mu = SummarizedExperiment::assays(dds_weighted_cr_gene)[["mu"]]
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_cr_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_cr_reference)) {
  native_weighted_cr_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_cr_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_cr_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_cr_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_cr_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_cr_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_cr_reference$row_meta$dispGeneIter,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_cr_rows,
    file.path(data_dir, "native_weighted_glm_mu_cr_reference.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_cr_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_cr_dispersion_mu.tsv")
  )
}
