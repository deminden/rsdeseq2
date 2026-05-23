#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required. Install it with BiocManager::install('DESeq2').")
  }
})

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath(file.path("scripts", "generate_deseq2_references.R"), mustWork = FALSE)
}

repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
out_dir <- file.path(repo_root, "results", "parity")
data_dir <- file.path(repo_root, "crates", "rsdeseq2", "tests", "data", "deseq2_reference")
if (dir.exists(data_dir)) {
  unlink(data_dir, recursive = TRUE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_matrix_tsv <- function(matrix_value, row_name_col, path) {
  row_ids <- rownames(matrix_value)
  if (is.null(row_ids)) {
    row_ids <- seq_len(nrow(matrix_value))
  }
  write_tsv(
    data.frame(
      setNames(data.frame(row_ids, check.names = FALSE), row_name_col),
      matrix_value,
      check.names = FALSE
    ),
    path
  )
}

counts <- matrix(
  c(
    10L, 12L, 20L, 24L,
    0L, 0L, 5L, 7L,
    100L, 80L, 90L, 120L,
    3L, 6L, 9L, 12L
  ),
  nrow = 4,
  byrow = TRUE
)
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("sample", seq_len(ncol(counts)))

col_data <- S4Vectors::DataFrame(
  condition = factor(c("A", "A", "B", "B")),
  row.names = colnames(counts)
)

full_design <- stats::model.matrix(~ condition, data = as.data.frame(col_data))
reduced_design <- stats::model.matrix(~ 1, data = as.data.frame(col_data))
rownames(full_design) <- rownames(col_data)
rownames(reduced_design) <- rownames(col_data)
fixed_dispersions <- c(0.10, 0.15, 0.05, 0.20)
names(fixed_dispersions) <- rownames(counts)
bioconductor_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

write_matrix_tsv(counts, "gene", file.path(data_dir, "counts.tsv"))
write_tsv(
  data.frame(
    sample = rownames(col_data),
    condition = as.character(col_data$condition),
    check.names = FALSE
  ),
  file.path(data_dir, "col_data.tsv")
)
write_matrix_tsv(full_design, "sample", file.path(data_dir, "design_full.tsv"))
write_matrix_tsv(reduced_design, "sample", file.path(data_dir, "design_reduced.tsv"))
write_tsv(
  data.frame(gene = names(fixed_dispersions), dispersion = fixed_dispersions),
  file.path(data_dir, "fixed_dispersions.tsv")
)

metadata <- data.frame(
  key = c(
    "DESeq2_version",
    "Bioconductor_version",
    "R_version",
    "platform",
    "reference_fixture",
    "fixed_glm_reference_mode",
    "normalization_factors_reference_mode",
    "native_nf_dispersion_reference_mode",
    "weighted_reference_mode",
    "native_weighted_glm_mu_cr_reference_mode",
    "native_weighted_glm_mu_reference_mode"
  ),
  value = c(
    as.character(utils::packageVersion("DESeq2")),
    bioconductor_version,
    paste(R.version$major, R.version$minor, sep = "."),
    R.version$platform,
    "tiny_two_group_four_gene",
    "DESeq2:::fitNbinomGLMs with supplied dispersions, default 1e-6 beta ridge, useQR=FALSE, useOptim=FALSE",
    "counts(dds, normalized=TRUE) with normalizationFactors(dds), which preempt sizeFactors(dds)",
    "DESeq2 roughDispEstimate, momentsDispEstimate, and estimateDispersionsGeneEst stored mu with normalizationFactors(dds)",
    "getBaseMeansAndVariances with raw weights, getAndCheckWeights row-normalized weights, and fitNbinomGLMs with supplied dispersions",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE) with observation weights for weighted Cox-Reid gene-wise dispersion anchors",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE), estimateDispersionsFit(fitType='mean'), estimateDispersionsMAP(useCR=FALSE), and full/reduced fitNbinomGLMs(useQR=FALSE,useOptim=FALSE) with observation weights"
  )
)
write_tsv(metadata, file.path(data_dir, "metadata.tsv"))

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

native_weighted_glm_mu_reference <- tryCatch({
  dds_weighted_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_weighted,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = FALSE,
    quiet = TRUE
  )
  dds_weighted_fit <- DESeq2:::estimateDispersionsFit(
    dds_weighted_gene,
    fitType = "mean",
    quiet = TRUE
  )
  disp_prior_var <- DESeq2:::estimateDispersionsPriorVar(
    dds_weighted_fit,
    modelMatrix = full_design
  )
  dds_weighted_map <- DESeq2:::estimateDispersionsMAP(
    dds_weighted_fit,
    dispPriorVar = disp_prior_var,
    modelMatrix = full_design,
    useCR = FALSE,
    quiet = TRUE
  )

  row_meta <- S4Vectors::mcols(dds_weighted_map)
  weights_fail <- if ("weightsFail" %in% names(row_meta)) {
    row_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  all_zero <- row_meta$allZero
  object_nz <- dds_weighted_map[!all_zero, , drop = FALSE]
  full_fit_nz <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = full_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))
  reduced_fit_nz <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = reduced_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))

  expand_vector <- function(value, missing_rows) {
    expanded <- rep(NA_real_, length(missing_rows))
    expanded[!missing_rows] <- value
    expanded
  }
  expand_logical <- function(value, missing_rows) {
    expanded <- rep(NA, length(missing_rows))
    expanded[!missing_rows] <- value
    expanded
  }
  expand_matrix <- function(value, missing_rows, dim_names) {
    expanded <- DESeq2:::buildMatrixWithNARows(value, missing_rows)
    dimnames(expanded) <- dim_names
    expanded
  }

  coefficient <- 2L
  wald_stat_nz <- full_fit_nz$betaMatrix[, coefficient] / full_fit_nz$betaSE[, coefficient]
  wald_pvalue_nz <- 2 * stats::pnorm(abs(wald_stat_nz), lower.tail = FALSE)
  lrt_stat_nz <- 2 * (full_fit_nz$logLike - reduced_fit_nz$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue_nz <- stats::pchisq(lrt_stat_nz, df = lrt_df, lower.tail = FALSE)

  list(
    row_meta = row_meta,
    weights_fail = weights_fail,
    disp_prior_var = disp_prior_var,
    var_log_disp_estimates = attr(DESeq2::dispersionFunction(dds_weighted_map), "varLogDispEsts"),
    dispersion_mu = SummarizedExperiment::assays(dds_weighted_map)[["mu"]],
    beta_intercept = expand_vector(full_fit_nz$betaMatrix[, 1], all_zero),
    beta_conditionB = expand_vector(full_fit_nz$betaMatrix[, 2], all_zero),
    beta_se_intercept = expand_vector(full_fit_nz$betaSE[, 1], all_zero),
    beta_se_conditionB = expand_vector(full_fit_nz$betaSE[, 2], all_zero),
    stat_conditionB = expand_vector(wald_stat_nz, all_zero),
    pvalue_conditionB = expand_vector(wald_pvalue_nz, all_zero),
    lrt_stat = expand_vector(lrt_stat_nz, all_zero),
    lrt_df = lrt_df,
    lrt_pvalue = expand_vector(lrt_pvalue_nz, all_zero),
    log_like = expand_vector(full_fit_nz$logLike, all_zero),
    reduced_log_like = expand_vector(reduced_fit_nz$logLike, all_zero),
    beta_converged = expand_logical(full_fit_nz$betaConv, all_zero),
    beta_iterations = expand_vector(full_fit_nz$betaIter, all_zero),
    reduced_beta_converged = expand_logical(reduced_fit_nz$betaConv, all_zero),
    reduced_beta_iterations = expand_vector(reduced_fit_nz$betaIter, all_zero),
    wald_mu = expand_matrix(full_fit_nz$mu, all_zero, dimnames(counts)),
    wald_hat = expand_matrix(full_fit_nz$hat_diagonals, all_zero, dimnames(counts))
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_reference)) {
  native_weighted_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_reference$row_meta$dispFit,
    dispPriorVar = native_weighted_glm_mu_reference$disp_prior_var,
    varLogDispEsts = native_weighted_glm_mu_reference$var_log_disp_estimates,
    dispMAP = native_weighted_glm_mu_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_reference$row_meta$dispersion,
    dispIter = native_weighted_glm_mu_reference$row_meta$dispIter,
    dispOutlier = native_weighted_glm_mu_reference$row_meta$dispOutlier,
    beta_intercept = native_weighted_glm_mu_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_reference$beta_se_conditionB,
    stat_conditionB = native_weighted_glm_mu_reference$stat_conditionB,
    pvalue_conditionB = native_weighted_glm_mu_reference$pvalue_conditionB,
    log_like = native_weighted_glm_mu_reference$log_like,
    beta_converged = native_weighted_glm_mu_reference$beta_converged,
    beta_iterations = native_weighted_glm_mu_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_rows,
    file.path(data_dir, "native_weighted_glm_mu_reference.tsv")
  )
  native_weighted_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_reference$row_meta$dispFit,
    dispMAP = native_weighted_glm_mu_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_reference$row_meta$dispersion,
    beta_intercept = native_weighted_glm_mu_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_reference$beta_se_conditionB,
    log_like_full = native_weighted_glm_mu_reference$log_like,
    log_like_reduced = native_weighted_glm_mu_reference$reduced_log_like,
    lrt_stat = native_weighted_glm_mu_reference$lrt_stat,
    pvalue = native_weighted_glm_mu_reference$lrt_pvalue,
    df = native_weighted_glm_mu_reference$lrt_df,
    full_converged = native_weighted_glm_mu_reference$beta_converged,
    reduced_converged = native_weighted_glm_mu_reference$reduced_beta_converged,
    full_iterations = native_weighted_glm_mu_reference$beta_iterations,
    reduced_iterations = native_weighted_glm_mu_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_lrt_rows,
    file.path(data_dir, "native_weighted_glm_mu_lrt_reference.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_wald_hat.tsv")
  )
}
trend_means <- c(10, 20, 40, 80, 160, 320)
trend_disps <- 0.05 + 2 / trend_means
trend_fit <- DESeq2:::parametricDispersionFit(trend_means, trend_disps)
trend_coefs <- attr(trend_fit, "coefficients")
write_tsv(
  data.frame(
    gene = paste0("trend_gene", seq_along(trend_means)),
    baseMean = trend_means,
    dispGeneEst = trend_disps,
    dispFit = trend_fit(trend_means),
    useForFit = trend_disps > 100 * 1e-8,
    asymptDisp = unname(trend_coefs[["asymptDisp"]]),
    extraPois = unname(trend_coefs[["extraPois"]]),
    check.names = FALSE
  ),
  file.path(data_dir, "parametric_trend_reference.tsv")
)

prior_residuals <- c(-2, -1, 0, 1, 2)
prior_disp_fit <- rep(1, length(prior_residuals))
prior_disp_gene_est <- exp(prior_residuals) * prior_disp_fit
prior_min_disp <- 1e-8
prior_n_samples <- 10L
prior_n_coefficients <- 2L
prior_residual_df <- prior_n_samples - prior_n_coefficients
prior_above_min <- prior_disp_gene_est >= 100 * prior_min_disp
prior_var_log <- mad(log(prior_disp_gene_est) - log(prior_disp_fit), na.rm = TRUE)^2
prior_expected_sampling <- trigamma(prior_residual_df / 2)
prior_var <- pmax(prior_var_log - prior_expected_sampling, 0.25)
write_tsv(
  data.frame(
    gene = paste0("prior_gene", seq_along(prior_residuals)),
    dispGeneEst = prior_disp_gene_est,
    dispFit = prior_disp_fit,
    aboveMinDisp = prior_above_min,
    nSamples = prior_n_samples,
    nCoefficients = prior_n_coefficients,
    residualDf = prior_residual_df,
    varLogDispEsts = prior_var_log,
    expectedLogDispVariance = prior_expected_sampling,
    dispPriorVar = prior_var,
    check.names = FALSE
  ),
  file.path(data_dir, "dispersion_prior_variance_reference.tsv")
)

map_counts <- matrix(c(10L, 30L, 10L, 30L), nrow = 1)
rownames(map_counts) <- "map_gene1"
colnames(map_counts) <- colnames(counts)
map_mu <- matrix(rep(20, 4), nrow = 1)
map_disp_gene_est <- 0.5
map_disp_fit <- 0.02
map_prior_var <- 0.05
map_var_log_disp_est <- 100
map_disp_init <- ifelse(map_disp_gene_est > 0.1 * map_disp_fit, map_disp_gene_est, map_disp_fit)
map_fit <- DESeq2:::fitDisp(
  ySEXP = map_counts,
  xSEXP = full_design,
  mu_hatSEXP = map_mu,
  log_alphaSEXP = log(map_disp_init),
  log_alpha_prior_meanSEXP = log(map_disp_fit),
  log_alpha_prior_sigmasqSEXP = map_prior_var,
  min_log_alphaSEXP = log(1e-8 / 10),
  kappa_0SEXP = 1,
  tolSEXP = 1e-6,
  maxitSEXP = 100,
  usePriorSEXP = TRUE,
  weightsSEXP = matrix(1, nrow = 1, ncol = 4),
  useWeightsSEXP = FALSE,
  weightThresholdSEXP = 1e-2,
  useCRSEXP = FALSE
)
map_disp_map <- exp(map_fit$log_alpha)
map_disp_outlier <- log(map_disp_gene_est) > log(map_disp_fit) + 2 * sqrt(map_var_log_disp_est)
map_dispersion <- ifelse(map_disp_outlier, map_disp_gene_est, map_disp_map)
write_tsv(
  data.frame(
    gene = "map_gene1",
    dispGeneEst = map_disp_gene_est,
    dispFit = map_disp_fit,
    dispInit = map_disp_init,
    dispPriorVar = map_prior_var,
    varLogDispEsts = map_var_log_disp_est,
    dispMAP = map_disp_map,
    dispersion = map_dispersion,
    dispIter = map_fit$iter,
    dispOutlier = map_disp_outlier,
    converged = map_fit$iter < 100,
    useCR = FALSE,
    check.names = FALSE
  ),
  file.path(data_dir, "map_dispersion_reference.tsv")
)

full_results <- tryCatch({
  dds_full <- DESeq2::DESeq(dds_ratio, quiet = TRUE)
  DESeq2::results(dds_full)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "results_wald_ratio_error.txt"))
  NULL
})

if (!is.null(full_results)) {
  write_tsv(
    data.frame(gene = rownames(full_results), as.data.frame(full_results), check.names = FALSE),
    file.path(data_dir, "results_wald_ratio.tsv")
  )
}

fixed_reference <- tryCatch({
  dds_fixed <- dds_ratio
  dds_fixed <- DESeq2::`dispersions<-`(dds_fixed, value = fixed_dispersions)
  mcols_fixed <- S4Vectors::mcols(dds_fixed)
  mcols_fixed$dispersion <- fixed_dispersions
  mcols_fixed$allZero <- rowSums(counts) == 0
  dds_fixed <- S4Vectors::`mcols<-`(dds_fixed, value = mcols_fixed)

  full_fit <- DESeq2:::fitNbinomGLMs(
    dds_fixed,
    modelMatrix = full_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )
  reduced_fit <- DESeq2:::fitNbinomGLMs(
    dds_fixed,
    modelMatrix = reduced_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )

  coefficient <- 2L
  wald_stat <- full_fit$betaMatrix[, coefficient] / full_fit$betaSE[, coefficient]
  wald_pvalue <- 2 * stats::pnorm(abs(wald_stat), lower.tail = FALSE)
  wald_t_df <- ncol(counts) - ncol(full_design)
  wald_t_pvalue <- 2 * stats::pt(abs(wald_stat), df = wald_t_df, lower.tail = FALSE)
  lrt_stat <- 2 * (full_fit$logLike - reduced_fit$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  assays_fixed <- SummarizedExperiment::assays(dds_fixed, withDimnames = FALSE)
  assays_fixed[["mu"]] <- full_fit$mu
  dds_fixed <- SummarizedExperiment::`assays<-`(
    dds_fixed,
    withDimnames = FALSE,
    value = assays_fixed
  )
  cooks <- DESeq2:::calculateCooksDistance(dds_fixed, full_fit$hat_diagonals, full_design)

  list(
    full_fit = full_fit,
    reduced_fit = reduced_fit,
    wald_stat = wald_stat,
    wald_pvalue = wald_pvalue,
    wald_t_df = wald_t_df,
    wald_t_pvalue = wald_t_pvalue,
    lrt_stat = lrt_stat,
    lrt_df = lrt_df,
    lrt_pvalue = lrt_pvalue,
    cooks = cooks
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "fixed_glm_reference_error.txt"))
  NULL
})

if (!is.null(fixed_reference)) {
  full_fit <- fixed_reference$full_fit
  reduced_fit <- fixed_reference$reduced_fit

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      beta_intercept = full_fit$betaMatrix[, 1],
      beta_conditionB = full_fit$betaMatrix[, 2],
      beta_se_intercept = full_fit$betaSE[, 1],
      beta_se_conditionB = full_fit$betaSE[, 2],
      stat_conditionB = fixed_reference$wald_stat,
      pvalue_conditionB = fixed_reference$wald_pvalue,
      log_like = full_fit$logLike,
      converged = full_fit$betaConv,
      iterations = full_fit$betaIter,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      stat_conditionB = fixed_reference$wald_stat,
      df_conditionB = fixed_reference$wald_t_df,
      pvalue_conditionB = fixed_reference$wald_t_pvalue,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_t_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      beta_intercept = full_fit$betaMatrix[, 1],
      beta_conditionB = full_fit$betaMatrix[, 2],
      beta_se_intercept = full_fit$betaSE[, 1],
      beta_se_conditionB = full_fit$betaSE[, 2],
      log_like_full = full_fit$logLike,
      log_like_reduced = reduced_fit$logLike,
      lrt_stat = fixed_reference$lrt_stat,
      pvalue = fixed_reference$lrt_pvalue,
      df = fixed_reference$lrt_df,
      full_converged = full_fit$betaConv,
      reduced_converged = reduced_fit$betaConv,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_lrt_reference.tsv")
  )

  write_matrix_tsv(full_fit$mu, "gene", file.path(data_dir, "fixed_mu_full.tsv"))
  write_matrix_tsv(full_fit$hat_diagonals, "gene", file.path(data_dir, "fixed_hat_full.tsv"))
  write_matrix_tsv(fixed_reference$cooks, "gene", file.path(data_dir, "fixed_cooks_full.tsv"))
}

weighted_fixed_reference <- tryCatch({
  dds_weighted_fixed <- DESeq2:::getBaseMeansAndVariances(dds_weighted)
  dds_weighted_fixed <- DESeq2::`dispersions<-`(dds_weighted_fixed, value = fixed_dispersions)
  mcols_weighted_fixed <- S4Vectors::mcols(dds_weighted_fixed)
  mcols_weighted_fixed$dispersion <- fixed_dispersions
  dds_weighted_fixed <- S4Vectors::`mcols<-`(dds_weighted_fixed, value = mcols_weighted_fixed)
  checked_weighted_fixed <- suppressWarnings(DESeq2:::getAndCheckWeights(
    dds_weighted_fixed,
    full_design
  )$object)
  fixed_all_zero <- S4Vectors::mcols(checked_weighted_fixed)$allZero
  object_nz <- checked_weighted_fixed[!fixed_all_zero, , drop = FALSE]

  full_fit <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = full_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))
  reduced_fit <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = reduced_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))

  coefficient <- 2L
  wald_stat <- full_fit$betaMatrix[, coefficient] / full_fit$betaSE[, coefficient]
  wald_pvalue <- 2 * stats::pnorm(abs(wald_stat), lower.tail = FALSE)
  lrt_stat <- 2 * (full_fit$logLike - reduced_fit$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  expand_numeric <- function(value) {
    expanded <- rep(NA_real_, length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }
  expand_integer <- function(value) {
    expanded <- integer(length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }
  expand_logical <- function(value) {
    expanded <- rep(FALSE, length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }

  weighted_meta <- S4Vectors::mcols(checked_weighted_fixed)
  weights_fail <- if ("weightsFail" %in% names(weighted_meta)) {
    weighted_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }

  list(
    full_fit = full_fit,
    reduced_fit = reduced_fit,
    wald_stat = expand_numeric(wald_stat),
    wald_pvalue = expand_numeric(wald_pvalue),
    lrt_stat = expand_numeric(lrt_stat),
    lrt_df = lrt_df,
    lrt_pvalue = expand_numeric(lrt_pvalue),
    row_meta = weighted_meta,
    weights_fail = weights_fail,
    beta_intercept = expand_numeric(full_fit$betaMatrix[, 1]),
    beta_conditionB = expand_numeric(full_fit$betaMatrix[, 2]),
    beta_se_intercept = expand_numeric(full_fit$betaSE[, 1]),
    beta_se_conditionB = expand_numeric(full_fit$betaSE[, 2]),
    log_like = expand_numeric(full_fit$logLike),
    converged = expand_logical(full_fit$betaConv),
    iterations = expand_integer(full_fit$betaIter),
    reduced_log_like = expand_numeric(reduced_fit$logLike),
    reduced_converged = expand_logical(reduced_fit$betaConv),
    reduced_iterations = expand_integer(reduced_fit$betaIter)
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "fixed_weighted_glm_reference_error.txt"))
  NULL
})

if (!is.null(weighted_fixed_reference)) {
  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      baseMean = weighted_fixed_reference$row_meta$baseMean,
      allZero = weighted_fixed_reference$row_meta$allZero,
      weightsFail = weighted_fixed_reference$weights_fail,
      beta_intercept = weighted_fixed_reference$beta_intercept,
      beta_conditionB = weighted_fixed_reference$beta_conditionB,
      beta_se_intercept = weighted_fixed_reference$beta_se_intercept,
      beta_se_conditionB = weighted_fixed_reference$beta_se_conditionB,
      stat_conditionB = weighted_fixed_reference$wald_stat,
      pvalue_conditionB = weighted_fixed_reference$wald_pvalue,
      log_like = weighted_fixed_reference$log_like,
      converged = weighted_fixed_reference$converged,
      iterations = weighted_fixed_reference$iterations,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_weighted_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      baseMean = weighted_fixed_reference$row_meta$baseMean,
      allZero = weighted_fixed_reference$row_meta$allZero,
      weightsFail = weighted_fixed_reference$weights_fail,
      beta_intercept = weighted_fixed_reference$beta_intercept,
      beta_conditionB = weighted_fixed_reference$beta_conditionB,
      beta_se_intercept = weighted_fixed_reference$beta_se_intercept,
      beta_se_conditionB = weighted_fixed_reference$beta_se_conditionB,
      log_like_full = weighted_fixed_reference$log_like,
      log_like_reduced = weighted_fixed_reference$reduced_log_like,
      lrt_stat = weighted_fixed_reference$lrt_stat,
      pvalue = weighted_fixed_reference$lrt_pvalue,
      df = weighted_fixed_reference$lrt_df,
      full_converged = weighted_fixed_reference$converged,
      reduced_converged = weighted_fixed_reference$reduced_converged,
      full_iterations = weighted_fixed_reference$iterations,
      reduced_iterations = weighted_fixed_reference$reduced_iterations,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_lrt_weighted_reference.tsv")
  )
}
session_info <- sub("[[:blank:]]+$", "", capture.output(sessionInfo()))
writeLines(session_info, file.path(data_dir, "sessionInfo.txt"))
invisible(file.copy(file.path(data_dir, "metadata.tsv"), file.path(out_dir, "metadata.tsv"), overwrite = TRUE))

message("Wrote DESeq2 references to: ", data_dir)
