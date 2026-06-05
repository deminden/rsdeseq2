#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- which(args == name)
  if (length(hit) == 0 || hit[length(hit)] == length(args)) {
    return(default)
  }
  args[hit[length(hit)] + 1]
}

study_root <- get_arg("--study-root", Sys.getenv("RSDESEQ2_REAL_DATA_ROOT", unset = NA_character_))
if (is.na(study_root) || !nzchar(study_root)) {
  stop("provide --study-root or set RSDESEQ2_REAL_DATA_ROOT", call. = FALSE)
}
diagnostics_path <- get_arg(
  "--diagnostics",
  "results/benchmarks/real_data_parity_non_lbfgsb_start_probe_diagnostics.tsv"
)
out_dir <- get_arg("--output-dir", "results/fixtures/se_covariance_hard_real_2026-06-05")
top_n <- as.integer(get_arg("--top-n", "64"))
min_se_abs <- as.numeric(get_arg("--min-se-abs", "1e-8"))
minmu <- as.numeric(get_arg("--minmu", "0.5"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  frame <- as.data.frame(x, check.names = FALSE)
  double_cols <- vapply(frame, is.double, logical(1))
  frame[double_cols] <- lapply(frame[double_cols], function(column) {
    formatted <- format(column, digits = 17, scientific = TRUE, trim = TRUE)
    formatted[is.na(column)] <- "NA"
    formatted
  })
  utils::write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

write_matrix_tsv <- function(matrix, path, row_name = "gene") {
  frame <- data.frame(row.names(matrix), matrix, check.names = FALSE)
  names(frame)[1] <- row_name
  write_tsv(frame, path)
}

read_counts <- function(path) {
  frame <- utils::read.delim(gzfile(path), check.names = FALSE)
  rownames(frame) <- frame[[1]]
  frame[[1]] <- NULL
  matrix <- as.matrix(frame)
  storage.mode(matrix) <- "integer"
  matrix
}

parse_contrast <- function(stem) {
  null_types <- c(
    "full_blocked_permutation",
    "blocked_permutation",
    "stratified_split",
    "random_split"
  )
  for (null_type in null_types) {
    marker <- paste0("_", null_type, "_rep")
    if (!grepl(marker, stem, fixed = TRUE)) {
      next
    }
    parts <- strsplit(stem, marker, fixed = TRUE)[[1]]
    return(list(tissue = parts[[1]], null_type = null_type, rep = as.integer(parts[[2]])))
  }
  stop("could not parse contrast stem: ", stem)
}

infer_design <- function(groups, null_type) {
  groups$condition <- factor(groups$condition, levels = c("A", "B"))
  groups$perm_block <- factor(groups$perm_block)
  full <- stats::model.matrix(~ perm_block + condition, groups)
  blocks <- levels(groups$perm_block)
  valid_blocks <- all(vapply(
    blocks,
    function(block) identical(sort(unique(as.character(groups$condition[groups$perm_block == block]))), c("A", "B")),
    logical(1)
  ))
  if (grepl("blocked_permutation", null_type) && valid_blocks && qr(full)$rank == ncol(full)) {
    return(list(design = ~ perm_block + condition, model_matrix = full, kind = "perm_block + condition"))
  }
  simple <- stats::model.matrix(~ condition, groups)
  list(design = ~ condition, model_matrix = simple, kind = "condition")
}

safe_float <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[is.nan(out)] <- NA_real_
  out
}

rank_rows <- function(diagnostics) {
  for (name in c("lfcSE_abs", "log2FoldChange_abs", "stat_abs", "rustBetaOptimIter", "refitRustBetaOptimIter")) {
    diagnostics[[name]] <- safe_float(diagnostics[[name]])
  }
  no_optim <- !is.finite(diagnostics$rustBetaOptimIter) &
    !is.finite(diagnostics$refitRustBetaOptimIter)
  no_replacement <- diagnostics$replace != "TRUE" & diagnostics$refitReplace != "TRUE"
  hard <- diagnostics[no_optim & no_replacement & is.finite(diagnostics$lfcSE_abs) &
    diagnostics$lfcSE_abs >= min_se_abs, , drop = FALSE]
  hard <- hard[order(-hard$lfcSE_abs, -hard$stat_abs), , drop = FALSE]
  hard[seq_len(min(nrow(hard), top_n)), , drop = FALSE]
}

diagnostics <- utils::read.delim(diagnostics_path, check.names = FALSE)
selected <- rank_rows(diagnostics)
if (nrow(selected) == 0) {
  stop("no non-optimizer SE/covariance rows passed the filters")
}

contrast <- unique(selected$contrast)
if (length(contrast) != 1) {
  stop("this focused generator expects one contrast; got: ", paste(contrast, collapse = ", "))
}
spec <- parse_contrast(contrast[[1]])

counts_path <- file.path(study_root, "01_inputs", sprintf("%s_raw_counts.tsv.gz", spec$tissue))
groups_path <- file.path(
  study_root,
  "02_null_splits",
  sprintf("%s_%s_rep%02d_groups.tsv", spec$tissue, spec$null_type, spec$rep)
)

counts_all <- read_counts(counts_path)
groups <- utils::read.delim(groups_path, check.names = FALSE)
groups <- groups[groups$retained == "TRUE" & groups$condition %in% c("A", "B"), , drop = FALSE]
counts <- counts_all[, groups$sample_id, drop = FALSE]
storage.mode(counts) <- "integer"

design_info <- infer_design(groups, spec$null_type)
col_data <- data.frame(
  row.names = groups$sample_id,
  condition = factor(groups$condition, levels = c("A", "B")),
  perm_block = factor(groups$perm_block),
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(counts, col_data, design = design_info$design)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
dds <- estimateDispersionsFit(dds, fitType = "parametric", quiet = TRUE)
dds <- estimateDispersionsMAP(dds, quiet = TRUE)
dispersion_function <- dispersionFunction(dds)
disp_prior_var <- attr(dispersion_function, "dispPriorVar")
var_log_disp_estimates <- attr(dispersion_function, "varLogDispEsts")

nonzero <- MatrixGenerics::rowSums(counts(dds)) > 0
dds_nz <- dds[nonzero, , drop = FALSE]
fit <- suppressWarnings(DESeq2:::fitNbinomGLMs(
  dds_nz,
  modelMatrix = design_info$model_matrix,
  useQR = FALSE,
  useOptim = FALSE,
  minmu = minmu
))
rownames(fit$betaMatrix) <- rownames(dds_nz)
rownames(fit$betaSE) <- rownames(dds_nz)
rownames(fit$mu) <- rownames(dds_nz)
rownames(fit$hat_diagonals) <- rownames(dds_nz)

genes <- selected$gene[selected$gene %in% rownames(dds_nz)]
if (length(genes) == 0) {
  stop("selected genes were not present in the nonzero DESeq2 fit")
}

selected <- selected[match(genes, selected$gene), , drop = FALSE]
row.names(selected) <- NULL
gene_idx <- match(genes, rownames(dds))
disp_init <- ifelse(
  mcols(dds)$dispGeneEst[gene_idx] > 0.1 * mcols(dds)$dispFit[gene_idx],
  mcols(dds)$dispGeneEst[gene_idx],
  mcols(dds)$dispFit[gene_idx]
)
map_mu <- assays(dds_nz)[["mu"]][genes, , drop = FALSE]
map_fit <- DESeq2:::fitDispWrapper(
  ySEXP = counts(dds_nz)[genes, , drop = FALSE],
  xSEXP = design_info$model_matrix,
  mu_hatSEXP = map_mu,
  log_alphaSEXP = log(disp_init),
  log_alpha_prior_meanSEXP = log(mcols(dds)$dispFit[gene_idx]),
  log_alpha_prior_sigmasqSEXP = disp_prior_var,
  min_log_alphaSEXP = log(1e-8 / 10),
  kappa_0SEXP = 1,
  tolSEXP = 1e-6,
  maxitSEXP = 100,
  usePriorSEXP = TRUE,
  weightsSEXP = matrix(1, nrow = length(genes), ncol = ncol(dds_nz)),
  useWeightsSEXP = FALSE,
  weightThresholdSEXP = 1e-2,
  useCRSEXP = TRUE
)

write_tsv(
  data.frame(
    key = c(
      "contrast",
      "design_kind",
      "n_selected",
      "min_se_abs",
      "minmu",
      "disp_prior_var",
      "var_log_disp_estimates"
    ),
    value = c(
      contrast,
      design_info$kind,
      length(genes),
      min_se_abs,
      minmu,
      disp_prior_var,
      var_log_disp_estimates
    )
  ),
  file.path(out_dir, "manifest.tsv")
)
write_tsv(selected, file.path(out_dir, "selected_diagnostics.tsv"))
write_tsv(
  data.frame(
    sample = rownames(design_info$model_matrix),
    design_info$model_matrix,
    check.names = FALSE
  ),
  file.path(out_dir, "design_matrix.tsv")
)
write_tsv(
  data.frame(sample = names(sizeFactors(dds)), sizeFactor = as.numeric(sizeFactors(dds))),
  file.path(out_dir, "size_factors.tsv")
)
write_matrix_tsv(counts(dds)[genes, , drop = FALSE], file.path(out_dir, "selected_counts.tsv"))
write_matrix_tsv(
  map_mu,
  file.path(out_dir, "map_input_mu.tsv")
)
write_tsv(
  data.frame(
    gene = genes,
    dispersion = mcols(dds)$dispersion[gene_idx],
    dispMAP = mcols(dds)$dispMAP[gene_idx],
    dispFit = mcols(dds)$dispFit[gene_idx],
    dispGeneEst = mcols(dds)$dispGeneEst[gene_idx],
    dispInit = disp_init,
    dispIter = mcols(dds)$dispIter[gene_idx],
    dispOutlier = mcols(dds)$dispOutlier[gene_idx]
  ),
  file.path(out_dir, "selected_dispersions.tsv")
)
write_tsv(
  data.frame(
    gene = genes,
    logAlpha = map_fit$log_alpha,
    iter = map_fit$iter,
    iterAccept = map_fit$iter_accept,
    lastChange = map_fit$last_change,
    initialLp = map_fit$initial_lp,
    initialDlp = map_fit$initial_dlp,
    lastLp = map_fit$last_lp,
    lastDlp = map_fit$last_dlp,
    lastD2lp = map_fit$last_d2lp
  ),
  file.path(out_dir, "map_line_search.tsv")
)
write_matrix_tsv(fit$betaMatrix[genes, , drop = FALSE], file.path(out_dir, "beta_no_optim.tsv"))
write_matrix_tsv(fit$betaSE[genes, , drop = FALSE], file.path(out_dir, "beta_se_no_optim.tsv"))
write_matrix_tsv(fit$mu[genes, , drop = FALSE], file.path(out_dir, "mu_no_optim.tsv"))
write_matrix_tsv(fit$hat_diagonals[genes, , drop = FALSE], file.path(out_dir, "hat_no_optim.tsv"))

message("wrote ", length(genes), " SE/covariance rows to ", out_dir)
