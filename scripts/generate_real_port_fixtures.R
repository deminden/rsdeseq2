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

study_root <- get_arg("--study-root", "/home/den/bio/rsfgsea/results/decor_method_study")
out_dir <- get_arg("--output-dir", "results/fixtures/real_optimizer_locfit_2026-05-31")
tissue <- get_arg("--tissue", "kidney")
null_type <- get_arg("--null-type", "full_blocked_permutation")
rep <- as.integer(get_arg("--rep", "1"))
fit_type <- get_arg("--fit-type", "local")
force_top_n <- as.integer(get_arg("--force-top-n", "256"))
min_disp <- as.numeric(get_arg("--min-disp", "1e-8"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "lbfgsb"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "locfit"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
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

infer_design <- function(groups) {
  blocks <- sort(unique(groups$perm_block[!is.na(groups$perm_block) & groups$perm_block != "NA"]))
  full <- stats::model.matrix(~ perm_block + condition, groups)
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

format_num <- function(x) {
  ifelse(is.na(x), NA_character_, format(x, digits = 17, scientific = TRUE))
}

contrast_stem <- sprintf("%s_%s_rep%02d", tissue, null_type, rep)
counts_path <- file.path(study_root, "01_inputs", sprintf("%s_raw_counts.tsv.gz", tissue))
groups_path <- file.path(study_root, "02_null_splits", sprintf("%s_groups.tsv", contrast_stem))
reference_results_path <- file.path(
  study_root,
  "02_deseq2_outputs",
  sprintf("%s_deseq2_results.tsv.gz", contrast_stem)
)

message("reading real-data split")
counts_all <- read_counts(counts_path)
groups <- utils::read.delim(groups_path, check.names = FALSE)
groups <- groups[groups$retained == "TRUE" & groups$condition %in% c("A", "B"), , drop = FALSE]
groups$condition <- factor(groups$condition, levels = c("A", "B"))
groups$perm_block <- factor(groups$perm_block)
counts <- counts_all[, groups$sample_id, drop = FALSE]
storage.mode(counts) <- "integer"

design_info <- infer_design(groups)
col_data <- data.frame(
  row.names = groups$sample_id,
  condition = groups$condition,
  perm_block = groups$perm_block,
  sex = groups$sex,
  age_bin = groups$age_bin,
  stringsAsFactors = FALSE
)

message("building DESeq2 object and estimating local dispersion trend")
dds <- DESeqDataSetFromMatrix(counts, col_data, design = design_info$design)
dds <- estimateSizeFactors(dds)
dds_gene <- estimateDispersionsGeneEst(dds, quiet = TRUE)
dds_fit <- estimateDispersionsFit(dds_gene, fitType = fit_type, quiet = TRUE)
disp_fun <- dispersionFunction(dds_fit)

row_meta_gene <- as.data.frame(S4Vectors::mcols(dds_gene), optional = TRUE)
row_meta_fit <- as.data.frame(S4Vectors::mcols(dds_fit), optional = TRUE)
base_mean <- row_meta_gene$baseMean
disp_gene_est <- row_meta_gene$dispGeneEst
disp_fit <- row_meta_fit$dispFit
locfit_use_for_fit <- is.finite(base_mean) & is.finite(disp_gene_est) &
  base_mean > 0 & disp_gene_est > 100 * min_disp
direct_local_fit <- DESeq2:::localDispersionFit(base_mean, disp_gene_est, minDisp = min_disp)
direct_local_pred <- rep(NA_real_, length(base_mean))
direct_eval <- is.finite(base_mean) & base_mean > 0
direct_local_pred[direct_eval] <- direct_local_fit(base_mean[direct_eval])

locfit_rows <- data.frame(
  gene = rownames(counts),
  baseMean = base_mean,
  dispGeneEst = disp_gene_est,
  useForFit = locfit_use_for_fit,
  dispFit = disp_fit,
  directLocalDispFit = direct_local_pred,
  abs_direct_delta = abs(disp_fit - direct_local_pred),
  log_baseMean = log(base_mean),
  log_dispGeneEst = log(disp_gene_est),
  log_dispFit = log(disp_fit),
  check.names = FALSE
)
write_tsv(locfit_rows, gzfile(file.path(out_dir, "locfit", "local_dispersion_all_rows.tsv.gz")))

fit_rows <- locfit_rows[locfit_rows$useForFit & is.finite(locfit_rows$log_baseMean) &
  is.finite(locfit_rows$log_dispGeneEst), , drop = FALSE]
write_tsv(
  fit_rows[, c("gene", "baseMean", "dispGeneEst", "log_baseMean", "log_dispGeneEst")],
  file.path(out_dir, "locfit", "local_dispersion_fit_points.tsv")
)
write_tsv(
  locfit_rows[order(-abs(locfit_rows$log_dispGeneEst - locfit_rows$log_dispFit)), ][seq_len(min(512, nrow(locfit_rows))), ],
  file.path(out_dir, "locfit", "local_dispersion_ranked_hard_rows.tsv")
)

message("estimating MAP dispersions and GLM optimizer surfaces")
dds_map <- estimateDispersionsMAP(dds_fit, quiet = TRUE)
row_meta_map <- as.data.frame(S4Vectors::mcols(dds_map), optional = TRUE)
glm_rows <- is.finite(row_meta_map$dispersion) & row_meta_map$dispersion > 0
dds_glm <- dds_map[glm_rows, ]
counts_glm <- counts[glm_rows, , drop = FALSE]
model_matrix <- design_info$model_matrix
coef_names <- colnames(model_matrix)
condition_coef <- grep("condition", coef_names)
if (length(condition_coef) == 0) {
  condition_coef <- ncol(model_matrix)
}
condition_coef <- condition_coef[length(condition_coef)]

fit_no_optim <- suppressWarnings(DESeq2:::fitNbinomGLMs(
  dds_glm,
  modelMatrix = model_matrix,
  useOptim = FALSE,
  useQR = FALSE
))
fit_optim <- suppressWarnings(DESeq2:::fitNbinomGLMs(
  dds_glm,
  modelMatrix = model_matrix,
  useOptim = TRUE,
  forceOptim = FALSE,
  useQR = FALSE
))
fit_force <- suppressWarnings(DESeq2:::fitNbinomGLMs(
  dds_glm,
  modelMatrix = model_matrix,
  useOptim = TRUE,
  forceOptim = TRUE,
  useQR = FALSE
))

beta_no <- fit_no_optim$betaMatrix
beta_opt <- fit_optim$betaMatrix
beta_force <- fit_force$betaMatrix
colnames(beta_no) <- coef_names
colnames(beta_opt) <- coef_names
colnames(beta_force) <- coef_names
rownames(beta_no) <- rownames(counts_glm)
rownames(beta_opt) <- rownames(counts_glm)
rownames(beta_force) <- rownames(counts_glm)
rownames(fit_no_optim$betaSE) <- rownames(counts_glm)
rownames(fit_optim$betaSE) <- rownames(counts_glm)
rownames(fit_force$betaSE) <- rownames(counts_glm)
rownames(fit_force$mu) <- rownames(counts_glm)
rownames(fit_force$hat_diagonals) <- rownames(counts_glm)

max_abs_by_row <- function(a, b) {
  apply(abs(a - b), 1, function(row) {
    if (all(is.na(row))) {
      NA_real_
    } else {
      max(row, na.rm = TRUE)
    }
  })
}

actual_delta <- max_abs_by_row(beta_no, beta_opt)
force_delta <- max_abs_by_row(beta_no, beta_force)
actual_optimizer <- is.finite(actual_delta) & actual_delta > 1e-10
rough_or_failed <- (!fit_no_optim$betaConv) | fit_no_optim$betaIter >= 100

all_summary <- data.frame(
  gene = rownames(counts_glm),
  enteredGlm = TRUE,
  baseMean = row_meta_map$baseMean[glm_rows],
  dispGeneEst = row_meta_map$dispGeneEst[glm_rows],
  dispFit = row_meta_map$dispFit[glm_rows],
  dispMAP = row_meta_map$dispMAP[glm_rows],
  dispersion = row_meta_map$dispersion[glm_rows],
  dispOutlier = row_meta_map$dispOutlier[glm_rows],
  noOptimConv = fit_no_optim$betaConv,
  optimConv = fit_optim$betaConv,
  forceConv = fit_force$betaConv,
  noOptimIter = fit_no_optim$betaIter,
  optimIter = fit_optim$betaIter,
  forceIter = fit_force$betaIter,
  actualOptimRouted = actual_optimizer,
  roughOrFailedNoOptim = rough_or_failed,
  maxAbsBetaNoOptimVsOptim = actual_delta,
  maxAbsBetaNoOptimVsForce = force_delta,
  conditionBetaNoOptim = beta_no[, condition_coef],
  conditionBetaOptim = beta_opt[, condition_coef],
  conditionBetaForce = beta_force[, condition_coef],
  conditionSeNoOptim = fit_no_optim$betaSE[, condition_coef],
  conditionSeOptim = fit_optim$betaSE[, condition_coef],
  conditionSeForce = fit_force$betaSE[, condition_coef],
  logLikeNoOptim = fit_no_optim$logLike,
  logLikeOptim = fit_optim$logLike,
  logLikeForce = fit_force$logLike,
  check.names = FALSE
)
write_tsv(all_summary, gzfile(file.path(out_dir, "lbfgsb", "all_gene_optimizer_summary.tsv.gz")))

actual_genes <- all_summary$gene[actual_optimizer | rough_or_failed]
force_order <- order(-force_delta)
force_genes <- all_summary$gene[force_order[seq_len(min(force_top_n, length(force_order)))]]
selected_genes <- unique(c(actual_genes, force_genes))
selected_genes <- selected_genes[!is.na(selected_genes)]

case_summary <- all_summary[match(selected_genes, all_summary$gene), , drop = FALSE]
case_summary$caseKind <- ifelse(case_summary$actualOptimRouted | case_summary$roughOrFailedNoOptim,
  "actual_or_rough_optimizer_row",
  "force_optimizer_probe"
)
write_tsv(case_summary, file.path(out_dir, "lbfgsb", "selected_gene_cases.tsv"))

write_matrix_tsv(counts_glm[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "selected_counts.tsv"))
write_matrix_tsv(model_matrix, file.path(out_dir, "lbfgsb", "design_matrix.tsv"), row_name = "sample")
write_tsv(
  data.frame(sample = names(sizeFactors(dds)), sizeFactor = sizeFactors(dds), check.names = FALSE),
  file.path(out_dir, "lbfgsb", "size_factors.tsv")
)
write_tsv(
  data.frame(sample = rownames(col_data), col_data, check.names = FALSE),
  file.path(out_dir, "lbfgsb", "sample_metadata.tsv")
)
write_matrix_tsv(beta_no[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_no_optim.tsv"))
write_matrix_tsv(beta_opt[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_use_optim.tsv"))
write_matrix_tsv(beta_force[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_force_optim.tsv"))
write_matrix_tsv(fit_no_optim$betaSE[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_se_no_optim.tsv"))
write_matrix_tsv(fit_optim$betaSE[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_se_use_optim.tsv"))
write_matrix_tsv(fit_force$betaSE[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "beta_se_force_optim.tsv"))
write_matrix_tsv(fit_force$mu[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "mu_force_optim.tsv"))
write_matrix_tsv(fit_force$hat_diagonals[selected_genes, , drop = FALSE], file.path(out_dir, "lbfgsb", "hat_force_optim.tsv"))

coef_rows <- do.call(rbind, lapply(seq_along(coef_names), function(i) {
  data.frame(
    gene = selected_genes,
    coefficient_index_1based = i,
    coefficient = coef_names[[i]],
    betaNoOptim = beta_no[selected_genes, i],
    betaUseOptim = beta_opt[selected_genes, i],
    betaForceOptim = beta_force[selected_genes, i],
    betaSeNoOptim = fit_no_optim$betaSE[selected_genes, i],
    betaSeUseOptim = fit_optim$betaSE[selected_genes, i],
    betaSeForceOptim = fit_force$betaSE[selected_genes, i],
    check.names = FALSE
  )
}))
write_tsv(coef_rows, file.path(out_dir, "lbfgsb", "selected_coefficients_long.tsv"))

reference_results <- utils::read.delim(gzfile(reference_results_path), check.names = FALSE)
write_tsv(reference_results[reference_results$gene %in% selected_genes, , drop = FALSE],
  file.path(out_dir, "lbfgsb", "selected_reference_results.tsv")
)

manifest <- data.frame(
  key = c(
    "study_root",
    "tissue",
    "null_type",
    "rep",
    "contrast",
    "fit_type",
    "design",
    "n_genes",
    "n_samples",
    "n_coefficients",
    "condition_coefficient_index_1based",
    "n_locfit_fit_points",
    "n_actual_or_rough_optimizer_rows",
    "n_force_probe_rows",
    "deseq2_version",
    "locfit_version",
    "r_version"
  ),
  value = c(
    basename(normalizePath(study_root)),
    tissue,
    null_type,
    as.character(rep),
    contrast_stem,
    fit_type,
    design_info$kind,
    as.character(nrow(counts)),
    as.character(ncol(counts)),
    as.character(ncol(model_matrix)),
    as.character(condition_coef),
    as.character(nrow(fit_rows)),
    as.character(sum(case_summary$caseKind == "actual_or_rough_optimizer_row")),
    as.character(sum(case_summary$caseKind == "force_optimizer_probe")),
    as.character(utils::packageVersion("DESeq2")),
    as.character(utils::packageVersion("locfit")),
    R.version.string
  ),
  check.names = FALSE
)
write_tsv(manifest, file.path(out_dir, "manifest.tsv"))

readme <- c(
  "# Real DESeq2 Port Fixtures",
  "",
  sprintf("Source contrast: `%s`.", contrast_stem),
  "",
  "This untracked bundle is intended for standalone pure-Rust ports of two numerically important DESeq2 dependencies:",
  "",
  "- `lbfgsb/`: GLM beta optimization surfaces from `DESeq2:::fitNbinomGLMs` on the real split.",
  "- `locfit/`: local dispersion-trend inputs and DESeq2/locfit outputs for every real gene row.",
  "",
  "The R/DESeq2 run is offline fixture generation only. Runtime computation in the Rust crate remains pure Rust.",
  "",
  "Important files:",
  "",
  "- `manifest.tsv`: source, versions, shape, and row counts.",
  "- `locfit/local_dispersion_all_rows.tsv.gz`: exhaustive local-trend input/output rows.",
  "- `locfit/local_dispersion_fit_points.tsv`: exact fit points used by the local trend.",
  "- `locfit/local_dispersion_ranked_hard_rows.tsv`: rows with largest log-dispersion residuals.",
  "- `lbfgsb/all_gene_optimizer_summary.tsv.gz`: exhaustive per-gene optimizer summary for this split.",
  "- `lbfgsb/selected_gene_cases.tsv`: all actual/rough optimizer rows plus ranked force-optimizer probes.",
  "- `lbfgsb/selected_counts.tsv`, `design_matrix.tsv`, `size_factors.tsv`: portable objective inputs.",
  "- `lbfgsb/beta_*`, `beta_se_*`, `mu_force_optim.tsv`, `hat_force_optim.tsv`: DESeq2 targets.",
  "",
  "The selected optimizer case set is intentionally compact enough to copy into a standalone optimizer repository, while `all_gene_optimizer_summary.tsv.gz` preserves the exhaustive real split scan."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote fixture bundle: ", normalizePath(out_dir))
