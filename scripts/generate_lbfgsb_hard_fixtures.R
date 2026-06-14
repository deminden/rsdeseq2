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
get_args <- function(name) {
  hit <- which(args == name)
  hit <- hit[hit < length(args)]
  args[hit + 1]
}

study_root <- get_arg("--study-root", Sys.getenv("RSDESEQ2_REAL_DATA_ROOT", unset = NA_character_))
if (is.na(study_root) || !nzchar(study_root)) {
  stop("provide --study-root or set RSDESEQ2_REAL_DATA_ROOT", call. = FALSE)
}
out_dir <- get_arg("--output-dir", "results/fixtures/lbfgsb_hard_real_2026-06-01")
fit_type <- get_arg("--fit-type", "parametric")
max_contrasts <- as.integer(get_arg("--max-contrasts", "16"))
per_contrast_top_n <- as.integer(get_arg("--per-contrast-top-n", "96"))
force_top_n <- as.integer(get_arg("--force-top-n", "192"))
min_beta_delta <- as.numeric(get_arg("--min-beta-delta", "1e-10"))
explicit_contrasts <- get_args("--contrast")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

parse_contrast <- function(spec) {
  parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
  if (length(parts) != 3) {
    stop("contrast must be tissue:null_type:rep, got: ", spec)
  }
  list(tissue = parts[[1]], null_type = parts[[2]], rep = as.integer(parts[[3]]))
}

contrast_stem <- function(spec) {
  sprintf("%s_%s_rep%02d", spec$tissue, spec$null_type, spec$rep)
}

list_available_contrasts <- function(root) {
  output_dir <- file.path(root, "02_deseq2_outputs")
  group_dir <- file.path(root, "02_null_splits")
  files <- list.files(output_dir, pattern = "_deseq2_results[.]tsv[.]gz$", full.names = FALSE)
  null_types <- c("full_blocked_permutation", "blocked_permutation", "stratified_split", "random_split")
  rows <- list()
  for (file in files) {
    stem <- sub("_deseq2_results[.]tsv[.]gz$", "", file)
    matched <- FALSE
    for (null_type in null_types) {
      marker <- paste0("_", null_type, "_rep")
      if (!grepl(marker, stem, fixed = TRUE)) {
        next
      }
      parts <- strsplit(stem, marker, fixed = TRUE)[[1]]
      tissue <- parts[[1]]
      rep <- as.integer(parts[[2]])
      groups <- file.path(group_dir, sprintf("%s_%s_rep%02d_groups.tsv", tissue, null_type, rep))
      counts <- file.path(root, "01_inputs", sprintf("%s_raw_counts.tsv.gz", tissue))
      if (file.exists(groups) && file.exists(counts)) {
        rows[[length(rows) + 1]] <- data.frame(
          tissue = tissue,
          null_type = null_type,
          rep = rep,
          stringsAsFactors = FALSE
        )
      }
      matched <- TRUE
      break
    }
    if (!matched) {
      warning("could not parse saved contrast output: ", file)
    }
  }
  do.call(rbind, rows)
}

default_contrasts <- function(root, max_n) {
  available <- list_available_contrasts(root)
  tissue_priority <- c(
    "kidney", "liver", "pancreas", "heart", "muscle", "blood", "lung", "brain",
    "thyroid", "testis", "adipose_subcutaneous", "artery_tibial", "spleen",
    "colon_transverse", "esophagus_mucosa", "breast_mammary", "nerve_tibial"
  )
  null_priority <- c("full_blocked_permutation", "blocked_permutation")
  rep_priority <- c(1, 5, 10, 15, 20, 2, 3, 8, 13, 17)
  selected <- list()
  for (tissue in tissue_priority) {
    for (null_type in null_priority) {
      for (rep in rep_priority) {
        hit <- available[
          available$tissue == tissue & available$null_type == null_type & available$rep == rep,
          ,
          drop = FALSE
        ]
        if (nrow(hit) > 0) {
          selected[[length(selected) + 1]] <- hit[1, , drop = FALSE]
          break
        }
      }
      if (length(selected) >= max_n) {
        return(do.call(rbind, selected))
      }
    }
    if (length(selected) >= max_n) {
      return(do.call(rbind, selected))
    }
  }
  selected_frame <- if (length(selected) > 0) do.call(rbind, selected) else available[FALSE, ]
  if (nrow(selected_frame) < max_n) {
    key <- paste(selected_frame$tissue, selected_frame$null_type, selected_frame$rep)
    for (idx in seq_len(nrow(available))) {
      candidate_key <- paste(available$tissue[idx], available$null_type[idx], available$rep[idx])
      if (candidate_key %in% key) {
        next
      }
      selected_frame <- rbind(selected_frame, available[idx, , drop = FALSE])
      key <- c(key, candidate_key)
      if (nrow(selected_frame) >= max_n) {
        break
      }
    }
  }
  selected_frame
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

max_abs_by_row <- function(a, b) {
  apply(abs(a - b), 1, function(row) {
    if (all(is.na(row))) {
      NA_real_
    } else {
      max(row, na.rm = TRUE)
    }
  })
}

safe_col <- function(frame, name, default = NA_real_) {
  if (name %in% names(frame)) {
    frame[[name]]
  } else {
    rep(default, nrow(frame))
  }
}

process_contrast <- function(spec) {
  stem <- contrast_stem(spec)
  message("processing ", stem)
  contrast_dir <- file.path(out_dir, "contrasts", stem)
  dir.create(file.path(contrast_dir, "lbfgsb"), recursive = TRUE, showWarnings = FALSE)

  counts_path <- file.path(study_root, "01_inputs", sprintf("%s_raw_counts.tsv.gz", spec$tissue))
  groups_path <- file.path(study_root, "02_null_splits", sprintf("%s_groups.tsv", stem))
  reference_results_path <- file.path(study_root, "02_deseq2_outputs", sprintf("%s_deseq2_results.tsv.gz", stem))

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
    sex = groups$sex,
    age_bin = groups$age_bin,
    stringsAsFactors = FALSE
  )

  dds <- DESeqDataSetFromMatrix(counts, col_data, design = design_info$design)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, fitType = fit_type, quiet = TRUE)
  dds <- estimateDispersionsMAP(dds, quiet = TRUE)

  row_meta <- as.data.frame(S4Vectors::mcols(dds), optional = TRUE)
  glm_rows <- is.finite(row_meta$dispersion) & row_meta$dispersion > 0
  dds_glm <- dds[glm_rows, ]
  counts_glm <- counts[glm_rows, , drop = FALSE]
  model_matrix <- design_info$model_matrix
  coef_names <- colnames(model_matrix)
  condition_coef <- grep("condition", coef_names)
  if (length(condition_coef) == 0) {
    condition_coef <- ncol(model_matrix)
  }
  condition_coef <- condition_coef[length(condition_coef)]

  fit_no <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    dds_glm,
    modelMatrix = model_matrix,
    useOptim = FALSE,
    useQR = FALSE
  ))
  fit_use <- suppressWarnings(DESeq2:::fitNbinomGLMs(
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

  beta_no <- fit_no$betaMatrix
  beta_use <- fit_use$betaMatrix
  beta_force <- fit_force$betaMatrix
  colnames(beta_no) <- coef_names
  colnames(beta_use) <- coef_names
  colnames(beta_force) <- coef_names
  rownames(beta_no) <- rownames(counts_glm)
  rownames(beta_use) <- rownames(counts_glm)
  rownames(beta_force) <- rownames(counts_glm)
  rownames(fit_no$betaSE) <- rownames(counts_glm)
  rownames(fit_use$betaSE) <- rownames(counts_glm)
  rownames(fit_force$betaSE) <- rownames(counts_glm)
  rownames(fit_force$mu) <- rownames(counts_glm)
  rownames(fit_force$hat_diagonals) <- rownames(counts_glm)

  actual_delta <- max_abs_by_row(beta_no, beta_use)
  force_delta <- max_abs_by_row(beta_no, beta_force)
  actual_optimizer <- is.finite(actual_delta) & actual_delta > min_beta_delta
  rough_or_failed <- (!fit_no$betaConv) | fit_no$betaIter >= 100
  hard_score <- pmax(
    ifelse(is.finite(force_delta), log10(1 + force_delta / 1e-12), 0),
    ifelse(is.finite(actual_delta), log10(1 + actual_delta / 1e-12), 0),
    fit_no$betaIter / 100,
    ifelse(!fit_no$betaConv, 20, 0),
    na.rm = TRUE
  )

  all_summary <- data.frame(
    contrast = stem,
    tissue = spec$tissue,
    null_type = spec$null_type,
    rep = spec$rep,
    gene = rownames(counts_glm),
    baseMean = row_meta$baseMean[glm_rows],
    dispGeneEst = safe_col(row_meta, "dispGeneEst")[glm_rows],
    dispFit = safe_col(row_meta, "dispFit")[glm_rows],
    dispMAP = safe_col(row_meta, "dispMAP")[glm_rows],
    dispersion = row_meta$dispersion[glm_rows],
    dispOutlier = safe_col(row_meta, "dispOutlier", FALSE)[glm_rows],
    noOptimConv = fit_no$betaConv,
    useOptimConv = fit_use$betaConv,
    forceOptimConv = fit_force$betaConv,
    noOptimIter = fit_no$betaIter,
    useOptimIter = fit_use$betaIter,
    forceOptimIter = fit_force$betaIter,
    actualOptimRouted = actual_optimizer,
    roughOrFailedNoOptim = rough_or_failed,
    maxAbsBetaNoOptimVsUseOptim = actual_delta,
    maxAbsBetaNoOptimVsForceOptim = force_delta,
    hardScore = hard_score,
    conditionBetaNoOptim = beta_no[, condition_coef],
    conditionBetaUseOptim = beta_use[, condition_coef],
    conditionBetaForceOptim = beta_force[, condition_coef],
    conditionSeNoOptim = fit_no$betaSE[, condition_coef],
    conditionSeUseOptim = fit_use$betaSE[, condition_coef],
    conditionSeForceOptim = fit_force$betaSE[, condition_coef],
    logLikeNoOptim = fit_no$logLike,
    logLikeUseOptim = fit_use$logLike,
    logLikeForceOptim = fit_force$logLike,
    check.names = FALSE
  )
  write_tsv(all_summary, gzfile(file.path(contrast_dir, "lbfgsb", "all_gene_optimizer_summary.tsv.gz")))

  must_keep <- all_summary$gene[actual_optimizer | rough_or_failed]
  force_order <- order(-force_delta)
  force_keep <- all_summary$gene[force_order[seq_len(min(force_top_n, length(force_order)))]]
  score_order <- order(-hard_score)
  score_keep <- all_summary$gene[score_order[seq_len(min(per_contrast_top_n, length(score_order)))]]
  selected_genes <- unique(c(must_keep, force_keep, score_keep))
  selected_genes <- selected_genes[!is.na(selected_genes)]
  selected_rows <- match(selected_genes, rownames(counts_glm))

  case_summary <- all_summary[match(selected_genes, all_summary$gene), , drop = FALSE]
  case_summary$caseKind <- ifelse(
    case_summary$actualOptimRouted | case_summary$roughOrFailedNoOptim,
    "actual_or_rough_optimizer_row",
    "force_optimizer_probe"
  )
  write_tsv(case_summary, file.path(contrast_dir, "lbfgsb", "selected_gene_cases.tsv"))

  write_matrix_tsv(counts_glm[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "selected_counts.tsv"))
  write_matrix_tsv(model_matrix, file.path(contrast_dir, "lbfgsb", "design_matrix.tsv"), row_name = "sample")
  write_tsv(
    data.frame(sample = names(sizeFactors(dds)), sizeFactor = sizeFactors(dds), check.names = FALSE),
    file.path(contrast_dir, "lbfgsb", "size_factors.tsv")
  )
  write_tsv(
    data.frame(sample = rownames(col_data), col_data, check.names = FALSE),
    file.path(contrast_dir, "lbfgsb", "sample_metadata.tsv")
  )

  selected_dispersion <- data.frame(
    gene = selected_genes,
    dispersion = row_meta$dispersion[glm_rows][selected_rows],
    dispMAP = safe_col(row_meta, "dispMAP")[glm_rows][selected_rows],
    dispFit = safe_col(row_meta, "dispFit")[glm_rows][selected_rows],
    check.names = FALSE
  )
  write_tsv(selected_dispersion, file.path(contrast_dir, "lbfgsb", "selected_dispersions.tsv"))
  write_matrix_tsv(beta_no[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_no_optim.tsv"))
  write_matrix_tsv(beta_use[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_use_optim.tsv"))
  write_matrix_tsv(beta_force[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_force_optim.tsv"))
  write_matrix_tsv(fit_no$betaSE[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_se_no_optim.tsv"))
  write_matrix_tsv(fit_use$betaSE[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_se_use_optim.tsv"))
  write_matrix_tsv(fit_force$betaSE[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "beta_se_force_optim.tsv"))
  write_matrix_tsv(fit_force$mu[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "mu_force_optim.tsv"))
  write_matrix_tsv(fit_force$hat_diagonals[selected_genes, , drop = FALSE], file.path(contrast_dir, "lbfgsb", "hat_force_optim.tsv"))

  design_long <- data.frame(sample = rownames(model_matrix), model_matrix, check.names = FALSE)
  input_long <- do.call(rbind, lapply(selected_genes, function(gene) {
    data.frame(
      gene = gene,
      sample = colnames(counts_glm),
      count = as.integer(counts_glm[gene, ]),
      sizeFactor = as.numeric(sizeFactors(dds)[colnames(counts_glm)]),
      dispersion = row_meta$dispersion[glm_rows][match(gene, rownames(counts_glm))],
      weight = 1,
      design_long[match(colnames(counts_glm), design_long$sample), -1, drop = FALSE],
      check.names = FALSE
    )
  }))
  write_tsv(input_long, gzfile(file.path(contrast_dir, "lbfgsb", "selected_optimizer_inputs_long.tsv.gz")))

  coef_rows <- do.call(rbind, lapply(seq_along(coef_names), function(i) {
    data.frame(
      contrast = stem,
      gene = selected_genes,
      coefficient_index_1based = i,
      coefficient = coef_names[[i]],
      betaNoOptim = beta_no[selected_genes, i],
      betaUseOptim = beta_use[selected_genes, i],
      betaForceOptim = beta_force[selected_genes, i],
      betaSeNoOptim = fit_no$betaSE[selected_genes, i],
      betaSeUseOptim = fit_use$betaSE[selected_genes, i],
      betaSeForceOptim = fit_force$betaSE[selected_genes, i],
      check.names = FALSE
    )
  }))
  write_tsv(coef_rows, file.path(contrast_dir, "lbfgsb", "selected_coefficients_long.tsv"))

  reference_results <- utils::read.delim(gzfile(reference_results_path), check.names = FALSE)
  write_tsv(reference_results[reference_results$gene %in% selected_genes, , drop = FALSE],
    file.path(contrast_dir, "lbfgsb", "selected_reference_results.tsv")
  )

  manifest <- data.frame(
    key = c(
      "contrast", "fit_type", "design", "n_genes", "n_glm_rows", "n_samples",
      "n_coefficients", "condition_coefficient_index_1based", "n_selected_genes",
      "n_actual_or_rough_optimizer_rows", "n_force_probe_rows"
    ),
    value = c(
      stem, fit_type, design_info$kind, nrow(counts), nrow(counts_glm), ncol(counts),
      ncol(model_matrix), condition_coef, length(selected_genes),
      sum(case_summary$caseKind == "actual_or_rough_optimizer_row"),
      sum(case_summary$caseKind == "force_optimizer_probe")
    ),
    check.names = FALSE
  )
  write_tsv(manifest, file.path(contrast_dir, "manifest.tsv"))

  case_summary
}

selected_specs <- if (length(explicit_contrasts) > 0) {
  do.call(rbind, lapply(explicit_contrasts, function(spec) {
    parsed <- parse_contrast(spec)
    data.frame(tissue = parsed$tissue, null_type = parsed$null_type, rep = parsed$rep, stringsAsFactors = FALSE)
  }))
} else {
  default_contrasts(study_root, max_contrasts)
}

write_tsv(selected_specs, file.path(out_dir, "selected_contrasts.tsv"))

all_cases <- list()
failures <- list()
for (idx in seq_len(nrow(selected_specs))) {
  spec <- list(
    tissue = selected_specs$tissue[idx],
    null_type = selected_specs$null_type[idx],
    rep = selected_specs$rep[idx]
  )
  result <- tryCatch(
    process_contrast(spec),
    error = function(error) {
      message("failed ", contrast_stem(spec), ": ", conditionMessage(error))
      failures[[length(failures) + 1]] <<- data.frame(
        contrast = contrast_stem(spec),
        error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )
  if (!is.null(result)) {
    all_cases[[length(all_cases) + 1]] <- result
  }
}

global_cases <- if (length(all_cases) > 0) do.call(rbind, all_cases) else data.frame()
if (nrow(global_cases) > 0) {
  global_cases <- global_cases[order(-global_cases$hardScore), , drop = FALSE]
  write_tsv(global_cases, file.path(out_dir, "global_selected_gene_cases.tsv"))
  write_tsv(global_cases[seq_len(min(512, nrow(global_cases))), , drop = FALSE], file.path(out_dir, "global_hardest_512.tsv"))
}
if (length(failures) > 0) {
  write_tsv(do.call(rbind, failures), file.path(out_dir, "failures.tsv"))
}

bundle_manifest <- data.frame(
  key = c(
    "study_root", "fit_type", "requested_contrasts", "successful_contrasts",
    "selected_gene_cases", "deseq2_version", "r_version"
  ),
  value = c(
    normalizePath(study_root), fit_type, nrow(selected_specs), length(all_cases),
    nrow(global_cases), as.character(utils::packageVersion("DESeq2")), R.version.string
  ),
  check.names = FALSE
)
write_tsv(bundle_manifest, file.path(out_dir, "manifest.tsv"))

readme <- c(
  "# Hard L-BFGS-B Real-Contrast Fixtures",
  "",
  "This ignored bundle is generated from saved real publication-study contrasts.",
  "It is intended for improving a standalone pure-Rust implementation of R-style L-BFGS-B behavior used by DESeq2's GLM beta fallback.",
  "",
  "Use `global_hardest_512.tsv` first. It ranks selected rows across all processed contrasts by actual optimizer routing, rough/non-converged IRLS starts, forced-optimizer beta movement, and iteration pressure.",
  "",
  "Per-contrast fixture folders live under `contrasts/<contrast>/lbfgsb/`.",
  "",
  "Most useful files per contrast:",
  "",
  "- `selected_optimizer_inputs_long.tsv.gz`: portable objective inputs for selected genes: counts, sample rows, size factors, dispersion, weight, and design columns.",
  "- `selected_counts.tsv`, `design_matrix.tsv`, `size_factors.tsv`, `selected_dispersions.tsv`: matrix-shaped inputs.",
  "- `beta_no_optim.tsv`: IRLS start surface.",
  "- `beta_use_optim.tsv`: DESeq2's normal `useOptim=TRUE` target.",
  "- `beta_force_optim.tsv`: exhaustive forced L-BFGS-B target.",
  "- `selected_coefficients_long.tsv`: coefficient-wise no/use/force beta and SE targets.",
  "- `mu_force_optim.tsv`, `hat_force_optim.tsv`: fitted mean and hat targets from forced optimizer fits.",
  "- `selected_reference_results.tsv`: saved high-level result rows for context.",
  "- `all_gene_optimizer_summary.tsv.gz`: exhaustive per-contrast scan summary.",
  "",
  "The R/DESeq2 run is offline fixture generation only; runtime computation in the Rust crate remains pure Rust."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote hard fixture bundle: ", normalizePath(out_dir))
