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
out_dir <- get_arg("--output-dir", "results/fixtures/locfit_hard_real_2026-06-01")
max_contrasts <- as.integer(get_arg("--max-contrasts", "24"))
jobs <- as.integer(get_arg("--jobs", "1"))
hard_top_n <- as.integer(get_arg("--hard-top-n", "512"))
min_disp <- as.numeric(get_arg("--min-disp", "1e-8"))
explicit_contrasts <- get_args("--contrast")

jobs <- max(1, min(jobs, 3))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  con <- if (is.character(path) && grepl("[.]gz$", path)) gzfile(path, "wt") else path
  if (inherits(con, "connection")) {
    on.exit(close(con), add = TRUE)
  }
  utils::write.table(
    x,
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
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
      break
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
  }
  do.call(rbind, selected)
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

safe_log <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log(x[ok])
  out
}

rank_hard_rows <- function(rows, top_n) {
  finite <- is.finite(rows$baseMean) & is.finite(rows$dispGeneEst) &
    is.finite(rows$dispFit) & rows$baseMean > 0 & rows$dispGeneEst > 0 & rows$dispFit > 0
  residual <- rep(0, nrow(rows))
  residual[finite] <- abs(log(rows$dispGeneEst[finite]) - log(rows$dispFit[finite]))
  direct_delta <- ifelse(is.finite(rows$absDirectDelta), rows$absDirectDelta, 0)
  threshold_distance <- abs(log10(pmax(rows$dispGeneEst / (100 * min_disp), 1e-300)))
  threshold_score <- ifelse(is.finite(threshold_distance), 1 / (1 + threshold_distance), 0)
  tail_score <- rep(0, nrow(rows))
  if (any(is.finite(rows$baseMean) & rows$baseMean > 0)) {
    log_mean <- safe_log(rows$baseMean)
    finite_mean <- is.finite(log_mean)
    q <- stats::quantile(log_mean[finite_mean], c(0.01, 0.99), na.rm = TRUE)
    tail_score[finite_mean & (log_mean <= q[[1]] | log_mean >= q[[2]])] <- 1
  }
  rows$hardScore <- pmax(
    residual,
    log10(1 + direct_delta / 1e-12),
    threshold_score,
    tail_score,
    ifelse(rows$useForFit, 0.1, 0),
    na.rm = TRUE
  )
  rows[order(-rows$hardScore), , drop = FALSE][seq_len(min(top_n, nrow(rows))), , drop = FALSE]
}

prediction_grid <- function(base_mean) {
  positive <- base_mean[is.finite(base_mean) & base_mean > 0]
  if (length(positive) == 0) {
    return(data.frame(gridKind = character(), baseMean = numeric()))
  }
  q_probs <- unique(c(seq(0, 1, by = 0.01), 0.001, 0.005, 0.995, 0.999))
  q_values <- as.numeric(stats::quantile(positive, q_probs, na.rm = TRUE, type = 7))
  log_values <- exp(seq(log(min(positive)), log(max(positive)), length.out = 256))
  edge_values <- c(min(positive) / 100, min(positive) / 10, min(positive), max(positive), max(positive) * 10, max(positive) * 100)
  values <- sort(unique(c(q_values, log_values, edge_values)))
  values <- values[is.finite(values) & values > 0]
  data.frame(
    gridKind = ifelse(values < min(positive) | values > max(positive), "extrapolation", "interpolation"),
    baseMean = values,
    check.names = FALSE
  )
}

predict_safe <- function(fit, values) {
  out <- rep(NA_real_, length(values))
  for (idx in seq_along(values)) {
    out[[idx]] <- tryCatch(fit(values[[idx]]), error = function(error) NA_real_)
  }
  out
}

process_contrast <- function(spec) {
  stem <- contrast_stem(spec)
  message("processing ", stem)
  contrast_dir <- file.path(out_dir, "contrasts", stem)
  locfit_dir <- file.path(contrast_dir, "locfit")
  dir.create(locfit_dir, recursive = TRUE, showWarnings = FALSE)

  counts_path <- file.path(study_root, "01_inputs", sprintf("%s_raw_counts.tsv.gz", spec$tissue))
  groups_path <- file.path(study_root, "02_null_splits", sprintf("%s_groups.tsv", stem))

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
  dds_gene <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds_fit <- estimateDispersionsFit(dds_gene, fitType = "local", quiet = TRUE)

  gene_meta <- as.data.frame(S4Vectors::mcols(dds_gene), optional = TRUE)
  fit_meta <- as.data.frame(S4Vectors::mcols(dds_fit), optional = TRUE)
  base_mean <- gene_meta$baseMean
  disp_gene_est <- gene_meta$dispGeneEst
  disp_fit <- fit_meta$dispFit
  use_for_fit <- is.finite(base_mean) & is.finite(disp_gene_est) &
    base_mean > 0 & disp_gene_est > 100 * min_disp
  direct_fit <- DESeq2:::localDispersionFit(base_mean, disp_gene_est, minDisp = min_disp)
  eval_rows <- is.finite(base_mean) & base_mean > 0
  direct_pred <- rep(NA_real_, length(base_mean))
  direct_pred[eval_rows] <- predict_safe(direct_fit, base_mean[eval_rows])

  all_rows <- data.frame(
    contrast = stem,
    tissue = spec$tissue,
    null_type = spec$null_type,
    rep = spec$rep,
    gene = rownames(counts),
    baseMean = base_mean,
    dispGeneEst = disp_gene_est,
    useForFit = use_for_fit,
    dispFit = disp_fit,
    directLocalDispFit = direct_pred,
    absDirectDelta = abs(disp_fit - direct_pred),
    logBaseMean = safe_log(base_mean),
    logDispGeneEst = safe_log(disp_gene_est),
    logDispFit = safe_log(disp_fit),
    thresholdRatio = disp_gene_est / (100 * min_disp),
    check.names = FALSE
  )
  write_tsv(all_rows, file.path(locfit_dir, "all_rows.tsv.gz"))

  fit_points <- all_rows[all_rows$useForFit & is.finite(all_rows$logBaseMean) & is.finite(all_rows$logDispGeneEst), , drop = FALSE]
  write_tsv(fit_points[, c("gene", "baseMean", "dispGeneEst", "logBaseMean", "logDispGeneEst")], file.path(locfit_dir, "fit_points.tsv"))

  hard_rows <- rank_hard_rows(all_rows, hard_top_n)
  write_tsv(hard_rows, file.path(locfit_dir, "hard_rows.tsv"))

  grid <- prediction_grid(base_mean)
  if (nrow(grid) > 0) {
    grid$directLocalDispFit <- predict_safe(direct_fit, grid$baseMean)
    grid$logBaseMean <- safe_log(grid$baseMean)
    grid$logDirectLocalDispFit <- safe_log(grid$directLocalDispFit)
  }
  write_tsv(grid, file.path(locfit_dir, "prediction_grid.tsv"))

  manifest <- data.frame(
    key = c(
      "contrast", "design", "n_genes", "n_samples", "n_fit_points",
      "n_hard_rows", "n_prediction_grid", "min_fit_base_mean", "max_fit_base_mean",
      "min_fit_disp_gene_est", "max_fit_disp_gene_est"
    ),
    value = c(
      stem, design_info$kind, nrow(counts), ncol(counts), nrow(fit_points),
      nrow(hard_rows), nrow(grid),
      min(fit_points$baseMean, na.rm = TRUE), max(fit_points$baseMean, na.rm = TRUE),
      min(fit_points$dispGeneEst, na.rm = TRUE), max(fit_points$dispGeneEst, na.rm = TRUE)
    ),
    check.names = FALSE
  )
  write_tsv(manifest, file.path(contrast_dir, "manifest.tsv"))

  list(
    summary = data.frame(
      contrast = stem,
      tissue = spec$tissue,
      null_type = spec$null_type,
      rep = spec$rep,
      design = design_info$kind,
      nGenes = nrow(counts),
      nSamples = ncol(counts),
      nFitPoints = nrow(fit_points),
      nHardRows = nrow(hard_rows),
      nPredictionGrid = nrow(grid),
      stringsAsFactors = FALSE
    ),
    hard = hard_rows
  )
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

specs <- lapply(seq_len(nrow(selected_specs)), function(idx) {
  list(tissue = selected_specs$tissue[idx], null_type = selected_specs$null_type[idx], rep = selected_specs$rep[idx])
})

results <- parallel::mclapply(
  specs,
  function(spec) {
    tryCatch(process_contrast(spec), error = function(error) {
      list(
        failure = data.frame(
          contrast = contrast_stem(spec),
          error = conditionMessage(error),
          stringsAsFactors = FALSE
        )
      )
    })
  },
  mc.cores = jobs
)

summaries <- lapply(results, `[[`, "summary")
summaries <- summaries[!vapply(summaries, is.null, logical(1))]
failures <- lapply(results, `[[`, "failure")
failures <- failures[!vapply(failures, is.null, logical(1))]
hard <- lapply(results, `[[`, "hard")
hard <- hard[!vapply(hard, is.null, logical(1))]

summary_rows <- if (length(summaries) > 0) do.call(rbind, summaries) else data.frame()
hard_rows <- if (length(hard) > 0) do.call(rbind, hard) else data.frame()
failure_rows <- if (length(failures) > 0) do.call(rbind, failures) else data.frame()

if (nrow(summary_rows) > 0) {
  write_tsv(summary_rows, file.path(out_dir, "completed_contrasts.tsv"))
}
if (nrow(hard_rows) > 0) {
  hard_rows <- hard_rows[order(-hard_rows$hardScore), , drop = FALSE]
  write_tsv(hard_rows, file.path(out_dir, "global_hard_rows.tsv"))
  write_tsv(hard_rows[seq_len(min(2048, nrow(hard_rows))), , drop = FALSE], file.path(out_dir, "global_hardest_2048.tsv"))
}
if (nrow(failure_rows) > 0) {
  write_tsv(failure_rows, file.path(out_dir, "failures.tsv"))
}

manifest <- data.frame(
  key = c(
    "study_root", "requested_contrasts", "successful_contrasts", "failed_contrasts",
    "jobs", "hard_rows", "deseq2_version", "locfit_version", "r_version"
  ),
  value = c(
    normalizePath(study_root), nrow(selected_specs), nrow(summary_rows), nrow(failure_rows),
    jobs, nrow(hard_rows), as.character(utils::packageVersion("DESeq2")),
    as.character(utils::packageVersion("locfit")), R.version.string
  ),
  check.names = FALSE
)
write_tsv(manifest, file.path(out_dir, "manifest.tsv"))

readme <- c(
  "# Hard Locfit Real-Contrast Fixtures",
  "",
  "This ignored bundle is generated from saved real publication-study contrasts.",
  "It targets the local dispersion trend used by DESeq2's `fitType=\"local\"` path.",
  "",
  "Start with `global_hardest_2048.tsv`, then drill into `contrasts/<contrast>/locfit/`.",
  "",
  "Per-contrast files:",
  "",
  "- `all_rows.tsv.gz`: every gene row with `baseMean`, `dispGeneEst`, DESeq2 `dispFit`, direct `localDispersionFit` prediction, log-space fields, and threshold metadata.",
  "- `fit_points.tsv`: exact positive rows used as local-fit training points.",
  "- `hard_rows.tsv`: ranked local-fit stress rows for this contrast.",
  "- `prediction_grid.tsv`: dense interpolation and extrapolation target grid on the real base-mean range.",
  "",
  "The R/DESeq2 run is offline fixture generation only; runtime computation in the Rust crate remains pure Rust."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote hard locfit fixture bundle: ", normalizePath(out_dir))
