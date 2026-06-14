#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- which(args == name)
  if (length(hit) == 0 || hit[length(hit)] == length(args)) {
    return(default)
  }
  args[hit[length(hit)] + 1]
}

source_root <- get_arg("--source-root", "results/fixtures/lbfgsb_hard_real_2026-06-01")
out_dir <- get_arg("--output-dir", "results/fixtures/lbfgsb_synthetic_stress_2026-06-06")
n_cases <- as.integer(get_arg("--n-cases", "512"))
n_seeds <- as.integer(get_arg("--n-seeds", "128"))
max_per_source <- as.integer(get_arg("--max-per-source", "16"))
lower <- as.numeric(get_arg("--lower", "-30"))
upper <- as.numeric(get_arg("--upper", "30"))
maxit <- as.integer(get_arg("--maxit", "100"))

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

bind_rows_fill <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  filled <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) {
      row[[name]] <- NA
    }
    row[, columns, drop = FALSE]
  })
  do.call(rbind, filled)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE)
}

read_matrix <- function(path) {
  frame <- read_tsv(path)
  row_names <- frame[[1]]
  values <- as.matrix(frame[-1])
  storage.mode(values) <- "double"
  rownames(values) <- row_names
  values
}

read_counts <- function(path) {
  frame <- read_tsv(path)
  row_names <- frame[[1]]
  values <- as.matrix(frame[-1])
  storage.mode(values) <- "integer"
  rownames(values) <- row_names
  values
}

read_named_values <- function(path, value_col) {
  frame <- read_tsv(path)
  stats::setNames(as.numeric(frame[[value_col]]), frame[[1]])
}

objective_value <- function(beta, counts, design, size_factors, dispersion, ridge_log2) {
  eta <- as.vector(design %*% beta)
  mu <- size_factors * 2^eta
  if (any(!is.finite(mu) | mu <= 0)) {
    return(1e300)
  }
  size <- 1 / dispersion
  term <- (counts + size) * log(size + mu) - size * log(size) - counts * log(mu)
  value <- sum(term) + 0.5 * sum(ridge_log2 * beta * beta)
  if (!is.finite(value)) 1e300 else value
}

numeric_gradient <- function(beta, fn, lower, upper, ndeps = 1e-3) {
  out <- numeric(length(beta))
  for (i in seq_along(beta)) {
    forward <- min(beta[[i]] + ndeps, upper)
    backward <- max(beta[[i]] - ndeps, lower)
    denom <- forward - backward
    if (!is.finite(denom) || denom <= 0) {
      out[[i]] <- NA_real_
      next
    }
    b_forward <- beta
    b_backward <- beta
    b_forward[[i]] <- forward
    b_backward[[i]] <- backward
    out[[i]] <- (fn(b_forward) - fn(b_backward)) / denom
  }
  out
}

projected_grad_norm <- function(beta, gradient, lower, upper) {
  projected <- gradient
  at_lower <- beta <= lower & gradient > 0
  at_upper <- beta >= upper & gradient < 0
  projected[at_lower | at_upper] <- 0
  max(abs(projected), na.rm = TRUE)
}

stable_case_signature <- function(candidate) {
  paste(
    candidate$contrast,
    candidate$source_gene,
    candidate$count_kind,
    candidate$dispersion_kind,
    candidate$start_kind,
    candidate$total_bin,
    candidate$zero_bin,
    candidate$dispersion_bin,
    candidate$start_bound_count,
    sep = "|"
  )
}

make_count_variants <- function(counts, design) {
  condition_col <- grep("condition", colnames(design))
  condition <- if (length(condition_col) > 0) design[, condition_col[[length(condition_col)]]] else rep(0, length(counts))
  list(
    real = counts,
    low_all = as.integer(round(counts * 0.25)),
    high_all = as.integer(round(counts * 4)),
    very_high_all = as.integer(round(counts * 16)),
    condition_low = as.integer(round(ifelse(condition > 0, counts * 0.1, counts))),
    condition_high = as.integer(round(ifelse(condition > 0, counts * 10, counts)))
  )
}

make_start_variants <- function(beta_no, beta_force) {
  clipped_no <- pmin(pmax(beta_no, lower), upper)
  clipped_force <- pmin(pmax(beta_force, lower), upper)
  delta <- clipped_no - clipped_force
  worst <- which.max(abs(delta))
  lower_edge <- clipped_no
  upper_edge <- clipped_no
  lower_edge[[worst]] <- lower + 1e-9
  upper_edge[[worst]] <- upper - 1e-9
  midpoint <- pmin(pmax((clipped_no + clipped_force) / 2, lower), upper)
  force_jitter <- pmin(pmax(clipped_force + sign(delta) * 1e-3, lower), upper)
  list(
    clipped_no_optim = clipped_no,
    midpoint_to_force = midpoint,
    near_lower_worst_coef = lower_edge,
    near_upper_worst_coef = upper_edge,
    force_target_jitter = force_jitter
  )
}

hardest <- read_tsv(file.path(source_root, "global_hardest_512.tsv"))
hardest <- hardest[order(-hardest$hardScore, -hardest$maxAbsBetaNoOptimVsForceOptim), , drop = FALSE]
hardest <- hardest[seq_len(min(n_seeds, nrow(hardest))), , drop = FALSE]

case_pool <- list()
case_inputs <- list()
case_coefficients <- list()
case_idx <- 0L

loaded <- new.env(parent = emptyenv())
load_contrast <- function(contrast) {
  if (exists(contrast, envir = loaded, inherits = FALSE)) {
    return(get(contrast, envir = loaded, inherits = FALSE))
  }
  root <- file.path(source_root, "contrasts", contrast, "lbfgsb")
  value <- list(
    counts = read_counts(file.path(root, "selected_counts.tsv")),
    design = read_matrix(file.path(root, "design_matrix.tsv")),
    size_factors = read_named_values(file.path(root, "size_factors.tsv"), "sizeFactor"),
    dispersion = read_named_values(file.path(root, "selected_dispersions.tsv"), "dispersion"),
    beta_no = read_matrix(file.path(root, "beta_no_optim.tsv")),
    beta_force = read_matrix(file.path(root, "beta_force_optim.tsv"))
  )
  assign(contrast, value, envir = loaded)
  value
}

for (seed_idx in seq_len(nrow(hardest))) {
  seed <- hardest[seed_idx, , drop = FALSE]
  contrast <- seed$contrast
  gene <- seed$gene
  data <- load_contrast(contrast)
  if (!gene %in% rownames(data$counts)) {
    next
  }
  counts <- as.integer(data$counts[gene, ])
  design <- data$design
  size_factors <- as.numeric(data$size_factors[colnames(data$counts)])
  if (any(!is.finite(size_factors))) {
    size_factors <- as.numeric(data$size_factors)
  }
  dispersion <- as.numeric(data$dispersion[[gene]])
  beta_no <- as.numeric(data$beta_no[gene, ])
  beta_force <- as.numeric(data$beta_force[gene, ])
  names(beta_no) <- colnames(data$beta_no)
  names(beta_force) <- colnames(data$beta_force)

  count_variants <- make_count_variants(counts, design)
  dispersion_variants <- c(real = 1, lower_dispersion = 0.25, higher_dispersion = 4, extreme_dispersion = 16)
  start_variants <- make_start_variants(beta_no, beta_force)
  ridge_log2 <- rep(1e-6, ncol(design))

  for (count_kind in names(count_variants)) {
    variant_counts <- pmin(pmax(count_variants[[count_kind]], 0L), .Machine$integer.max)
    for (dispersion_kind in names(dispersion_variants)) {
      variant_dispersion <- dispersion * dispersion_variants[[dispersion_kind]]
      if (!is.finite(variant_dispersion) || variant_dispersion <= 0) {
        next
      }
      for (start_kind in names(start_variants)) {
        start <- start_variants[[start_kind]]
        total <- sum(variant_counts)
        zero_fraction <- mean(variant_counts == 0)
        case_idx <- case_idx + 1L
        candidate <- data.frame(
          pool_id = sprintf("pool_%06d", case_idx),
          contrast = contrast,
          source_gene = gene,
          seed_rank = seed_idx,
          source_hard_score = seed$hardScore,
          source_force_beta_delta = seed$maxAbsBetaNoOptimVsForceOptim,
          count_kind = count_kind,
          dispersion_kind = dispersion_kind,
          start_kind = start_kind,
          n_samples = length(variant_counts),
          n_coefficients = ncol(design),
          total_count = total,
          zero_fraction = zero_fraction,
          dispersion = variant_dispersion,
          total_bin = floor(log10(total + 1)),
          zero_bin = floor(zero_fraction * 10),
          dispersion_bin = floor(log10(variant_dispersion) * 2),
          start_bound_count = sum(abs(start) >= upper - 1e-8),
          start_objective = objective_value(start, variant_counts, design, size_factors, variant_dispersion, ridge_log2),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        candidate$variant_risk <- 0
        candidate$variant_risk <- candidate$variant_risk + ifelse(grepl("edge|near_", start_kind), 4, 0)
        candidate$variant_risk <- candidate$variant_risk + ifelse(count_kind %in% c("condition_low", "condition_high"), 3, 0)
        candidate$variant_risk <- candidate$variant_risk + ifelse(dispersion_kind %in% c("higher_dispersion", "extreme_dispersion"), 3, 0)
        candidate$selection_score <- as.numeric(candidate$source_hard_score) +
          log10(1 + as.numeric(candidate$source_force_beta_delta)) +
          candidate$variant_risk +
          0.01 * candidate$zero_bin
        candidate$signature <- stable_case_signature(candidate)
        case_pool[[length(case_pool) + 1L]] <- candidate
        case_inputs[[candidate$pool_id]] <- list(
          counts = variant_counts,
          design = design,
          size_factors = size_factors,
          dispersion = variant_dispersion,
          ridge_log2 = ridge_log2
        )
        case_coefficients[[candidate$pool_id]] <- list(
          coefficient = colnames(design),
          beta_start = start,
          beta_no_optim = beta_no,
          beta_force_optim = beta_force
        )
      }
    }
  }
}

pool <- do.call(rbind, case_pool)
pool <- pool[is.finite(pool$start_objective), , drop = FALSE]
pool <- pool[order(-pool$selection_score, pool$signature), , drop = FALSE]
pool <- pool[!duplicated(pool$signature), , drop = FALSE]
pool$stratum <- paste(pool$count_kind, pool$dispersion_kind, pool$start_kind, sep = "|")
strata <- split(pool, pool$stratum)
strata <- lapply(strata, function(frame) frame[order(-frame$selection_score), , drop = FALSE])
strata <- strata[order(vapply(strata, function(frame) max(frame$selection_score), numeric(1)), decreasing = TRUE)]
selected_rows <- list()
seen_pool_id <- character()
source_counts <- new.env(parent = emptyenv())
while (length(selected_rows) < n_cases && length(strata) > 0) {
  progressed <- FALSE
  for (name in names(strata)) {
    frame <- strata[[name]]
    while (nrow(frame) > 0) {
      source_key <- paste(frame$contrast[[1]], frame$source_gene[[1]], sep = "|")
      used <- if (exists(source_key, envir = source_counts, inherits = FALSE)) {
        get(source_key, envir = source_counts, inherits = FALSE)
      } else {
        0L
      }
      if (!(frame$pool_id[[1]] %in% seen_pool_id) && used < max_per_source) {
        break
      }
      frame <- frame[-1, , drop = FALSE]
    }
    if (nrow(frame) == 0) {
      strata[[name]] <- frame
      next
    }
    selected_rows[[length(selected_rows) + 1L]] <- frame[1, , drop = FALSE]
    seen_pool_id <- c(seen_pool_id, frame$pool_id[[1]])
    source_key <- paste(frame$contrast[[1]], frame$source_gene[[1]], sep = "|")
    used <- if (exists(source_key, envir = source_counts, inherits = FALSE)) {
      get(source_key, envir = source_counts, inherits = FALSE)
    } else {
      0L
    }
    assign(source_key, used + 1L, envir = source_counts)
    strata[[name]] <- frame[-1, , drop = FALSE]
    progressed <- TRUE
    if (length(selected_rows) >= n_cases) {
      break
    }
  }
  strata <- strata[vapply(strata, nrow, integer(1)) > 0]
  if (!progressed) {
    break
  }
}
selected <- bind_rows_fill(selected_rows)
selected$case_id <- sprintf("lbfgsb_stress_%04d", seq_len(nrow(selected)))

case_rows <- list()
sample_rows <- list()
coefficient_rows <- list()

for (row_idx in seq_len(nrow(selected))) {
  row <- selected[row_idx, , drop = FALSE]
  pool_id <- row$pool_id
  input <- case_inputs[[pool_id]]
  coefs <- case_coefficients[[pool_id]]
  fn <- function(beta) {
    objective_value(
      beta,
      input$counts,
      input$design,
      input$size_factors,
      input$dispersion,
      input$ridge_log2
    )
  }
  opt <- tryCatch(
    stats::optim(
      par = coefs$beta_start,
      fn = fn,
      method = "L-BFGS-B",
      lower = rep(lower, length(coefs$beta_start)),
      upper = rep(upper, length(coefs$beta_start)),
      control = list(maxit = maxit, factr = 1e7, pgtol = 0, lmm = 5)
    ),
    error = function(error) {
      list(
        par = rep(NA_real_, length(coefs$beta_start)),
        value = NA_real_,
        counts = c("function" = NA_integer_, "gradient" = NA_integer_),
        convergence = 999L,
        message = conditionMessage(error)
      )
    }
  )
  grad <- if (all(is.finite(opt$par))) numeric_gradient(opt$par, fn, lower, upper) else rep(NA_real_, length(coefs$beta_start))
  row$case_id <- selected$case_id[[row_idx]]
  row$optim_value <- opt$value
  row$optim_convergence <- opt$convergence
  row$optim_message <- if (is.null(opt$message)) "" else opt$message
  row$optim_fn_count <- unname(opt$counts[["function"]])
  row$optim_gr_count <- unname(opt$counts[["gradient"]])
  row$projected_gradient_norm <- projected_grad_norm(opt$par, grad, lower, upper)
  row$max_abs_start_to_optim <- max(abs(coefs$beta_start - opt$par), na.rm = TRUE)
  row$max_abs_force_to_optim <- max(abs(coefs$beta_force_optim - opt$par), na.rm = TRUE)
  case_rows[[length(case_rows) + 1L]] <- row

  sample_frame <- data.frame(
    case_id = row$case_id,
    sample_index_1based = seq_along(input$counts),
    count = input$counts,
    size_factor = input$size_factors,
    dispersion = input$dispersion,
    weight = 1,
    input$design,
    check.names = FALSE
  )
  sample_rows[[length(sample_rows) + 1L]] <- sample_frame

  coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
    case_id = row$case_id,
    coefficient_index_1based = seq_along(coefs$coefficient),
    coefficient = coefs$coefficient,
    lower = lower,
    upper = upper,
    ridge_log2 = input$ridge_log2,
    beta_start = coefs$beta_start,
    beta_no_optim = coefs$beta_no_optim,
    beta_force_optim = coefs$beta_force_optim,
    optim_par = opt$par,
    optim_numeric_gradient = grad,
    check.names = FALSE
  )
}

cases <- bind_rows_fill(case_rows)
samples <- bind_rows_fill(sample_rows)
coefficients <- bind_rows_fill(coefficient_rows)

write_tsv(cases, file.path(out_dir, "cases.tsv"))
write_tsv(coefficients, file.path(out_dir, "coefficients.tsv"))
write_tsv(samples, gzfile(file.path(out_dir, "samples.tsv.gz")))
write_tsv(pool, file.path(out_dir, "candidate_pool.tsv"))

manifest <- data.frame(
  key = c(
    "source_root", "n_requested_cases", "n_written_cases", "n_seed_rows",
    "n_candidate_pool", "max_per_source", "lower", "upper", "maxit", "factr",
    "pgtol", "lmm", "r_version"
  ),
  value = c(
    source_root, n_cases, nrow(cases), nrow(hardest),
    nrow(pool), max_per_source, lower, upper, maxit, 1e7, 0, 5,
    R.version.string
  ),
  check.names = FALSE
)
write_tsv(manifest, file.path(out_dir, "manifest.tsv"))

readme <- c(
  "# Synthetic L-BFGS-B Stress Fixtures",
  "",
  "This ignored bundle expands the worst real DESeq2 GLM beta optimizer rows into a compact, nonredundant set of bounded L-BFGS-B objective cases.",
  "",
  "The generator starts from `global_hardest_512.tsv` in the real hard-row bundle, takes the highest-scoring source rows, builds nearby count/dispersion/start variants, removes duplicate optimizer-shape signatures, and runs base R `optim(..., method = \"L-BFGS-B\")` on the selected cases.",
  "",
  "Files:",
  "",
  "- `cases.tsv`: one row per fixture with source row, synthetic variant labels, start objective, R optim value, convergence code, function counts, projected gradient norm, and distance to source targets.",
  "- `coefficients.tsv`: coefficient-level lower/upper bounds, ridge, start vector, source no-optim/force-optim beta, R optim par, and numeric gradient at the final par.",
  "- `samples.tsv.gz`: sample-level count, size factor, dispersion, weight, and design row for each case.",
  "- `candidate_pool.tsv`: all nonduplicate candidates considered before selecting the final cases.",
  "- `manifest.tsv`: generation settings.",
  "",
  "Good first targets for an optimizer port are cases with `optim_convergence == 0`, high `max_abs_start_to_optim`, high `optim_fn_count`, and small `projected_gradient_norm`."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote ", nrow(cases), " stress fixtures to ", normalizePath(out_dir))
