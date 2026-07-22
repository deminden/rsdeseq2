#!/usr/bin/env Rscript

# Regenerate DESeq2 real-data references without modifying the source study.
#
# The output study exposes the immutable source inputs and group definitions via
# directory symlinks, while all newly computed DESeq2 outputs and provenance are
# written below --output-study-root.

`%||%` <- function(left, right) if (is.null(left)) right else left

usage <- function(status = 0L) {
  stream <- if (status == 0L) stdout() else stderr()
  cat(
    paste0(
      "Usage:\n",
      "  Rscript scripts/generate_real_data_references.R \\\n",
      "    --source-study-root PATH --output-study-root PATH \\\n",
      "    [--contrast tissue:null_type:rep ...] [--tissue tissue ...]\n\n",
      "Selection:\n",
      "  --contrast SPEC          Regenerate one contrast; may be repeated.\n",
      "  --tissue NAME            Regenerate normalized counts for one tissue;\n",
      "                           may be repeated. With --all-existing-groups,\n",
      "                           restrict discovered contrasts to these tissues.\n",
      "  --all-existing-groups    Regenerate every source *_groups.tsv contrast\n",
      "                           (or every one for the selected --tissue values).\n\n",
      "Safety and execution:\n",
      "  --require-r-version VER  Refuse to run under a different R version.\n",
      "  --contrast-size-factors MODE\n",
      "                           Use 'estimate' (the source study's cached path;\n",
      "                           default) or 'full' (reuse full-tissue factors).\n",
      "  --workers N              Forked contrast workers per tissue (default: 1).\n",
      "  --resume                 Reuse complete files already in the output root.\n",
      "  --force                  Replace existing files only in the output root.\n",
      "  --help                   Show this help.\n\n",
      "At least one --contrast, --tissue, or --all-existing-groups option is required.\n"
    ),
    file = stream
  )
  quit(save = "no", status = status, runLast = FALSE)
}

parse_cli <- function(args) {
  cfg <- list(
    source_study_root = NULL,
    output_study_root = NULL,
    contrasts = character(),
    tissues = character(),
    all_existing_groups = FALSE,
    required_r_version = NULL,
    contrast_size_factors = "estimate",
    workers = 1L,
    resume = FALSE,
    force = FALSE
  )

  take_value <- function(i, option) {
    if (i >= length(args) || startsWith(args[[i + 1L]], "--")) {
      stop(option, " requires a value", call. = FALSE)
    }
    args[[i + 1L]]
  }

  i <- 1L
  while (i <= length(args)) {
    option <- args[[i]]
    if (option == "--help") {
      usage(0L)
    } else if (option == "--source-study-root") {
      cfg$source_study_root <- take_value(i, option)
      i <- i + 1L
    } else if (option == "--output-study-root") {
      cfg$output_study_root <- take_value(i, option)
      i <- i + 1L
    } else if (option == "--contrast") {
      cfg$contrasts <- c(cfg$contrasts, take_value(i, option))
      i <- i + 1L
    } else if (option == "--tissue") {
      cfg$tissues <- c(cfg$tissues, take_value(i, option))
      i <- i + 1L
    } else if (option == "--all-existing-groups") {
      cfg$all_existing_groups <- TRUE
    } else if (option == "--require-r-version") {
      cfg$required_r_version <- take_value(i, option)
      i <- i + 1L
    } else if (option == "--contrast-size-factors") {
      cfg$contrast_size_factors <- tolower(take_value(i, option))
      i <- i + 1L
    } else if (option == "--workers") {
      cfg$workers <- suppressWarnings(as.integer(take_value(i, option)))
      i <- i + 1L
    } else if (option == "--resume") {
      cfg$resume <- TRUE
    } else if (option == "--force") {
      cfg$force <- TRUE
    } else {
      stop("unknown option: ", option, call. = FALSE)
    }
    i <- i + 1L
  }

  if (is.null(cfg$source_study_root) || !nzchar(cfg$source_study_root)) {
    stop("--source-study-root is required", call. = FALSE)
  }
  if (is.null(cfg$output_study_root) || !nzchar(cfg$output_study_root)) {
    stop("--output-study-root is required", call. = FALSE)
  }
  if (length(cfg$contrasts) == 0L && length(cfg$tissues) == 0L &&
      !cfg$all_existing_groups) {
    stop(
      "select work with --contrast, --tissue, or --all-existing-groups",
      call. = FALSE
    )
  }
  if (is.na(cfg$workers) || cfg$workers < 1L) {
    stop("--workers must be a positive integer", call. = FALSE)
  }
  if (!cfg$contrast_size_factors %in% c("estimate", "full")) {
    stop("--contrast-size-factors must be 'estimate' or 'full'", call. = FALSE)
  }
  if (cfg$resume && cfg$force) {
    stop("--resume and --force are mutually exclusive", call. = FALSE)
  }
  cfg
}

canonical_path <- function(path, must_work) {
  normalizePath(path.expand(path), winslash = "/", mustWork = must_work)
}

path_is_within <- function(path, parent) {
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

validate_component <- function(value, label) {
  value <- tolower(trimws(value))
  if (!nzchar(value) || !grepl("^[a-z0-9_]+$", value)) {
    stop(label, " must contain only letters, digits, and underscores: ", value,
         call. = FALSE)
  }
  value
}

empty_contrasts <- function() {
  data.frame(
    tissue = character(),
    null_type = character(),
    rep = integer(),
    group_path = character(),
    stringsAsFactors = FALSE
  )
}

contrast_stem <- function(tissue, null_type, rep) {
  sprintf("%s_%s_rep%02d", tissue, null_type, rep)
}

parse_contrast <- function(spec, group_dir) {
  parts <- strsplit(spec, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 3L) {
    stop("contrast must be tissue:null_type:rep, got: ", spec, call. = FALSE)
  }
  tissue <- validate_component(parts[[1L]], "contrast tissue")
  null_type <- validate_component(parts[[2L]], "contrast null type")
  if (!grepl("^[0-9]+$", parts[[3L]])) {
    stop("contrast rep must be a positive integer, got: ", spec, call. = FALSE)
  }
  rep <- as.integer(parts[[3L]])
  if (is.na(rep) || rep < 1L) {
    stop("contrast rep must be a positive integer, got: ", spec, call. = FALSE)
  }
  stem <- contrast_stem(tissue, null_type, rep)
  data.frame(
    tissue = tissue,
    null_type = null_type,
    rep = rep,
    group_path = file.path(group_dir, paste0(stem, "_groups.tsv")),
    stringsAsFactors = FALSE
  )
}

available_tissues <- function(input_dir) {
  files <- list.files(
    input_dir,
    pattern = "_raw_counts[.]tsv[.]gz$",
    full.names = FALSE
  )
  sort(sub("_raw_counts[.]tsv[.]gz$", "", files))
}

discover_group_contrasts <- function(group_dir, tissues, tissue_filter = character()) {
  files <- sort(list.files(
    group_dir,
    pattern = "_groups[.]tsv$",
    full.names = FALSE
  ))
  if (length(files) == 0L) {
    return(empty_contrasts())
  }

  # Tissue and null-design names both contain underscores. Resolve the tissue
  # from the known raw-count filenames, preferring the longest matching prefix.
  tissues <- tissues[order(nchar(tissues), decreasing = TRUE)]
  rows <- list()
  for (file in files) {
    candidates <- tissues[vapply(
      tissues,
      function(tissue) startsWith(file, paste0(tissue, "_")),
      logical(1)
    )]
    if (length(candidates) == 0L) {
      warning("ignoring group file without a matching raw-count tissue: ", file)
      next
    }
    tissue <- candidates[[1L]]
    if (length(tissue_filter) > 0L && !tissue %in% tissue_filter) {
      next
    }
    remainder <- substring(file, nchar(tissue) + 2L)
    match <- regexec("^(.+)_rep([0-9]+)_groups[.]tsv$", remainder)
    fields <- regmatches(remainder, match)[[1L]]
    if (length(fields) != 3L) {
      warning("ignoring unrecognized group filename: ", file)
      next
    }
    null_type <- validate_component(fields[[2L]], "discovered null type")
    rep <- suppressWarnings(as.integer(fields[[3L]]))
    if (is.na(rep) || rep < 1L) {
      warning("ignoring group filename with invalid rep: ", file)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      tissue = tissue,
      null_type = null_type,
      rep = rep,
      group_path = file.path(group_dir, file),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) empty_contrasts() else do.call(rbind, rows)
}

assert_not_symlink <- function(path, label) {
  target <- Sys.readlink(path)
  if (length(target) > 0L && !is.na(target) && nzchar(target)) {
    stop(label, " must not be a symlink: ", path, " -> ", target, call. = FALSE)
  }
}

expose_source_directory <- function(source_dir, output_dir) {
  link_target <- Sys.readlink(output_dir)
  destination_exists <- file.exists(output_dir) ||
    (length(link_target) > 0L && !is.na(link_target) && nzchar(link_target))
  if (destination_exists) {
    if (!dir.exists(output_dir) ||
        !identical(canonical_path(output_dir, TRUE), canonical_path(source_dir, TRUE))) {
      stop(
        "output exposure already exists but does not resolve to the source: ",
        output_dir,
        call. = FALSE
      )
    }
    return(invisible(output_dir))
  }
  if (!isTRUE(file.symlink(canonical_path(source_dir, TRUE), output_dir))) {
    stop("could not expose source directory with a symlink: ", output_dir,
         call. = FALSE)
  }
  invisible(output_dir)
}

promote_file <- function(temporary, destination, force) {
  assert_not_symlink(destination, "output file")
  if (dir.exists(destination)) {
    stop("output file path is a directory: ", destination, call. = FALSE)
  }
  exists <- file.exists(destination)
  if (exists && !force) {
    stop("output already exists; pass --force to replace it: ", destination,
         call. = FALSE)
  }

  backup <- NULL
  if (exists) {
    backup <- tempfile(
      pattern = paste0(".", basename(destination), ".backup-"),
      tmpdir = dirname(destination)
    )
    if (!file.rename(destination, backup)) {
      stop("could not move existing output aside: ", destination, call. = FALSE)
    }
  }

  if (!file.rename(temporary, destination)) {
    if (!is.null(backup)) {
      file.rename(backup, destination)
    }
    stop("could not install generated output: ", destination, call. = FALSE)
  }
  if (!is.null(backup)) {
    unlink(backup)
  }
  invisible(destination)
}

atomic_write_table <- function(frame, path, force, gzip = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), ".partial-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  connection <- if (gzip) gzfile(temporary, "wt") else file(temporary, "wt")
  tryCatch(
    utils::write.table(
      frame,
      connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    ),
    finally = close(connection)
  )
  promote_file(temporary, path, force)
}

validate_gzip_tsv <- function(path, expected_columns, expected_rows) {
  if (!file.exists(path)) {
    return(list(complete = FALSE, reason = "file does not exist"))
  }
  assert_not_symlink(path, "resumable output file")
  info <- file.info(path)
  if (is.na(info$size[[1L]]) || info$size[[1L]] <= 0) {
    return(list(complete = FALSE, reason = "file is empty"))
  }

  connection <- NULL
  warning_messages <- character()
  error_message <- NULL
  row_count <- 0L
  tryCatch(
    withCallingHandlers({
      connection <- gzfile(path, "rt")
      header <- readLines(connection, n = 1L, warn = TRUE)
      if (length(header) != 1L) {
        stop("missing header", call. = FALSE)
      }
      observed_columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
      if (!identical(observed_columns, expected_columns)) {
        stop("header does not match the expected columns", call. = FALSE)
      }
      repeat {
        lines <- readLines(connection, n = 50000L, warn = TRUE)
        row_count <- row_count + length(lines)
        if (row_count > expected_rows) {
          stop("table has more rows than expected", call. = FALSE)
        }
        if (length(lines) == 0L) {
          break
        }
      }
      close(connection)
      connection <- NULL
      if (row_count != expected_rows) {
        stop(
          sprintf("table has %d rows; expected %d", row_count, expected_rows),
          call. = FALSE
        )
      }
    }, warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }),
    error = function(error) {
      error_message <<- conditionMessage(error)
    }
  )
  if (!is.null(connection)) {
    try(close(connection), silent = TRUE)
  }
  if (length(warning_messages) > 0L) {
    return(list(
      complete = FALSE,
      reason = paste(unique(warning_messages), collapse = "; ")
    ))
  }
  if (!is.null(error_message)) {
    return(list(complete = FALSE, reason = error_message))
  }
  list(complete = TRUE, reason = "complete")
}

read_counts <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection), add = TRUE)
  frame <- utils::read.delim(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(frame) < 2L || nrow(frame) < 1L) {
    stop("raw-count file is empty or has no sample columns: ", path, call. = FALSE)
  }
  genes <- as.character(frame[[1L]])
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("raw-count file has missing, empty, or duplicate gene identifiers: ", path,
         call. = FALSE)
  }
  frame[[1L]] <- NULL
  if (any(!nzchar(names(frame))) || anyDuplicated(names(frame))) {
    stop("raw-count file has empty or duplicate sample identifiers: ", path,
         call. = FALSE)
  }
  matrix <- as.matrix(frame)
  suppressWarnings(storage.mode(matrix) <- "numeric")
  rownames(matrix) <- genes
  matrix
}

prepare_deseq2_counts <- function(count_matrix, label, quiet = FALSE) {
  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "numeric"
  rounded <- round(count_matrix)

  bad_value <- !is.finite(rounded) | rounded < 0 |
    rounded > .Machine$integer.max
  keep_valid <- rowSums(bad_value) == 0L
  if (any(!keep_valid) && !quiet) {
    examples <- paste(utils::head(rownames(rounded)[!keep_valid], 5L), collapse = ", ")
    warning(sprintf(
      "%s: dropped %d genes with non-finite, negative, or > integer-range counts%s",
      label,
      sum(!keep_valid),
      if (nzchar(examples)) paste0(" (e.g. ", examples, ")") else ""
    ))
  }
  rounded <- rounded[keep_valid, , drop = FALSE]

  keep_count <- rowSums(rounded) >= 1
  rounded <- rounded[keep_count, , drop = FALSE]
  if (nrow(rounded) == 0L) {
    stop(label, ": count validation removed every gene", call. = FALSE)
  }
  if (!quiet) {
    message(sprintf("  %s: %d validated genes x %d samples", label,
                    nrow(rounded), ncol(rounded)))
  }
  rounded
}

read_groups <- function(path, sample_names) {
  groups <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    na.strings = "NA"
  )
  required <- c("sample_id", "condition", "retained")
  missing <- setdiff(required, names(groups))
  if (length(missing) > 0L) {
    stop("group file is missing columns ", paste(missing, collapse = ", "),
         ": ", path, call. = FALSE)
  }
  retained <- toupper(as.character(groups$retained)) == "TRUE" &
    as.character(groups$condition) %in% c("A", "B")
  retained[is.na(retained)] <- FALSE
  groups <- groups[retained, , drop = FALSE]
  groups$sample_id <- as.character(groups$sample_id)
  groups$condition <- as.character(groups$condition)
  if (nrow(groups) == 0L) {
    stop("group file has no retained A/B samples: ", path, call. = FALSE)
  }
  if (anyNA(groups$sample_id) || any(!nzchar(groups$sample_id)) ||
      anyDuplicated(groups$sample_id)) {
    stop("retained group rows have missing, empty, or duplicate sample IDs: ", path,
         call. = FALSE)
  }
  absent <- setdiff(groups$sample_id, sample_names)
  if (length(absent) > 0L) {
    stop(
      "group file contains samples absent from raw counts: ",
      paste(utils::head(absent, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  condition_counts <- table(groups$condition)
  if (!all(c("A", "B") %in% names(condition_counts)) ||
      any(condition_counts[c("A", "B")] < 2L)) {
    stop("contrast requires at least two retained samples in each condition: ", path,
         call. = FALSE)
  }
  groups
}

infer_design <- function(groups, null_type) {
  col_data <- data.frame(
    row.names = groups$sample_id,
    condition = factor(groups$condition, levels = c("A", "B"))
  )
  design <- ~ condition
  design_name <- "condition"

  blocked <- grepl("blocked_permutation$", null_type)
  if (blocked && "perm_block" %in% names(groups)) {
    perm_block <- factor(groups$perm_block)
    can_try_blocked <- !anyNA(perm_block) && nlevels(perm_block) > 0L
    if (can_try_blocked) {
      candidate_data <- col_data
      candidate_data$perm_block <- droplevels(perm_block)
      candidate_matrix <- tryCatch(
        stats::model.matrix(~ perm_block + condition, data = candidate_data),
        error = function(error) NULL
      )
      block_counts <- table(candidate_data$perm_block, candidate_data$condition)
      block_valid <- nrow(block_counts) > 0L &&
        all(rowSums(block_counts > 0L) == 2L)
      design_valid <- !is.null(candidate_matrix) &&
        qr(candidate_matrix)$rank == ncol(candidate_matrix)
      if (block_valid && design_valid) {
        col_data <- candidate_data
        design <- ~ perm_block + condition
        design_name <- "perm_block + condition"
      }
    }
  }
  list(col_data = col_data, formula = design, name = design_name)
}

write_normalized_counts <- function(matrix, path, force) {
  frame <- data.frame(
    gene = rownames(matrix),
    as.data.frame(matrix, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  atomic_write_table(frame, path, force = force, gzip = TRUE)
}

write_results <- function(results, path, force) {
  frame <- as.data.frame(results, optional = TRUE)
  frame$gene <- rownames(frame)
  atomic_write_table(frame, path, force = force, gzip = TRUE)
}

process_contrast <- function(spec, full_counts, full_size_factors,
                             contrast_size_factors, output_dir, force, resume,
                             reuse_outputs) {
  stem <- contrast_stem(spec$tissue[[1L]], spec$null_type[[1L]], spec$rep[[1L]])
  message("    ", stem)
  groups <- read_groups(spec$group_path[[1L]], colnames(full_counts))
  counts <- full_counts[, groups$sample_id, drop = FALSE]

  design <- infer_design(groups, spec$null_type[[1L]])
  path <- file.path(output_dir, paste0(stem, "_deseq2_results.tsv.gz"))
  if (reuse_outputs && file.exists(path)) {
    validation <- validate_gzip_tsv(
      path,
      c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene"),
      nrow(counts)
    )
    if (validation$complete) {
      message("      reusing validated output")
      return(list(
        path = canonical_path(path, TRUE),
        tissue = spec$tissue[[1L]],
        null_type = spec$null_type[[1L]],
        rep = spec$rep[[1L]],
        design = design$name,
        n_genes = nrow(counts),
        n_samples = ncol(counts)
      ))
    }
    message("      replacing incomplete output: ", validation$reason)
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = design$col_data,
    design = design$formula
  )

  if (contrast_size_factors == "estimate") {
    # The source study loads its saved normalization cache before fitting these
    # contrasts, so no in-memory full-tissue DESeqDataSet is available. Its
    # recorded execution path therefore calls DESeq() on each retained split,
    # including split-specific size-factor estimation and the same conditional
    # Cook's-distance outlier refitting.
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
  } else {
    size_factors <- full_size_factors[colnames(counts)]
    if (length(size_factors) != ncol(counts) || anyNA(size_factors) ||
        any(!is.finite(size_factors)) || any(size_factors <= 0)) {
      stop(stem, ": could not select valid full-tissue size factors", call. = FALSE)
    }
    DESeq2::sizeFactors(dds) <- unname(size_factors)
    dds <- DESeq2::estimateDispersions(dds, fitType = "parametric", quiet = TRUE)
    dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
  }
  results <- DESeq2::results(dds, contrast = c("condition", "B", "A"))

  write_results(results, path, force || resume)
  list(
    path = canonical_path(path, TRUE),
    tissue = spec$tissue[[1L]],
    null_type = spec$null_type[[1L]],
    rep = spec$rep[[1L]],
    design = design$name,
    n_genes = nrow(results),
    n_samples = ncol(counts)
  )
}

file_record <- function(record_type, name, value = "", tissue = "",
                        null_type = "", rep = NA_integer_, path = "",
                        include_md5 = FALSE) {
  size_bytes <- NA_real_
  mtime_utc <- ""
  md5 <- ""
  if (nzchar(path) && file.exists(path)) {
    info <- file.info(path)
    size_bytes <- unname(info$size[[1L]])
    mtime_utc <- format(info$mtime[[1L]], "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    if (include_md5) {
      md5 <- unname(as.character(tools::md5sum(path)[[1L]]))
    }
  }
  data.frame(
    record_type = record_type,
    name = name,
    value = as.character(value),
    tissue = as.character(tissue),
    null_type = as.character(null_type),
    rep = as.integer(rep),
    path = as.character(path),
    size_bytes = size_bytes,
    mtime_utc = mtime_utc,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

metadata_record <- function(name, value) {
  file_record("metadata", name, value = value)
}

read_reference_manifest <- function(path) {
  manifest <- tryCatch(
    utils::read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      na.strings = "NA"
    ),
    error = function(error) {
      stop("cannot read resume manifest: ", conditionMessage(error), call. = FALSE)
    }
  )
  required <- c(
    "record_type", "name", "value", "tissue", "null_type", "rep",
    "path", "size_bytes", "mtime_utc", "md5"
  )
  if (!identical(names(manifest), required)) {
    stop("resume manifest has an unexpected header", call. = FALSE)
  }
  manifest
}

manifest_metadata_value <- function(manifest, name) {
  hit <- manifest$record_type == "metadata" & manifest$name == name
  hit[is.na(hit)] <- FALSE
  if (sum(hit) != 1L) {
    stop("resume manifest must contain exactly one metadata row for ", name,
         call. = FALSE)
  }
  as.character(manifest$value[hit][[1L]])
}

manifest_input_key <- function(frame) {
  rep <- ifelse(is.na(frame$rep), "", as.character(as.integer(frame$rep)))
  paste(
    frame$name,
    ifelse(is.na(frame$tissue), "", frame$tissue),
    ifelse(is.na(frame$null_type), "", frame$null_type),
    rep,
    frame$path,
    sep = "\034"
  )
}

validate_resume_manifest <- function(path, expected_metadata, expected_inputs) {
  assert_not_symlink(path, "resume manifest")
  manifest <- read_reference_manifest(path)
  for (name in names(expected_metadata)) {
    observed <- manifest_metadata_value(manifest, name)
    expected <- as.character(expected_metadata[[name]])
    if (!identical(observed, expected)) {
      stop(
        "resume manifest mismatch for ", name, ": expected '", expected,
        "', found '", observed, "'",
        call. = FALSE
      )
    }
  }
  status_hit <- manifest$record_type == "metadata" & manifest$name == "run_status"
  status_hit[is.na(status_hit)] <- FALSE
  if (any(status_hit)) {
    status <- manifest_metadata_value(manifest, "run_status")
    if (!status %in% c("in_progress", "complete")) {
      stop("resume manifest has invalid run_status: ", status, call. = FALSE)
    }
  }

  observed_inputs <- manifest[manifest$record_type == "input", , drop = FALSE]
  expected_inputs <- expected_inputs[expected_inputs$record_type == "input", , drop = FALSE]
  observed_keys <- manifest_input_key(observed_inputs)
  expected_keys <- manifest_input_key(expected_inputs)
  if (anyDuplicated(observed_keys) || anyDuplicated(expected_keys) ||
      !setequal(observed_keys, expected_keys)) {
    stop("resume manifest selected-input set does not match this run", call. = FALSE)
  }
  observed_order <- match(expected_keys, observed_keys)
  observed_md5 <- as.character(observed_inputs$md5[observed_order])
  expected_md5 <- as.character(expected_inputs$md5)
  if (anyNA(observed_md5) || anyNA(expected_md5) ||
      any(!nzchar(observed_md5)) || any(!nzchar(expected_md5)) ||
      !identical(observed_md5, expected_md5)) {
    stop("resume manifest input hashes do not match this run", call. = FALSE)
  }
  invisible(manifest)
}

cfg <- tryCatch(
  parse_cli(commandArgs(trailingOnly = TRUE)),
  error = function(error) {
    message("error: ", conditionMessage(error))
    usage(2L)
  }
)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required but is not installed", call. = FALSE)
}

r_version <- paste(R.version$major, R.version$minor, sep = ".")
if (!is.null(cfg$required_r_version) &&
    !identical(r_version, cfg$required_r_version)) {
  stop(
    "R version mismatch: required ", cfg$required_r_version,
    ", running ", r_version,
    call. = FALSE
  )
}

source_root <- canonical_path(cfg$source_study_root, TRUE)
if (!dir.exists(source_root)) {
  stop("source study root is not a directory: ", source_root, call. = FALSE)
}
output_root_candidate <- canonical_path(cfg$output_study_root, FALSE)
if (path_is_within(output_root_candidate, source_root) ||
    path_is_within(source_root, output_root_candidate)) {
  stop(
    "source and output study roots must be separate, non-nested directories",
    call. = FALSE
  )
}
if (file.exists(output_root_candidate) && !dir.exists(output_root_candidate)) {
  stop("output study root is not a directory: ", output_root_candidate,
       call. = FALSE)
}
dir.create(output_root_candidate, recursive = TRUE, showWarnings = FALSE)
output_root <- canonical_path(output_root_candidate, TRUE)
if (path_is_within(output_root, source_root) || path_is_within(source_root, output_root)) {
  stop(
    "source and output study roots resolve to nested directories",
    call. = FALSE
  )
}

source_input_dir <- file.path(source_root, "01_inputs")
source_group_dir <- file.path(source_root, "02_null_splits")
if (!dir.exists(source_input_dir) || !dir.exists(source_group_dir)) {
  stop("source study must contain 01_inputs and 02_null_splits directories",
       call. = FALSE)
}

known_tissues <- available_tissues(source_input_dir)
if (length(known_tissues) == 0L) {
  stop("source 01_inputs has no *_raw_counts.tsv.gz files", call. = FALSE)
}
requested_tissues <- unique(vapply(
  cfg$tissues,
  validate_component,
  character(1),
  label = "tissue"
))
unknown_tissues <- setdiff(requested_tissues, known_tissues)
if (length(unknown_tissues) > 0L) {
  stop("no raw counts found for tissue(s): ",
       paste(unknown_tissues, collapse = ", "), call. = FALSE)
}

explicit_contrasts <- if (length(cfg$contrasts) == 0L) {
  empty_contrasts()
} else {
  do.call(rbind, lapply(cfg$contrasts, parse_contrast, group_dir = source_group_dir))
}
discovered_contrasts <- if (cfg$all_existing_groups) {
  discover_group_contrasts(
    source_group_dir,
    known_tissues,
    tissue_filter = requested_tissues
  )
} else {
  empty_contrasts()
}
contrasts <- rbind(explicit_contrasts, discovered_contrasts)
if (nrow(contrasts) > 0L) {
  keys <- contrast_stem(contrasts$tissue, contrasts$null_type, contrasts$rep)
  contrasts <- contrasts[!duplicated(keys), , drop = FALSE]
  contrasts <- contrasts[order(
    contrasts$tissue,
    contrasts$null_type,
    contrasts$rep
  ), , drop = FALSE]
  rownames(contrasts) <- NULL
}
if (cfg$all_existing_groups && nrow(discovered_contrasts) == 0L) {
  stop("--all-existing-groups found no matching group files", call. = FALSE)
}
if (nrow(contrasts) > 0L) {
  missing_groups <- contrasts$group_path[!file.exists(contrasts$group_path)]
  if (length(missing_groups) > 0L) {
    stop("selected group file does not exist: ", missing_groups[[1L]],
         call. = FALSE)
  }
  unknown_contrast_tissues <- setdiff(unique(contrasts$tissue), known_tissues)
  if (length(unknown_contrast_tissues) > 0L) {
    stop("no raw counts found for contrast tissue(s): ",
         paste(unknown_contrast_tissues, collapse = ", "), call. = FALSE)
  }
}
selected_tissues <- sort(unique(c(requested_tissues, contrasts$tissue)))

expose_source_directory(source_input_dir, file.path(output_root, "01_inputs"))
expose_source_directory(source_group_dir, file.path(output_root, "02_null_splits"))
output_deseq_dir <- file.path(output_root, "02_deseq2_outputs")
if (file.exists(output_deseq_dir)) {
  if (!dir.exists(output_deseq_dir)) {
    stop("02_deseq2_outputs exists but is not a directory", call. = FALSE)
  }
  assert_not_symlink(output_deseq_dir, "output DESeq2 directory")
} else {
  dir.create(output_deseq_dir, recursive = TRUE, showWarnings = FALSE)
}

session <- sessionInfo()
bioconductor_version <- if (requireNamespace("BiocVersion", quietly = TRUE)) {
  as.character(packageVersion("BiocVersion"))
} else if (requireNamespace("BiocManager", quietly = TRUE)) {
  suppressWarnings(as.character(BiocManager::version()))
} else {
  "unavailable"
}
metadata <- list(
  manifest_schema = "rsdeseq2-real-data-reference-v1",
  run_status = "in_progress",
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  command = paste(vapply(commandArgs(), shQuote, character(1)), collapse = " "),
  source_study_root = source_root,
  output_study_root = output_root,
  r_version = r_version,
  r_version_string = R.version.string,
  deseq2_version = as.character(packageVersion("DESeq2")),
  bioconductor_version = bioconductor_version,
  biocparallel_version = if (requireNamespace("BiocParallel", quietly = TRUE)) {
    as.character(packageVersion("BiocParallel"))
  } else {
    "unavailable"
  },
  s4vectors_version = as.character(packageVersion("S4Vectors")),
  summarizedexperiment_version = as.character(packageVersion("SummarizedExperiment")),
  r_platform = R.version$platform,
  os = paste(unname(Sys.info()[c("sysname", "release", "machine")]), collapse = " "),
  blas = session$BLAS %||% "unavailable",
  lapack = session$LAPACK %||% "unavailable",
  locale = paste(session$locale, collapse = ";"),
  requested_workers = cfg$workers,
  contrast_size_factors = cfg$contrast_size_factors,
  logical_cores = parallel::detectCores(logical = TRUE),
  physical_cores = parallel::detectCores(logical = FALSE),
  omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = "unset"),
  openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = "unset"),
  mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", unset = "unset"),
  blis_num_threads = Sys.getenv("BLIS_NUM_THREADS", unset = "unset"),
  selected_tissues = paste(selected_tissues, collapse = ";"),
  selected_contrasts = if (nrow(contrasts) == 0L) "" else paste(
    contrast_stem(contrasts$tissue, contrasts$null_type, contrasts$rep),
    collapse = ";"
  ),
  all_existing_groups = cfg$all_existing_groups,
  initial_force = cfg$force,
  resume = cfg$resume,
  force = cfg$force
)

input_manifest_rows <- list()
for (tissue in selected_tissues) {
  path <- canonical_path(
    file.path(source_input_dir, paste0(tissue, "_raw_counts.tsv.gz")),
    TRUE
  )
  input_manifest_rows[[length(input_manifest_rows) + 1L]] <- file_record(
    "input", "raw_counts", tissue = tissue, path = path, include_md5 = TRUE
  )
}
if (nrow(contrasts) > 0L) {
  for (index in seq_len(nrow(contrasts))) {
    spec <- contrasts[index, , drop = FALSE]
    input_manifest_rows[[length(input_manifest_rows) + 1L]] <- file_record(
      "input",
      "groups",
      tissue = spec$tissue,
      null_type = spec$null_type,
      rep = spec$rep,
      path = canonical_path(spec$group_path, TRUE),
      include_md5 = TRUE
    )
  }
}
input_manifest <- do.call(rbind, input_manifest_rows)

manifest_path <- file.path(output_root, "reference_manifest.tsv")
target_output_paths <- c(
  file.path(output_deseq_dir, paste0(selected_tissues, "_norm_counts.tsv.gz")),
  if (nrow(contrasts) == 0L) character() else file.path(
    output_deseq_dir,
    paste0(contrast_stem(contrasts$tissue, contrasts$null_type, contrasts$rep),
           "_deseq2_results.tsv.gz")
  )
)
manifest_exists <- file.exists(manifest_path)
reuse_outputs <- cfg$resume
critical_metadata_names <- c(
  "manifest_schema", "source_study_root", "output_study_root", "r_version",
  "r_version_string", "deseq2_version", "bioconductor_version",
  "biocparallel_version", "s4vectors_version", "summarizedexperiment_version",
  "r_platform", "os", "blas", "lapack", "locale", "requested_workers",
  "contrast_size_factors", "omp_num_threads", "openblas_num_threads",
  "mkl_num_threads", "blis_num_threads", "selected_tissues",
  "selected_contrasts", "all_existing_groups"
)
if (cfg$resume) {
  if (manifest_exists) {
    previous_manifest <- validate_resume_manifest(
      manifest_path,
      metadata[critical_metadata_names],
      input_manifest
    )
    status_hit <- previous_manifest$record_type == "metadata" &
      previous_manifest$name == "run_status"
    status_hit[is.na(status_hit)] <- FALSE
    previous_status <- if (any(status_hit)) {
      manifest_metadata_value(previous_manifest, "run_status")
    } else {
      "complete"
    }
    initial_force_hit <- previous_manifest$record_type == "metadata" &
      previous_manifest$name == "initial_force"
    initial_force_hit[is.na(initial_force_hit)] <- FALSE
    previous_initial_force <- if (any(initial_force_hit)) {
      identical(manifest_metadata_value(previous_manifest, "initial_force"), "TRUE")
    } else {
      force_hit <- previous_manifest$record_type == "metadata" &
        previous_manifest$name == "force"
      force_hit[is.na(force_hit)] <- FALSE
      any(force_hit) && identical(
        manifest_metadata_value(previous_manifest, "force"),
        "TRUE"
      )
    }
    metadata$initial_force <- previous_initial_force
    if (identical(previous_status, "in_progress") && previous_initial_force) {
      reuse_outputs <- FALSE
      message(
        "The interrupted run began with --force; regenerating every selected ",
        "output rather than reusing files with ambiguous provenance"
      )
    }
  } else if (any(file.exists(target_output_paths))) {
    stop(
      "cannot safely resume existing outputs without a provenance manifest; ",
      "use --force or a new output root",
      call. = FALSE
    )
  }
} else if (manifest_exists && !cfg$force) {
  stop("output manifest already exists; pass --resume or --force", call. = FALSE)
} else if (!cfg$force && any(file.exists(target_output_paths))) {
  stop("selected output already exists; pass --resume or --force", call. = FALSE)
}

preflight_manifest <- do.call(rbind, c(
  lapply(names(metadata), function(name) metadata_record(name, metadata[[name]])),
  list(input_manifest)
))
atomic_write_table(
  preflight_manifest,
  manifest_path,
  force = manifest_exists && (cfg$resume || cfg$force),
  gzip = FALSE
)

message("Source study: ", source_root)
message("Output study: ", output_root)
message("R: ", r_version, "; DESeq2: ", as.character(packageVersion("DESeq2")))
message("Selected tissues: ", paste(selected_tissues, collapse = ", "))
message("Selected contrasts: ", nrow(contrasts))

output_records <- list()
for (tissue in selected_tissues) {
  message("\n[", tissue, "]")
  count_path <- file.path(source_input_dir, paste0(tissue, "_raw_counts.tsv.gz"))
  raw_counts <- read_counts(count_path)
  full_counts <- prepare_deseq2_counts(raw_counts, paste0(tissue, " full tissue"))
  full_col_data <- data.frame(
    row.names = colnames(full_counts),
    condition = factor(rep("all", ncol(full_counts)))
  )
  full_dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = full_counts,
    colData = full_col_data,
    design = ~ 1
  )
  full_dds <- DESeq2::estimateSizeFactors(full_dds)
  full_size_factors <- DESeq2::sizeFactors(full_dds)
  if (is.null(names(full_size_factors))) {
    names(full_size_factors) <- colnames(full_counts)
  }
  normalized <- DESeq2::counts(full_dds, normalized = TRUE)
  normalized_path <- file.path(output_deseq_dir, paste0(tissue, "_norm_counts.tsv.gz"))
  if (reuse_outputs && file.exists(normalized_path)) {
    validation <- validate_gzip_tsv(
      normalized_path,
      c("gene", colnames(normalized)),
      nrow(normalized)
    )
    if (validation$complete) {
      message("  reusing validated normalized-count output")
    } else {
      message("  replacing incomplete normalized-count output: ", validation$reason)
      write_normalized_counts(normalized, normalized_path, force = TRUE)
    }
  } else {
    write_normalized_counts(normalized, normalized_path, cfg$force || cfg$resume)
  }
  output_records[[length(output_records) + 1L]] <- list(
    path = canonical_path(normalized_path, TRUE),
    tissue = tissue,
    null_type = "",
    rep = NA_integer_,
    design = "~ 1 (size factors only)",
    n_genes = nrow(normalized),
    n_samples = ncol(normalized)
  )

  tissue_contrasts <- contrasts[contrasts$tissue == tissue, , drop = FALSE]
  if (nrow(tissue_contrasts) == 0L) {
    next
  }
  message("  Computing ", nrow(tissue_contrasts), " contrast(s)")
  indices <- seq_len(nrow(tissue_contrasts))
  worker <- function(index) {
    process_contrast(
      tissue_contrasts[index, , drop = FALSE],
      full_counts,
      full_size_factors,
      cfg$contrast_size_factors,
      output_deseq_dir,
      cfg$force,
      cfg$resume,
      reuse_outputs
    )
  }
  use_forking <- cfg$workers > 1L && length(indices) > 1L &&
    .Platform$OS.type != "windows"
  contrast_records <- if (use_forking) {
    parallel::mclapply(
      indices,
      worker,
      mc.cores = min(cfg$workers, length(indices)),
      mc.preschedule = TRUE
    )
  } else {
    if (cfg$workers > 1L && .Platform$OS.type == "windows") {
      warning("--workers > 1 uses one worker on Windows because fork workers are unavailable")
    }
    lapply(indices, worker)
  }
  failed <- vapply(contrast_records, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(
      "contrast worker failed: ",
      paste(as.character(contrast_records[failed]), collapse = "; "),
      call. = FALSE
    )
  }
  output_records <- c(output_records, contrast_records)
}

metadata$run_status <- "complete"
metadata$completed_at_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
manifest_rows <- lapply(names(metadata), function(name) {
  metadata_record(name, metadata[[name]])
})
manifest_rows[[length(manifest_rows) + 1L]] <- input_manifest
for (record in output_records) {
  kind <- if (nzchar(record$null_type)) "deseq2_results" else "normalized_counts"
  manifest_rows[[length(manifest_rows) + 1L]] <- file_record(
    "output",
    kind,
    value = sprintf(
      "design=%s;n_genes=%d;n_samples=%d",
      record$design,
      record$n_genes,
      record$n_samples
    ),
    tissue = record$tissue,
    null_type = record$null_type,
    rep = record$rep,
    path = record$path,
    include_md5 = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
atomic_write_table(manifest, manifest_path, force = TRUE, gzip = FALSE)
message("\nReference generation complete")
message("Manifest: ", canonical_path(manifest_path, TRUE))
