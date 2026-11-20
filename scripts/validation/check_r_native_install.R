## Exercise an installed package without optional Bioconductor packages.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
library(rsdeseq2, lib.loc = normalizePath(args[[1L]], mustWork = TRUE))

dll <- getLoadedDLLs()[["rsdeseq2"]]
stopifnot(!is.null(dll), identical(dll[["dynamicLookup"]], FALSE))
routines <- getDLLRegisteredRoutines(dll)$.Call
signatures <- c(
  rsdeseq2_validate_counts = 1L,
  rsdeseq2_estimate_size_factors = 4L,
  rsdeseq2_normalized_counts = 3L,
  rsdeseq2_base_metadata = 4L,
  rsdeseq2_fit = 11L,
  rsdeseq2_transform = 16L
)
for (routine in names(signatures)) {
  stopifnot(identical(routines[[routine]]$numParameters, signatures[[routine]]))
}

counts <- matrix(c(10, 20, 0, 20, 40, 5), nrow = 3L)
factors <- estimateSizeFactorsRust(counts)
stopifnot(isTRUE(all.equal(factors, c(1 / sqrt(2), sqrt(2)), tolerance = 1e-12)))
normalized <- normalizedCountsRust(counts, c(1, 2))
metadata <- baseMetadataRust(counts, c(1, 2))
stopifnot(
  identical(normalized, matrix(c(10, 20, 0, 10, 20, 2.5), nrow = 3L)),
  identical(metadata$baseMean, c(10, 20, 1.25)),
  identical(metadata$allZero, rep(FALSE, 3L))
)

rejected <- tryCatch(.Call(routines$rsdeseq2_validate_counts, counts + 0.5), error = identity)
stopifnot(
  inherits(rejected, "error"),
  grepl("valid u32 count", conditionMessage(rejected), fixed = TRUE),
  identical(.Call(routines$rsdeseq2_validate_counts, counts), counts)
)

set.seed(42)
counts <- matrix(as.double(stats::rnbinom(320L, mu = 50, size = 5)), 40L, 8L)
design <- cbind(1, rep(c(0, 1), each = 4L))
fit <- function(workers) {
  .Call(routines$rsdeseq2_fit, counts, design, NULL, rep(1, 8L), NULL, NULL,
    2L, 0L, 0L, 10, workers)
}
single <- fit(1L)
parallel <- fit(2L)
stopifnot(
  identical(dim(parallel$mu), dim(counts)),
  all(is.finite(parallel$mu)),
  identical(parallel$beta, single$beta),
  identical(parallel$mu, single$mu)
)
# Use per-cell factors and an interior zero row to exercise row selection in
# both estimated and frozen rlog. A separate fit sees only the active rows.
zero <- 13L
counts[zero, ] <- 0
active <- seq_len(nrow(counts))[-zero]
cell_factors <- outer(seq_len(nrow(counts)), seq_len(ncol(counts)),
  function(gene, sample) 0.8 + (gene + 3 * sample) %% 11 / 20)
frozen_intercepts <- 4 + seq_len(nrow(counts)) / 100
rlog_fit <- function(rows, workers, intercept = NULL, prior = NULL) {
  disp_fit <- if (is.null(intercept)) rep(0.12, length(rows)) else NULL
  .Call(routines$rsdeseq2_transform, counts[rows, , drop = FALSE], design, NULL,
    cell_factors[rows, , drop = FALSE], NULL, 2L, 1L, FALSE, 0L, workers,
    intercept, prior, -1L, NULL, disp_fit, NULL)
}

masked <- rlog_fit(seq_len(nrow(counts)), 2L)
compact <- rlog_fit(active, 1L)
stopifnot(
  all(is.finite(masked$transformed)),
  identical(masked$transformed[active, ], compact$transformed),
  identical(masked$transformed[zero, ], rep(0, ncol(counts))),
  identical(masked$intercept[active], compact$intercept),
  identical(masked$intercept[zero], 0),
  is.finite(masked$samplePriorVariance),
  masked$samplePriorVariance > 0,
  identical(masked$samplePriorVariance, compact$samplePriorVariance)
)

masked <- rlog_fit(seq_len(nrow(counts)), 2L, frozen_intercepts, 0.5)
compact <- rlog_fit(active, 1L, frozen_intercepts[active], 0.5)
stopifnot(
  all(is.finite(masked$transformed)),
  identical(masked$transformed[active, ], compact$transformed),
  identical(masked$transformed[zero, ], rep(frozen_intercepts[zero], ncol(counts))),
  identical(masked$intercept, frozen_intercepts),
  identical(masked$samplePriorVariance, 0.5)
)
cat("Installed R native routines passed on", R.version$platform, "\n")
