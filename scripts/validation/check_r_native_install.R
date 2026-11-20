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
  rsdeseq2_fit = 11L
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
cat("Installed R native routines passed on", R.version$platform, "\n")
