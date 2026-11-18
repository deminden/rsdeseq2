# Run the package's testthat suite against an installed rsdeseq2 and fail if
# any test is skipped or fails.
#
# In the reference lane every optional package (DESeq2, apeglm, ashr, ...) is
# installed, so a skip there means a parity test silently did not run.
#
# Usage:
#   Rscript --vanilla scripts/validation/check_no_skipped_r_tests.R LIBRARY TESTDIR
# LIBRARY is a library containing rsdeseq2 (for example rsdeseq2.Rcheck);
# TESTDIR is r-pkg/rsdeseq2/tests/testthat.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
.libPaths(c(normalizePath(args[[1]]), .libPaths()))
# Some tests start a child R process; it must see the same library.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
suppressPackageStartupMessages({
  library(testthat)
  library(rsdeseq2)
})
results <- as.data.frame(
  test_dir(
    args[[2]],
    package = "rsdeseq2",
    load_package = "installed",
    reporter = "summary",
    stop_on_failure = FALSE
  )
)
skipped <- results[results$skipped, c("file", "test"), drop = FALSE]
failed <- results[results$failed > 0 | results$error, c("file", "test"), drop = FALSE]
cat(sprintf(
  "\n%d tests, %d expectations, %d skipped, %d failed\n",
  nrow(results), sum(results$nb), nrow(skipped), nrow(failed)
))
if (nrow(skipped) > 0L) {
  cat("Skipped tests (a parity check did not run):\n")
  print(skipped, row.names = FALSE)
}
if (nrow(failed) > 0L) {
  cat("Failed tests:\n")
  print(failed, row.names = FALSE)
}
if (nrow(skipped) > 0L || nrow(failed) > 0L) {
  quit(status = 1L)
}
