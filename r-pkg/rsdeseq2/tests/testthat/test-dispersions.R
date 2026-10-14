testthat::test_that("R dispersion entry point returns the Rust fit dispersions", {
  counts <- .rsdeseq2_test_counts()
  design <- stats::model.matrix(~ factor(rep(c("A", "B"), each = 4)))

  fit <- DESeq(list(counts = counts, design = design), fitType = "mean")
  dispersions <- estimateDispersionsRust(counts, design, fitType = "mean")

  testthat::expect_named(dispersions, rownames(counts))
  testthat::expect_equal(unname(dispersions), unname(fit$results$dispersion))
  testthat::expect_true(all(is.finite(dispersions) & dispersions >= 0))
})
