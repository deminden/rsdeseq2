testthat::test_that("default two-group Wald and LRT fits select DESeq2 dispersion means", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  condition <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  col_data <- data.frame(condition = condition, row.names = colnames(counts))
  design <- stats::model.matrix(~ condition, col_data)
  factors <- stats::setNames(exp(seq(log(0.4), log(2.2), length.out = ncol(counts))), colnames(counts))
  for (test in c("Wald", "LRT")) {
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
    DESeq2::sizeFactors(dds) <- factors
    reference_args <- list(object = dds, fitType = "mean", test = test,
      minReplicatesForReplace = Inf, quiet = TRUE)
    if (test == "LRT") reference_args$reduced <- ~ 1
    reference <- do.call(DESeq2::DESeq, reference_args)
    fit <- rsdeseq2::DESeq(list(counts = counts, design = design), fitType = "mean",
      test = test, sizeFactors = factors,
      reduced = if (test == "LRT") stats::model.matrix(~ 1, col_data) else NULL)
    expected <- DESeq2::results(reference, name = "condition_B_vs_A",
      cooksCutoff = FALSE, independentFiltering = FALSE)
    actual <- rsdeseq2::results(fit, name = "conditionB",
      cooksCutoff = FALSE, independentFiltering = FALSE)
    testthat::expect_equal(unname(fit$dispersions), unname(DESeq2::dispersions(reference)), tolerance = 1e-7)
    testthat::expect_equal(unname(fit$dispFit), unname(S4Vectors::mcols(reference)$dispFit), tolerance = 1e-7)
    for (column in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
      testthat::expect_equal(unname(actual[[column]]), unname(expected[[column]]), tolerance = 1e-6)
    }
  }
})
