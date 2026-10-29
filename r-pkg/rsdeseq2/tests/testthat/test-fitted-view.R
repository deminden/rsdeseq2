test_that("fitted views map gene order and reject mutations of fitted inputs", {
  ct <- .rsdeseq2_test_counts()
  cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)), row.names = colnames(ct))
  fit <- DESeq(DESeq2::DESeqDataSetFromMatrix(ct, cd, ~ condition),
               fitType = "mean", minReplicatesForReplace = Inf)
  full <- as.data.frame(results(fit, independentFiltering = FALSE))
  subset <- as.data.frame(results(fit[c(3, 1), ], independentFiltering = FALSE))
  expect_identical(rownames(subset), rownames(ct)[c(3, 1)])
  expect_equal(subset$pvalue, full$pvalue[c(3, 1)], tolerance = 0)
  expect_equal(subset$padj, p.adjust(subset$pvalue, "BH"), tolerance = 0)
  changed <- fit
  SummarizedExperiment::assay(changed, "counts")[1, 1] <- ct[1, 1] + 1L
  expect_error(results(changed), "counts changed")
  expect_error(resultsNames(changed), "counts changed")
  expect_error(lfcShrink(changed, coef = 2, type = "normal"), "counts changed")
  expect_error(results(fit[, 8:1]), "samples changed")
  changed <- fit
  DESeq2::sizeFactors(changed) <- DESeq2::sizeFactors(changed) * 2
  expect_error(results(changed), "size factors changed")
  changed <- fit
  SummarizedExperiment::colData(changed)$condition <- rev(cd$condition)
  expect_error(results(changed), "design changed")
})
