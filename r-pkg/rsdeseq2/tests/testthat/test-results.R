testthat::test_that("R wrapper results are explicit future work", {
  testthat::expect_error(results(NULL))
})

testthat::test_that("R wrapper assembles primitive DESeq2-shaped result tables", {
  out <- resultsTableRust(
    baseMean = c(g1 = 10, g2 = 20, g3 = 30),
    log2FoldChange = c(1, -0.5, NA),
    lfcSE = c(0.2, 0.3, NA),
    stat = c(5, -1.5, NA),
    pvalue = c(0.001, 0.05, NA),
    dispersion = c(0.1, 0.2, 0.3),
    converged = c(TRUE, TRUE, FALSE)
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(
    colnames(out),
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged")
  )
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.001, 0.05, NA), method = "BH"))
})

testthat::test_that("R wrapper result table helper validates result vectors", {
  testthat::expect_error(resultsTableRust(c(1, NA)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1, 2)))
  testthat::expect_error(resultsTableRust(c(1, 2), padj = c(0.1, 0.2)))
  testthat::expect_error(resultsTableRust(c(1, 2), converged = c(TRUE)))

  out <- resultsTableRust(c(1, 2), rowNames = c("a", "b"))
  testthat::expect_equal(rownames(out), c("a", "b"))
  testthat::expect_true(all(is.na(out$pvalue)))
  testthat::expect_true(all(is.na(out$padj)))
})

testthat::test_that("R wrapper applies primitive independent filtering", {
  table <- resultsTableRust(
    baseMean = c(g1 = 0, g2 = 1, g3 = 2),
    stat = c(1, 1, 1),
    pvalue = c(0.01, 0.02, 0.03)
  )

  filtered <- applyIndependentFilteringRust(
    table,
    alpha = 0.5,
    theta = c(1, 1)
  )
  metadata <- attr(filtered, "independentFiltering")

  testthat::expect_equal(filtered$filtered, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(filtered$padj, c(NA, NA, 0.03))
  testthat::expect_true(metadata$enabled)
  testthat::expect_equal(metadata$selected, 1L)
  testthat::expect_equal(metadata$filterTheta, 1)
  testthat::expect_equal(metadata$filterThreshold, 2)
  testthat::expect_equal(metadata$numRej, c(1L, 1L))
})

testthat::test_that("R wrapper can disable primitive independent filtering", {
  table <- resultsTableRust(
    baseMean = c(0, 1, 2),
    pvalue = c(0.01, NA, 0.03)
  )

  filtered <- applyIndependentFilteringRust(table, enabled = FALSE)
  metadata <- attr(filtered, "independentFiltering")

  testthat::expect_equal(filtered$padj, stats::p.adjust(c(0.01, NA, 0.03), method = "BH"))
  testthat::expect_true(all(is.na(filtered$filtered)))
  testthat::expect_false(metadata$enabled)
  testthat::expect_equal(metadata$alpha, 0.1)
})

testthat::test_that("R wrapper independent filtering validates primitive inputs", {
  table <- resultsTableRust(c(1, 2), pvalue = c(0.1, 0.2))

  testthat::expect_error(applyIndependentFilteringRust(table[, c("baseMean", "pvalue")]))
  testthat::expect_error(applyIndependentFilteringRust(table, alpha = 1))
  testthat::expect_error(applyIndependentFilteringRust(table, theta = c(0.5)))
  testthat::expect_error(applyIndependentFilteringRust(table, theta = c(0, 2)))
})

testthat::test_that("R wrapper applies Cook's cutoff masking and BH adjustment", {
  out <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02, g3 = 0.5),
    maxCooks = c(0, 10, NA),
    cooksCutoff = 5
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02, g3 = 0.5),
    maxCooks = c(0, 10, NA),
    cooksCutoff = 5,
    native = TRUE
  )

  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(out$cooksOutlier, c(FALSE, TRUE, NA))
  testthat::expect_equal(out$pvalue, c(0.01, NA, 0.5))
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.01, NA, 0.5), method = "BH"))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper treats NaN Cook's values as missing", {
  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02, 0.03),
    maxCooks = c(NaN, NA_real_, 10),
    cooksCutoff = 5
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02, 0.03),
    maxCooks = c(NaN, NA_real_, 10),
    cooksCutoff = 5,
    native = TRUE
  )

  testthat::expect_equal(out$cooksOutlier, c(NA, NA, TRUE))
  testthat::expect_equal(out$pvalue, c(0.01, 0.02, NA))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper can disable Cook's cutoff masking", {
  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(100, 200),
    cooksCutoff = FALSE
  )

  testthat::expect_equal(out$pvalue, c(0.01, 0.02))
  testthat::expect_true(all(is.na(out$cooksOutlier)))
})

testthat::test_that("R wrapper low-count Cook's heuristic spares eligible rows", {
  counts <- matrix(c(1, 5, 6, 7, 9, 5, 6, 7), nrow = 2, byrow = TRUE)
  cooks <- matrix(c(10, 0.1, 0.2, 0.3, 10, 0.1, 0.2, 0.3), nrow = 2, byrow = TRUE)

  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(10, 10),
    cooksCutoff = 5,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = TRUE
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(10, 10),
    cooksCutoff = 5,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = TRUE,
    native = TRUE
  )

  testthat::expect_equal(out$cooksOutlier, c(FALSE, TRUE))
  testthat::expect_equal(out$pvalue, c(0.01, NA))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper Cook's cutoff helper validates primitive inputs", {
  testthat::expect_error(applyCooksCutoffRust(c(-0.1), c(1), 5))
  testthat::expect_error(applyCooksCutoffRust(c(0.1), c(Inf), 5))
  testthat::expect_error(applyCooksCutoffRust(c(0.1), c(1), TRUE))
  testthat::expect_error(
    applyCooksCutoffRust(
      c(0.1),
      c(1),
      5,
      counts = matrix(1, nrow = 1, ncol = 2),
      cooks = matrix(1, nrow = 1, ncol = 1),
      lowCountHeuristic = TRUE
    )
  )
})
