test_that("the standalone BH engine matches R on names, ties, tails and missingness", {
  set.seed(437)
  p <- c(0, 1, NA_real_, NaN, 1e-320, 1e-200, runif(2002), rep(.01, 12))
  names(p) <- paste0("g", seq_along(p))
  expect_identical(rsdeseq2:::.rsdeseq2_bh_adjust(p), stats::p.adjust(p, "BH"))
  for (p in list(numeric(), c(NA_real_, NA_real_), .1, c(NA_real_, .1)))
    expect_identical(rsdeseq2:::.rsdeseq2_bh_adjust(p), stats::p.adjust(p, "BH"))
})

test_that("standalone filtering preserves native DESeq2 adjustment and threshold decisions", {
  set.seed(438)
  n <- 2002L
  filter <- c(rep(0, 102), seq_len(n - 102L))
  p <- runif(n)^ifelse(seq_len(n) > 800, 5, 1)
  p[c(3, 120, 1750)] <- NA_real_
  table <- resultsTableRust(filter, pvalue = p)
  referenceInput <- DESeq2::DESeqResults(S4Vectors::DataFrame(table))
  for (alpha in c(.01, .1)) for (theta in list(NULL, c(0, .15, .3, .5, .75, .95), c(1, 1))) {
    args <- list(res = referenceInput, independentFiltering = TRUE, alpha = alpha, pAdjustMethod = "BH")
    if (!is.null(theta)) args$theta <- theta
    expected <- do.call(getFromNamespace("pvalueAdjustment", "DESeq2"), args)
    actual <- applyIndependentFilteringRust(table, alpha = alpha, theta = theta)
    meta <- attr(actual, "independentFiltering"); ref <- S4Vectors::metadata(expected)
    expect_identical(actual$padj, as.numeric(expected$padj))
    expect_equal(as.numeric(meta$numRej), ref$filterNumRej$numRej, tolerance = 0)
    expect_identical(meta$theta, ref$filterNumRej$theta)
    expect_identical(meta$filterTheta, unname(ref$filterTheta))
    expect_identical(meta$filterThreshold, unname(ref$filterThreshold))
    # Same absolute bound as the existing standalone R-lowess reference tests.
    expect_lte(max(abs(meta$lo.fit$y - ref$lo.fit$y)), 1e-9)
    if (is.null(theta)) expect_gt(meta$selected, 1L)
  }
})
