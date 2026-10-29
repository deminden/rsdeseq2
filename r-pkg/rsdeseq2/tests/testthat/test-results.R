# Result-contract coverage grouped by responsibility.

testthat::test_that("results rejects objects that were not fitted by a supported route", {
  testthat::expect_error(results(NULL))
})

testthat::test_that("R results computes primitive Wald coefficient and contrast tables", {
  beta <- matrix(
    c(
      1, 2, 5,
      2, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "batch_Y_vs_X"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  cov[, 2, 3] <- 0.1
  cov[, 3, 2] <- 0.1
  fit <- waldFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
    betaSE = matrix(
      c(
        0.1, 0.5, 1,
        0.2, 0.5, 1
      ),
      nrow = 2,
      byrow = TRUE
    )
  )

  named <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)
  contrast <- results(
    fit,
    contrast = list("batch_Y_vs_X", "condition_B_vs_A"),
    independentFiltering = FALSE
  )
  numeric <- results(
    fit,
    contrast = c(0, -1, 1),
    independentFiltering = FALSE
  )

  testthat::expect_equal(named$log2FoldChange, c(2, 4))
  testthat::expect_equal(named$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(attr(named, "resultName"), "condition_B_vs_A")
  testthat::expect_equal(attr(named, "comparison"), "coefficient condition_B_vs_A")
  testthat::expect_equal(attr(named, "contrastAllZero"), "none")
  testthat::expect_equal(contrast$log2FoldChange, c(3, 3))
  testthat::expect_equal(contrast$lfcSE, rep(sqrt(1 + 0.25 - 2 * 0.1), 2))
  testthat::expect_equal(numeric$log2FoldChange, contrast$log2FoldChange)
  testthat::expect_equal(numeric$lfcSE, contrast$lfcSE)
})

testthat::test_that("R primitive Wald fit results validate required contrast diagnostics", {
  beta <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(NULL, c("Intercept", "condition_B_vs_A")))
  fit <- waldFitRust(c(1, 2), beta = beta, betaSE = matrix(c(0.1, 0.2, 0.1, 0.2), nrow = 2))

  contrastWins <- results(
    fit,
    name = "Intercept",
    contrast = c(0, 1),
    independentFiltering = FALSE
  )
  testthat::expect_equal(contrastWins$log2FoldChange, c(3, 4))
  testthat::expect_equal(attr(contrastWins, "resultName"), "contrast")
  testthat::expect_error(results(fit, contrast = c(1, -1), independentFiltering = FALSE))
  testthat::expect_error(waldFitRust(c(1), beta = beta))
})

testthat::test_that("R results supports thresholded Wald alternatives for primitive fits", {
  beta <- matrix(
    c(
      0, 2,
      0, -2,
      0, 0.2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )
  betaSE <- matrix(
    c(
      0.1, 0.5,
      0.1, 0.5,
      0.1, 0.5
    ),
    nrow = 3,
    byrow = TRUE
  )
  fit <- waldFitRust(c(10, 20, 30), beta = beta, betaSE = betaSE)

  greaterAbs <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greaterAbs",
    independentFiltering = FALSE
  )
  greaterAbs2014 <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greaterAbs2014",
    independentFiltering = FALSE
  )
  lessAbs <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "lessAbs",
    independentFiltering = FALSE
  )
  greater <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greater",
    independentFiltering = FALSE
  )
  less <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "less",
    independentFiltering = FALSE
  )

  testthat::expect_equal(greaterAbs$stat[[1]], 4)
  testthat::expect_equal(greaterAbs$pvalue[[1]], 0.02275013293476686, tolerance = 1e-11)
  testthat::expect_equal(greaterAbs2014$stat[[1]], 2)
  testthat::expect_equal(greaterAbs2014$pvalue[[1]], 0.04550026389635842, tolerance = 1e-11)
  testthat::expect_equal(lessAbs$stat[[3]], 1.6)
  testthat::expect_equal(lessAbs$pvalue[[3]], 0.054799291699557995, tolerance = 1e-11)
  testthat::expect_equal(greater$stat[[1]], 2)
  testthat::expect_equal(greater$pvalue[[1]], 0.02275013194817921, tolerance = 1e-11)
  testthat::expect_equal(less$stat[[2]], -2)
  testthat::expect_equal(less$pvalue[[2]], 0.02275013194817921, tolerance = 1e-11)
  testthat::expect_equal(attr(greater, "lfcThreshold"), 1)
  testthat::expect_equal(attr(greater, "altHypothesis"), "greater")
})

testthat::test_that("R results supports thresholded contrast tests with t tails", {
  beta <- matrix(
    c(0, 2, 4),
    nrow = 1,
    dimnames = list("g1", c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(1, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  fit <- waldFitRust(c(g1 = 10), beta = beta, betaCovariance = cov)

  out <- results(
    fit,
    contrast = list("condition_C_vs_A", "condition_B_vs_A"),
    lfcThreshold = 1,
    altHypothesis = "greater",
    useT = TRUE,
    degreesOfFreedom = 10,
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 2)
  testthat::expect_equal(out$lfcSE[[1]], sqrt(0.5))
  testthat::expect_equal(out$stat[[1]], (2 - 1) / sqrt(0.5))
  testthat::expect_equal(
    out$pvalue[[1]],
    stats::pt((2 - 1) / sqrt(0.5), df = 10, lower.tail = FALSE),
    tolerance = 1e-15
  )
  testthat::expect_equal(attr(out, "degreesOfFreedom"), 10)
})

testthat::test_that("R results validates thresholded Wald options", {
  beta <- matrix(c(1, 2), nrow = 1, dimnames = list("g1", c("Intercept", "condition_B_vs_A")))
  betaSE <- matrix(c(0.1, 0.5), nrow = 1)
  fit <- waldFitRust(c(g1 = 10), beta = beta, betaSE = betaSE)

  testthat::expect_error(results(fit, lfcThreshold = -1))
  testthat::expect_error(results(fit, lfcThreshold = 0, altHypothesis = "lessAbs"))
  testthat::expect_error(results(fit, useT = "yes", degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = c(TRUE, FALSE), degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = NA, degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = TRUE))
  testthat::expect_error(results(fit, useT = TRUE, degreesOfFreedom = 10, altHypothesis = "greaterAbsUPSHOT"))
})

testthat::test_that("R wrapper assembles primitive DESeq2-shaped result tables", {
  out <- resultsTableRust(
    baseMean = c(g1 = 10, g2 = 20, g3 = 30),
    log2FoldChange = c(g3 = NA, g1 = 1, g2 = -0.5),
    lfcSE = c(g2 = 0.3, g3 = NA, g1 = 0.2),
    stat = c(g1 = 5, g3 = NA, g2 = -1.5),
    pvalue = c(g2 = 0.05, g1 = 0.001, g3 = NA),
    dispersion = c(g3 = 0.3, g2 = 0.2, g1 = 0.1),
    converged = c(g2 = TRUE, g1 = TRUE, g3 = FALSE),
    contrast = c(Intercept = 0, condition_B_vs_A = 1)
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(
    colnames(out),
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged")
  )
  testthat::expect_equal(out$log2FoldChange, c(1, -0.5, NA))
  testthat::expect_equal(out$lfcSE, c(0.2, 0.3, NA))
  testthat::expect_equal(out$stat, c(5, -1.5, NA))
  testthat::expect_equal(out$dispersion, c(0.1, 0.2, 0.3))
  testthat::expect_equal(out$converged, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(attr(out, "contrast"), c(Intercept = 0, condition_B_vs_A = 1))
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.001, 0.05, NA), method = "BH"))
})

testthat::test_that("R wrapper result table helper validates result vectors", {
  testthat::expect_error(resultsTableRust(c(1, NA)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1, 2)))
  testthat::expect_error(resultsTableRust(c(1, 2), padj = c(0.1, 0.2)))
  testthat::expect_error(resultsTableRust(c(1, 2), converged = c(TRUE)))
  testthat::expect_error(resultsTableRust(c(1, 2), contrast = numeric()))
  testthat::expect_error(resultsTableRust(c(1, 2), contrast = c(0, NA)))
  testthat::expect_error(
    resultsTableRust(
      c(g1 = 1, g2 = 2),
      pvalue = c(g1 = 0.1, gx = 0.2)
    ),
    "pvalue names must contain all beta row names"
  )

  out <- resultsTableRust(c(1, 2), rowNames = c("a", "b"))
  testthat::expect_equal(rownames(out), c("a", "b"))
  testthat::expect_true(all(is.na(out$pvalue)))
  testthat::expect_true(all(is.na(out$padj)))
})
