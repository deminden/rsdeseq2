testthat::test_that("borrowed Rust projection matches independent matrix arithmetic", {
  beta <- matrix(c(1, 3, -2, 4, 5, 2), 3, 2, byrow = TRUE)
  covariance <- array(0, c(3, 2, 2))
  covariance[1, , ] <- matrix(c(4, .5, .25, 9), 2, byrow = TRUE)
  covariance[2, , ] <- matrix(c(1, -.5, .25, 4), 2, byrow = TRUE)
  covariance[3, , ] <- matrix(c(9, 1, .5, 16), 2, byrow = TRUE)
  fit <- list(beta = beta, betaSE = NULL, betaCovariance = covariance)
  for (contrast in list(c(-1, 1), c(.25, -2), c(0, -1))) {
    effect <- as.numeric(beta %*% contrast)
    se <- sqrt(vapply(1:3, function(gene)
      as.numeric(t(contrast) %*% covariance[gene, , ] %*% contrast), numeric(1)))
    raw <- rsdeseq2:::.rsdeseq2_lrt_from_primitive_fit(fit, contrast)
    testthat::expect_equal(raw$log2FoldChange, effect, tolerance = 1e-14)
    testthat::expect_equal(raw$lfcSE, se, tolerance = 1e-14)
    options <- list(lfcThreshold = .75, altHypothesis = "greater", useT = TRUE,
      degreesOfFreedom = c(2, 10, 100))
    actual <- rsdeseq2:::.rsdeseq2_wald_from_primitive_fit(fit, contrast, options)
    statistic <- pmax((effect - .75) / se, 0)
    testthat::expect_equal(actual$stat, statistic, tolerance = 1e-14)
    testthat::expect_equal(actual$pvalue,
      stats::pt((effect - .75) / se, c(2, 10, 100), lower.tail = FALSE), tolerance = 1e-13)
  }
})

testthat::test_that("Wald contrast invalid rows follow standalone policy while raw LRT projection remains available", {
  fit <- list(beta = cbind(1, 2:5), betaSE = cbind(1, c(.5, 0, NaN, 2)),
    betaCovariance = array(NaN, c(4, 2, 2)))
  raw <- rsdeseq2:::.rsdeseq2_lrt_from_primitive_fit(fit, c(0, -2))
  testthat::expect_identical(raw$log2FoldChange, c(-4, -6, -8, -10))
  testthat::expect_equal(raw$lfcSE, c(1, 0, NA_real_, 4), tolerance = 0)
  options <- list(lfcThreshold = 0, altHypothesis = "greaterAbs", useT = TRUE,
    degreesOfFreedom = c(10, 10, 10, NA_real_))
  wald <- rsdeseq2:::.rsdeseq2_wald_from_primitive_fit(fit, c(0, -2), options)
  testthat::expect_identical(wald$log2FoldChange, c(-4, NA_real_, NA_real_, -10))
  testthat::expect_identical(wald$lfcSE, c(1, NA_real_, NA_real_, 4))
  testthat::expect_identical(wald$stat, c(-4, NA_real_, NA_real_, -2.5))
  testthat::expect_identical(is.na(wald$pvalue), c(FALSE, TRUE, TRUE, TRUE))
})
