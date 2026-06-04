testthat::test_that("R wrapper Wald test packages primitive Wald fit inputs", {
  beta <- matrix(
    c(
      1, 2,
      2, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A"))
  )
  betaSE <- matrix(c(0.1, 0.5, 0.1, 0.5), nrow = 2, byrow = TRUE)

  fit <- nbinomWaldTestRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = betaSE
  )
  out <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)

  testthat::expect_s3_class(fit, "rsdeseq2PrimitiveWaldFit")
  testthat::expect_equal(out$log2FoldChange, c(2, 4))
  testthat::expect_equal(out$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(attr(out, "comparison"), "coefficient condition_B_vs_A")
})

testthat::test_that("R wrapper Wald test accepts list-like primitive objects", {
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = matrix(
      c(
        1, 2, 3,
        2, 4, 7
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
    )
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object$betaCovariance <- cov

  fit <- nbinomWaldTestRust(object = object)
  out <- results(
    fit,
    contrast = c("condition", "C", "B"),
    reference = "A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(out$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(attr(out, "resultName"), "condition_C_vs_B")
})

testthat::test_that("R wrapper Wald test validates primitive inputs", {
  testthat::expect_error(nbinomWaldTestRust())
  testthat::expect_error(nbinomWaldTestRust(baseMean = 1, beta = matrix(1, nrow = 2)))
})
