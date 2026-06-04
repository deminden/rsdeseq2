testthat::test_that("R wrapper LRT packages primitive LRT fit inputs", {
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

  fit <- nbinomLRTRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = betaSE,
    lrtStat = c(6, 8),
    lrtDf = 1
  )
  out <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)

  testthat::expect_s3_class(fit, "rsdeseq2PrimitiveLrtFit")
  testthat::expect_equal(out$log2FoldChange, c(2, 4))
  testthat::expect_equal(out$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(out$stat, c(6, 8))
  testthat::expect_equal(out$pvalue, stats::pchisq(c(6, 8), df = 1, lower.tail = FALSE))
  testthat::expect_equal(attr(out, "comparison"), "coefficient condition_B_vs_A; LRT full vs reduced")
})

testthat::test_that("R wrapper LRT results reports selected contrast LFC with LRT p-values", {
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
    ),
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01)
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object$betaCovariance <- cov

  fit <- nbinomLRTRust(object = object)
  out <- results(
    fit,
    contrast = c("condition", "C", "B"),
    reference = "A",
    independentFiltering = FALSE
  )
  contrastWins <- results(
    fit,
    name = "condition_B_vs_A",
    contrast = c("condition", "C", "B"),
    reference = "A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(out$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(out$stat, c(6, 8))
  testthat::expect_equal(out$pvalue, c(0.02, 0.01))
  testthat::expect_equal(attr(out, "resultName"), "condition_C_vs_B")
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(contrastWins$log2FoldChange, out$log2FoldChange)
  testthat::expect_equal(attr(contrastWins, "resultName"), "condition_C_vs_B")
})

testthat::test_that("R wrapper LRT applies contrast all-zero only to reported fold change", {
  beta <- matrix(
    c(
      1, 2, 3,
      1, 2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  counts <- matrix(
    c(
      0, 0, 0, 0, 50, 60,
      10, 12, 30, 36, 50, 60
    ),
    nrow = 2,
    byrow = TRUE
  )
  modelMatrix <- matrix(
    c(
      1, 0, 0,
      1, 0, 0,
      1, 1, 0,
      1, 1, 0,
      1, 0, 1,
      1, 0, 1
    ),
    nrow = 6,
    byrow = TRUE
  )
  fit <- lrtFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts,
    sampleLevels = c("A", "A", "B", "B", "C", "C"),
    modelMatrix = modelMatrix
  )

  character <- results(
    fit,
    contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )
  numeric <- results(
    fit,
    contrast = c(0, 1, 0),
    independentFiltering = FALSE
  )

  testthat::expect_equal(character$log2FoldChange[[1]], 0)
  testthat::expect_equal(character$stat[[1]], 6)
  testthat::expect_equal(character$pvalue[[1]], 0.02)
  testthat::expect_equal(numeric$log2FoldChange[[1]], 2)
  testthat::expect_equal(attr(character, "contrastAllZero"), "character")
  testthat::expect_equal(attr(numeric, "contrastAllZero"), "numeric")
})

testthat::test_that("R wrapper LRT validates primitive inputs", {
  testthat::expect_error(nbinomLRTRust())
  testthat::expect_error(nbinomLRTRust(baseMean = 1, beta = matrix(1, nrow = 2), lrtStat = 1))
  testthat::expect_error(lrtFitRust(baseMean = 1, beta = matrix(1, nrow = 1), lrtStat = 1))
})

testthat::test_that("R wrapper LRT validates character contrast sample levels", {
  beta <- matrix(
    c(
      0, 2,
      0, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A"))
  )
  fit <- lrtFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(c(0.1, 0.5, 0.1, 0.5), nrow = 2, byrow = TRUE),
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = matrix(c(0, 0, 10, 12), nrow = 2, byrow = TRUE),
    sampleLevels = c("A", "A")
  )

  testthat::expect_error(
    results(fit, contrast = c("condition", "B", "A"), independentFiltering = FALSE),
    "sampleLevels to contain numerator and denominator"
  )
})
