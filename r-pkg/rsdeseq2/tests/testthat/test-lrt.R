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

  relevelFit <- lrtFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts,
    sampleLevels = c("A", "A", "B", "B", "C", "C"),
    factorLevels = list(condition = c("B", "A", "C")),
    factorReferences = c(condition = "A"),
    modelMatrix = modelMatrix
  )
  relevel <- results(
    relevelFit,
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(relevel$log2FoldChange, c(1, 1))
  testthat::expect_equal(attr(relevel, "resultName"), "condition_C_vs_B")
  testthat::expect_equal(attr(relevel, "contrast")[["condition_B_vs_A"]], -1)
  testthat::expect_equal(attr(relevel, "contrast")[["condition_C_vs_A"]], 1)
})

testthat::test_that("R wrapper LRT aligns colData row names for character contrast all-zero handling", {
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
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  fit <- lrtFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts,
    colData = data.frame(
      condition = factor(c("C", "A", "B", "C", "A", "B"), levels = c("A", "B", "C")),
      row.names = c("s5", "s1", "s3", "s6", "s2", "s4")
    )
  )

  out <- results(
    fit,
    contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$stat[[1]], 6)
  testthat::expect_equal(out$pvalue[[1]], 0.02)
  testthat::expect_equal(out$log2FoldChange[[2]], 2)
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")

  aliasFit <- lrtFitRust(
    baseMean = c(10, 20),
    beta = `colnames<-`(beta, c("Intercept", "cell.type_B_vs_A", "cell.type_C_vs_A")),
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts,
    colData = data.frame(
      "cell type" = factor(c("C", "A", "B", "C", "A", "B"), levels = c("A", "B", "C")),
      check.names = FALSE,
      row.names = c("s5", "s1", "s3", "s6", "s2", "s4")
    )
  )
  cleaned <- results(
    aliasFit,
    contrast = c("cell.type", "B", "A"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(cleaned$log2FoldChange[[1]], 0)
  testthat::expect_equal(cleaned$stat[[1]], 6)
  testthat::expect_equal(cleaned$pvalue[[1]], 0.02)
  testthat::expect_equal(attr(cleaned, "contrastAllZero"), "character")
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
    "sampleLevels or colData to contain numerator level B"
  )
})
