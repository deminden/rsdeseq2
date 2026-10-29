# Count, sample, and fitted-state alignment for result queries.

testthat::test_that("R wrapper aligns named model matrix columns with cleaned result aliases", {
  beta <- matrix(
    c(
      1, 2, 0,
      1, 2, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_as.double", "batch Y"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  counts <- matrix(
    c(
      50, 60, 0, 0,
      10, 12, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3", "s4"))
  )
  modelMatrix <- matrix(
    c(
      1, 1, 0,
      0, 1, 0,
      1, 1, 0,
      0, 1, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("s3", "s1", "s4", "s2"), c("cell.type_as.double", "Intercept", "batch.Y"))
  )
  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    counts = counts,
    modelMatrix = modelMatrix
  )

  out <- results(
    fit,
    contrast = c(0, 1, -1),
    independentFiltering = FALSE
  )

  testthat::expect_equal(colnames(fit$modelMatrix), c("Intercept", "cell type_as.double", "batch Y"))
  testthat::expect_equal(rownames(fit$modelMatrix), colnames(counts))
  testthat::expect_equal(unname(fit$modelMatrix[, "cell type_as.double"]), c(0, 0, 1, 1))
  testthat::expect_equal(out$log2FoldChange, c(0, 2))
  testthat::expect_equal(out$stat[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
  testthat::expect_equal(attr(out, "contrastAllZero"), "numeric")

  duplicateColumns <- modelMatrix
  colnames(duplicateColumns) <- c("Intercept", "cell.type_as.double", "cell.type_as.double")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      counts = counts,
      modelMatrix = duplicateColumns
    ),
    "modelMatrix column names must be unique when supplied"
  )
})

testthat::test_that("R wrapper aligns named betaSE columns with cleaned result aliases", {
  beta <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_as.double", "batch Y"))
  )
  betaSE <- matrix(
    c(
      0.31, 0.11, 0.21,
      0.32, 0.12, 0.22
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("batch.Y", "Intercept", "cell.type_as.double"))
  )

  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = betaSE
  )

  testthat::expect_equal(colnames(fit$betaSE), colnames(beta))
  testthat::expect_equal(unname(fit$betaSE[, "Intercept"]), c(0.11, 0.12))
  testthat::expect_equal(unname(fit$betaSE[, "cell type_as.double"]), c(0.21, 0.22))
  testthat::expect_equal(unname(fit$betaSE[, "batch Y"]), c(0.31, 0.32))

  duplicateColumns <- betaSE
  colnames(duplicateColumns) <- c("Intercept", "cell.type_as.double", "cell.type_as.double")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaSE = duplicateColumns
    ),
    "betaSE column names must be unique when supplied"
  )
})

testthat::test_that("R wrapper aligns named betaCovariance coefficient axes with cleaned result aliases", {
  beta <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_as.double", "batch Y"))
  )
  axisNames <- c("batch.Y", "Intercept", "cell.type_as.double")
  betaCovariance <- array(
    0,
    dim = c(2, 3, 3),
    dimnames = list(c("g1", "g2"), axisNames, axisNames)
  )
  betaCovariance[, "Intercept", "Intercept"] <- c(0.11, 0.12)
  betaCovariance[, "cell.type_as.double", "cell.type_as.double"] <- c(0.21, 0.22)
  betaCovariance[, "batch.Y", "batch.Y"] <- c(0.31, 0.32)
  betaCovariance[, "cell.type_as.double", "batch.Y"] <- c(0.04, 0.05)
  betaCovariance[, "batch.Y", "cell.type_as.double"] <- c(0.04, 0.05)

  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = betaCovariance
  )

  testthat::expect_equal(dimnames(fit$betaCovariance)[[2L]], colnames(beta))
  testthat::expect_equal(dimnames(fit$betaCovariance)[[3L]], colnames(beta))
  testthat::expect_equal(
    unname(fit$betaCovariance[, "cell type_as.double", "batch Y"]),
    c(0.04, 0.05)
  )
  testthat::expect_equal(
    unname(fit$betaCovariance[, "batch Y", "batch Y"]),
    c(0.31, 0.32)
  )

  duplicateAxis <- betaCovariance
  dimnames(duplicateAxis)[[2L]] <- c("Intercept", "cell.type_as.double", "cell.type_as.double")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = duplicateAxis
    ),
    "betaCovariance coefficient names must be unique when supplied"
  )
})

testthat::test_that("R results coerces DataFrame-like colData before sample alignment", {
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
      10, 12, 20, 24, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    counts = counts,
    colData = matrix(
      c("C", "A", "B", "C", "A", "B"),
      ncol = 1L,
      dimnames = list(
        c("s5", "s1", "s3", "s6", "s2", "s4"),
        "condition"
      )
    )
  )

  fit <- nbinomWaldTestRust(object = object)
  out <- results(fit, contrast = c("condition", "B", "A"), independentFiltering = FALSE)

  testthat::expect_true(is.data.frame(fit$colData))
  testthat::expect_equal(rownames(fit$colData), colnames(counts))
  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")

  bad <- object
  bad$colData <- environment()
  testthat::expect_error(
    nbinomWaldTestRust(object = bad),
    "colData must be a data.frame or coercible with as.data.frame"
  )
})

testthat::test_that("R results aligns named sampleLevels to count columns", {
  beta <- matrix(
    c(
      1, 2, 3,
      1, 2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  counts <- matrix(
    c(
      0, 0, 0, 0, 10, 10,
      0, 5, 0, 8, 10, 10
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  sampleLevels <- c(s5 = "C", s3 = "B", s1 = "A", s6 = "C", s4 = "B", s2 = "A")
  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(1, nrow = 2, ncol = 3),
    counts = counts,
    sampleLevels = sampleLevels
  )

  out <- results(fit, contrast = c("condition", "B", "A"), independentFiltering = FALSE)

  testthat::expect_equal(unname(fit$sampleLevels), c("A", "A", "B", "B", "C", "C"))
  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      counts = counts,
      sampleLevels = c(s1 = "A", s2 = "A", s3 = "B", s4 = "B", s5 = "C", s5 = "C")
    ),
    "sampleLevels names must be unique"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      counts = counts,
      sampleLevels = c(s1 = "A", s2 = "A", s3 = "B", s4 = "B", s5 = "C", sx = "C")
    ),
    "sampleLevels names must contain all count sample names"
  )
})

testthat::test_that("R primitive fits align named baseMean to beta rows", {
  beta <- matrix(
    c(
      1, 2,
      1, 3,
      1, 4
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )
  fit <- waldFitRust(
    baseMean = c(g3 = 30, g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(1, nrow = 3, ncol = 2)
  )

  out <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)

  testthat::expect_equal(unname(fit$baseMean), c(10, 20, 30))
  testthat::expect_equal(names(fit$baseMean), c("g1", "g2", "g3"))
  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(out$baseMean, c(10, 20, 30))
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20, g2 = 30),
      beta = beta
    ),
    "baseMean names must be unique"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20, gx = 30),
      beta = beta
    ),
    "baseMean names must contain all beta row names"
  )
})

testthat::test_that("R primitive fits align named row-shaped inputs to beta rows", {
  beta <- matrix(
    c(
      0, 2,
      0, 4,
      0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )
  betaSE <- matrix(
    c(
      0.1, 0.2,
      0.3, 0.4,
      0.5, 0.6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g3", "g1", "g2"), colnames(beta))
  )
  betaCovariance <- array(
    0,
    dim = c(3, 2, 2),
    dimnames = list(c("g2", "g3", "g1"), colnames(beta), colnames(beta))
  )
  betaCovariance[, "condition_B_vs_A", "condition_B_vs_A"] <- c(0.36, 0.04, 0.16)
  counts <- matrix(
    c(
      0, 0, 10, 10,
      0, 5, 10, 10,
      0, 0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g2", "g3", "g1"), paste0("s", 1:4))
  )

  fit <- waldFitRust(
    baseMean = c(g3 = 30, g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    counts = counts,
    sampleLevels = c("A", "A", "B", "B")
  )
  out <- results(fit, contrast = c("condition", "B", "A"), independentFiltering = FALSE)

  testthat::expect_equal(rownames(fit$betaSE), c("g1", "g2", "g3"))
  testthat::expect_equal(unname(fit$betaSE[, "condition_B_vs_A"]), c(0.4, 0.6, 0.2))
  testthat::expect_equal(dimnames(fit$betaCovariance)[[1L]], c("g1", "g2", "g3"))
  testthat::expect_equal(
    unname(fit$betaCovariance[, "condition_B_vs_A", "condition_B_vs_A"]),
    c(0.16, 0.36, 0.04)
  )
  testthat::expect_equal(rownames(fit$counts), c("g1", "g2", "g3"))
  testthat::expect_equal(out$log2FoldChange, c(2, 4, 6))

  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20, g3 = 30),
      beta = beta,
      betaSE = betaSE[c("g1", "g2", "g2"), , drop = FALSE]
    ),
    "betaSE row names must be unique"
  )
  badCovariance <- betaCovariance
  dimnames(badCovariance)[[1L]] <- c("g1", "g2", "gx")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20, g3 = 30),
      beta = beta,
      betaCovariance = badCovariance
    ),
    "betaCovariance row names must contain all beta row names"
  )
  badCounts <- counts
  rownames(badCounts) <- c("g1", "g2", "gx")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20, g3 = 30),
      beta = beta,
      counts = badCounts
    ),
    "counts row names must contain all beta row names"
  )
})

testthat::test_that("R primitive LRT fits align named test vectors to beta rows", {
  beta <- matrix(
    c(
      0, 2,
      0, 4,
      0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )

  fit <- lrtFitRust(
    baseMean = c(g3 = 30, g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(1, nrow = 3, ncol = 2),
    lrtStat = c(g2 = 8, g3 = 12, g1 = 4),
    lrtDf = c(g3 = 3, g1 = 1, g2 = 2)
  )
  explicit <- lrtFitRust(
    baseMean = c(g3 = 30, g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(1, nrow = 3, ncol = 2),
    lrtStat = c(g2 = 8, g3 = 12, g1 = 4),
    lrtPvalue = c(g3 = 0.003, g1 = 0.04, g2 = 0.008)
  )
  scalarDf <- lrtFitRust(
    baseMean = c(g3 = 30, g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(1, nrow = 3, ncol = 2),
    lrtStat = c(g2 = 8, g3 = 12, g1 = 4),
    lrtDf = c(model = 2)
  )

  testthat::expect_equal(unname(fit$lrtStat), c(4, 8, 12))
  testthat::expect_equal(unname(fit$lrtDf), c(1, 2, 3))
  testthat::expect_identical(names(fit$lrtPvalue), rownames(beta))
  testthat::expect_identical(names(scalarDf$lrtPvalue), rownames(beta))
  testthat::expect_equal(unname(fit$lrtPvalue), stats::pchisq(c(4, 8, 12), df = c(1, 2, 3), lower.tail = FALSE))
  testthat::expect_equal(unname(explicit$lrtPvalue), c(0.04, 0.008, 0.003))
  testthat::expect_equal(unname(scalarDf$lrtDf), c(2, 2, 2))
  testthat::expect_equal(unname(scalarDf$lrtPvalue), stats::pchisq(c(4, 8, 12), df = 2, lower.tail = FALSE))

  testthat::expect_error(
    lrtFitRust(
      baseMean = c(g1 = 10, g2 = 20, g3 = 30),
      beta = beta,
      lrtStat = c(g1 = 4, g2 = 8, g2 = 12),
      lrtDf = 1
    ),
    "lrtStat names must be unique"
  )
  testthat::expect_error(
    lrtFitRust(
      baseMean = c(g1 = 10, g2 = 20, g3 = 30),
      beta = beta,
      lrtStat = c(g1 = 4, g2 = 8, gx = 12),
      lrtDf = 1
    ),
    "lrtStat names must contain all beta row names"
  )
})

testthat::test_that("R object extraction preserves mcols row names for primitive alignment", {
  beta <- matrix(
    c(
      0, 2,
      0, 4,
      0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )
  object <- list(
    beta = beta,
    betaSE = matrix(1, nrow = 3, ncol = 2),
    mcols = data.frame(
      baseMean = c(30, 10, 20),
      stat = c(12, 4, 8),
      pvalue = c(0.003, 0.04, 0.008),
      row.names = c("g3", "g1", "g2")
    )
  )

  wald <- nbinomWaldTestRust(object = object)
  lrt <- nbinomLRTRust(object = object)
  out <- results(lrt, name = "condition_B_vs_A", independentFiltering = FALSE)

  testthat::expect_equal(unname(wald$baseMean), c(10, 20, 30))
  testthat::expect_equal(names(wald$baseMean), c("g1", "g2", "g3"))
  testthat::expect_equal(unname(lrt$lrtStat), c(4, 8, 12))
  testthat::expect_equal(unname(lrt$lrtPvalue), c(0.04, 0.008, 0.003))
  testthat::expect_equal(out$baseMean, c(10, 20, 30))
  testthat::expect_equal(out$stat, c(4, 8, 12))
  testthat::expect_equal(out$pvalue, c(0.04, 0.008, 0.003))
})

testthat::test_that("R mcols row-name preservation covers vector matrix and array payloads", {
  rowNames <- c("g1", "g2")
  helperEnvironment <- if ("rsdeseq2" %in% loadedNamespaces()) {
    asNamespace("rsdeseq2")
  } else {
    .GlobalEnv
  }
  applyMcolsRowNames <- get(
    ".rsdeseq2_apply_mcols_row_names",
    envir = helperEnvironment,
    inherits = TRUE
  )
  vector <- applyMcolsRowNames(c(10, 20), rowNames)
  matrix <- applyMcolsRowNames(matrix(1:4, nrow = 2), rowNames)
  array <- applyMcolsRowNames(array(1:8, dim = c(2, 2, 2)), rowNames)

  testthat::expect_equal(names(vector), rowNames)
  testthat::expect_equal(rownames(matrix), rowNames)
  testthat::expect_equal(dimnames(array)[[1L]], rowNames)

  named <- c(g2 = 20, g1 = 10)
  testthat::expect_identical(applyMcolsRowNames(named, rowNames), named)
})
