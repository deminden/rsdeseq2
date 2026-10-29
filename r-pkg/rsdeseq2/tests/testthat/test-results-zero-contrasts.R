# Contrast all-zero rules across factor levels and observation weights.

testthat::test_that("R results does not apply contrast all-zero handling for coefficient names", {
  beta <- matrix(
    c(
      0, 2,
      0, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A"))
  )
  betaSE <- matrix(c(0.1, 0.5, 0.1, 0.5), nrow = 2, byrow = TRUE)
  counts <- matrix(
    c(
      0, 0,
      10, 12
    ),
    nrow = 2,
    byrow = TRUE
  )
  waldFit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = betaSE,
    counts = counts
  )
  lrtFit <- lrtFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = betaSE,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts
  )

  wald <- results(waldFit, name = "condition_B_vs_A", independentFiltering = FALSE)
  lrt <- results(lrtFit, name = "condition_B_vs_A", independentFiltering = FALSE)

  testthat::expect_equal(wald$log2FoldChange, c(2, 4))
  testthat::expect_equal(
    wald$pvalue[[1]],
    stats::pnorm(abs(2 / 0.5), lower.tail = FALSE) * 2,
    tolerance = 1e-15
  )
  testthat::expect_equal(lrt$log2FoldChange, c(2, 4))
  testthat::expect_equal(lrt$pvalue, c(0.02, 0.01))
  testthat::expect_equal(attr(wald, "contrastAllZero"), "none")
  testthat::expect_equal(attr(lrt, "contrastAllZero"), "none")
})

testthat::test_that("R results applies character and numeric contrast all-zero semantics", {
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
  fit <- waldFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
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
  list <- results(
    fit,
    contrast = list("condition_C_vs_A", "condition_B_vs_A"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(character$log2FoldChange[[1]], 0)
  testthat::expect_equal(character$stat[[1]], 0)
  testthat::expect_equal(character$pvalue[[1]], 1)
  testthat::expect_true(numeric$pvalue[[1]] < 1)
  testthat::expect_equal(attr(character, "contrastAllZero"), "character")
  testthat::expect_equal(attr(numeric, "contrastAllZero"), "numeric")
  testthat::expect_equal(list$log2FoldChange[[1]], 1)
  testthat::expect_true(list$pvalue[[1]] < 1)

  expandedBeta <- matrix(
    c(
      1, 2, 3,
      1, 2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("conditionA", "conditionB", "conditionC"))
  )
  expandedCov <- array(0, dim = c(2, 3, 3))
  expandedCov[, 1, 1] <- 0.25
  expandedCov[, 2, 2] <- 0.25
  expandedCov[, 3, 3] <- 0.25
  expandedModelMatrix <- diag(3)[c(1, 1, 2, 2, 3, 3), ]
  expandedFit <- waldFitRust(
    baseMean = c(10, 20),
    beta = expandedBeta,
    betaCovariance = expandedCov,
    counts = counts,
    sampleLevels = c("A", "A", "B", "B", "C", "C"),
    modelMatrix = expandedModelMatrix
  )
  expandedList <- results(
    expandedFit,
    contrast = list("conditionB", "conditionA"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(expandedList$log2FoldChange[[1]], 0)
  testthat::expect_equal(expandedList$pvalue[[1]], 1)

  cleanedLevelFit <- waldFitRust(
    baseMean = c(10, 20),
    beta = matrix(
      c(
        1, 2,
        1, 2
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", "cell.type_B.cell_vs_A.cell"))
    ),
    betaCovariance = array(c(0, 0, 0, 0.25, 0, 0, 0, 0.25), dim = c(2, 2, 2)),
    counts = counts,
    sampleLevels = c("A cell", "A cell", "B cell", "B cell", "C cell", "C cell"),
    modelMatrix = modelMatrix[, 1:2, drop = FALSE]
  )
  cleanedLevelOut <- results(
    cleanedLevelFit,
    contrast = c("cell.type", "B.cell", "A.cell"),
    reference = "A.cell",
    independentFiltering = FALSE
  )
  testthat::expect_equal(cleanedLevelOut$log2FoldChange[[1]], 0)
  testthat::expect_equal(cleanedLevelOut$pvalue[[1]], 1)
  testthat::expect_equal(attr(cleanedLevelOut, "contrastAllZero"), "character")
  testthat::expect_equal(attr(cleanedLevelOut, "resultName"), "cell.type_B.cell_vs_A.cell")
  testthat::expect_equal(attr(cleanedLevelOut, "comparison"), "factor-level contrast: cell.type B cell vs A cell")
})

testthat::test_that("R results uses colData factor values for character contrast all-zero handling", {
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
    colData = data.frame(
      condition = factor(c("C", "A", "B", "C", "A", "B"), levels = c("A", "B", "C")),
      row.names = c("s5", "s1", "s3", "s6", "s2", "s4")
    )
  )

  out <- results(
    object,
    contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(rownames(object$colData), c("s5", "s1", "s3", "s6", "s2", "s4"))

  aliasObject <- object
  colnames(aliasObject$beta) <- c("Intercept", "cell.type_B_vs_A", "cell.type_C_vs_A")
  names(aliasObject$colData) <- "cell type"
  cleaned <- results(
    aliasObject,
    contrast = c("cell.type", "B", "A"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(cleaned$log2FoldChange[[1]], 0)
  testthat::expect_equal(cleaned$pvalue[[1]], 1)
  testthat::expect_equal(attr(cleaned, "contrastAllZero"), "character")

  aliasObject$factorLevels <- list("cell type" = c("A", "B", "C"))
  reference <- results(
    aliasObject,
    contrast = c("cell.type", "A", "B"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(attr(reference, "contrast")[["cell.type_B_vs_A"]], -1)

  levelAliasObject <- object
  colnames(levelAliasObject$beta) <- c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell")
  names(levelAliasObject$colData) <- "cell type"
  levelAliasObject$colData[["cell type"]] <- factor(
    paste(as.character(levelAliasObject$colData[["cell type"]]), "cell"),
    levels = c("A cell", "B cell", "C cell")
  )
  levelAliasObject$factorLevels <- list("cell type" = c("A cell", "B cell", "C cell"))
  cleanedLevels <- results(
    levelAliasObject,
    contrast = c("cell.type", "B.cell", "A.cell"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(cleanedLevels$log2FoldChange[[1]], 0)
  testthat::expect_equal(cleanedLevels$pvalue[[1]], 1)
  testthat::expect_equal(attr(cleanedLevels, "contrastAllZero"), "character")
  testthat::expect_equal(attr(cleanedLevels, "resultName"), "cell.type_B.cell_vs_A.cell")
  testthat::expect_equal(attr(cleanedLevels, "comparison"), "factor-level contrast: cell type B cell vs A cell")
  testthat::expect_equal(attr(cleanedLevels, "contrast")[["cell type_B cell_vs_A cell"]], 1)
  cleanedReference <- results(
    levelAliasObject,
    contrast = c("cell.type", "C.cell", "B.cell"),
    reference = "A.cell",
    independentFiltering = FALSE
  )
  testthat::expect_equal(attr(cleanedReference, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(attr(cleanedReference, "comparison"), "factor-level contrast: cell type C cell vs B cell")
  testthat::expect_equal(attr(cleanedReference, "contrast")[["cell type_C cell_vs_A cell"]], 1)
  testthat::expect_equal(attr(cleanedReference, "contrast")[["cell type_B cell_vs_A cell"]], -1)

  ambiguousLevelObject <- object
  colnames(ambiguousLevelObject$beta) <- c("Intercept", "cell.type_B.cell_vs_A.cell", "cell.type_C.cell_vs_A.cell")
  names(ambiguousLevelObject$colData) <- "cell type"
  ambiguousLevelObject$colData[["cell type"]] <- c("A cell", "A cell", "B cell", "B-cell", "C cell", "C cell")
  testthat::expect_error(
    results(
      ambiguousLevelObject,
      contrast = c("cell.type", "B.cell", "A.cell"),
      independentFiltering = FALSE
    ),
    "factorLevels entries must resolve to unique R-cleaned level aliases"
  )
  ambiguousReferenceObject <- levelAliasObject
  ambiguousReferenceObject$factorLevels <- list("cell type" = c("A cell", "A-cell", "B cell", "C cell"))
  testthat::expect_error(
    results(
      ambiguousReferenceObject,
      contrast = c("cell.type", "C.cell", "B.cell"),
      reference = "A.cell",
      independentFiltering = FALSE
    ),
    "factorLevels entries must resolve to unique R-cleaned level aliases"
  )

  aliasObject$colData <- data.frame(
    "cell type" = factor(c("A", "A", "B", "B", "C", "C")),
    "cell-type" = factor(c("A", "A", "B", "B", "C", "C")),
    check.names = FALSE
  )
  aliasObject$factorLevels <- NULL
  testthat::expect_error(
    results(aliasObject, contrast = c("cell.type", "B", "A"), independentFiltering = FALSE),
    "factorLevels names must resolve to unique R-cleaned aliases"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = aliasObject$beta,
      betaCovariance = aliasObject$betaCovariance,
      counts = aliasObject$counts,
      colData = aliasObject$colData
    ),
    "factorLevels names must resolve to unique R-cleaned aliases"
  )
})
