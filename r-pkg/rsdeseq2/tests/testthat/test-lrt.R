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
  filtered <- results(fit, name = "condition_B_vs_A")

  testthat::expect_s3_class(fit, "rsdeseq2PrimitiveLrtFit")
  testthat::expect_equal(out$log2FoldChange, c(2, 4))
  testthat::expect_equal(out$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(out$stat, c(6, 8))
  testthat::expect_equal(out$pvalue, stats::pchisq(c(6, 8), df = 1, lower.tail = FALSE))
  testthat::expect_equal(attr(out, "comparison"), "coefficient condition_B_vs_A; LRT full vs reduced")
  testthat::expect_equal(attr(filtered, "resultName"), "condition_B_vs_A")
  testthat::expect_equal(
    attr(filtered, "comparison"),
    "coefficient condition_B_vs_A; LRT full vs reduced"
  )
  testthat::expect_equal(attr(filtered, "contrast")[["Intercept"]], 0)
  testthat::expect_equal(attr(filtered, "contrast")[["condition_B_vs_A"]], 1)
  testthat::expect_true(attr(filtered, "independentFiltering")$enabled)
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

testthat::test_that("R wrapper LRT preserves transformed factor metadata from objects", {
  transformed <- "ordered(condition, levels = c(\"B\", \"A\"))"
  transformedAlias <- make.names(transformed)
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = matrix(
      c(
        1, 2,
        1, 2
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", paste0(transformed, "_A_vs_B")))
    ),
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = matrix(
      c(
        0, 0, 0, 0, 50, 60,
        10, 12, 20, 24, 50, 60
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), paste0("s", 1:6))
    )
  )
  cov <- array(0, dim = c(2, 2, 2))
  cov[, 2, 2] <- 0.25
  object$betaCovariance <- cov
  object$colData <- data.frame(check.names = FALSE, row.names = colnames(object$counts))
  object$colData[[transformed]] <- factor(c("A", "A", "B", "B", "C", "C"), levels = c("B", "A", "C"))

  fit <- nbinomLRTRust(object = object)
  out <- results(
    fit,
    contrast = c(transformedAlias, "A", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$stat[[1]], 6)
  testthat::expect_equal(out$pvalue[[1]], 0.02)
  testthat::expect_equal(attr(out, "resultName"), paste0(transformedAlias, "_A_vs_B"))
  testthat::expect_equal(
    attr(out, "comparison"),
    paste0("factor-level contrast: ", transformed, " A vs B; LRT full vs reduced")
  )
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(fit$factorLevels[[transformed]], c("B", "A", "C"))
})

testthat::test_that("R wrapper LRT preserves cleaned character contrast reference metadata", {
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = matrix(
      c(
        1, 2, 3,
        2, 4, 7
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
    ),
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = matrix(
      c(
        0, 0, 0, 0, 50, 60,
        10, 12, 20, 24, 50, 60
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), paste0("s", 1:6))
    )
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object$betaCovariance <- cov
  object$colData <- data.frame(check.names = FALSE, row.names = colnames(object$counts))
  object$colData[["cell type"]] <- factor(
    c("A cell", "A cell", "B cell", "B cell", "C cell", "C cell"),
    levels = c("A cell", "B cell", "C cell")
  )

  fit <- nbinomLRTRust(object = object)
  out <- results(
    fit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    reference = "A.cell",
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(out$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(out$stat, c(6, 8))
  testthat::expect_equal(out$pvalue, c(0.02, 0.01))
  testthat::expect_equal(attr(out, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(attr(out, "comparison"), "factor-level contrast: cell type C cell vs B cell; LRT full vs reduced")
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(fit$factorLevels[["cell type"]], c("A cell", "B cell", "C cell"))
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

testthat::test_that("R wrapper LRT extracts factor reference aliases from fitted objects", {
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = matrix(
      c(
        1, 2, 3,
        1, 2, 3
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
    ),
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    factorLevels = list("cell type" = c("B cell", "A cell", "C cell")),
    factorReferenceLevels = c("cell.type" = "A.cell")
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object$betaCovariance <- cov

  fit <- nbinomLRTRust(object = object)
  out <- results(
    fit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(fit$factorReferences[["cell type"]], "A cell")
  testthat::expect_equal(out$log2FoldChange, c(1, 1))
  testthat::expect_equal(out$stat, c(6, 8))
  testthat::expect_equal(out$pvalue, c(0.02, 0.01))
  testthat::expect_equal(attr(out, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(
    attr(out, "comparison"),
    "factor-level contrast: cell type C cell vs B cell; LRT full vs reduced"
  )
  testthat::expect_equal(attr(out, "contrast")[["cell type_B cell_vs_A cell"]], -1)
  testthat::expect_equal(attr(out, "contrast")[["cell type_C cell_vs_A cell"]], 1)
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

testthat::test_that("R wrapper LRT resolves ordered formula-transform character contrasts", {
  transformed <- "ordered(condition, levels = c(\"B\", \"A\"))"
  transformedAlias <- make.names(transformed)
  beta <- matrix(
    c(
      1, 2,
      1, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", paste0(transformed, "_A_vs_B")))
  )
  cov <- array(0, dim = c(2, 2, 2))
  cov[, 2, 2] <- 0.25
  counts <- matrix(
    c(
      0, 0, 0, 0, 50, 60,
      10, 12, 20, 24, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  colData <- data.frame(
    check.names = FALSE,
    row.names = colnames(counts)
  )
  colData[[transformed]] <- factor(c("A", "A", "B", "B", "C", "C"), levels = c("B", "A", "C"))
  fit <- lrtFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01),
    counts = counts,
    colData = colData,
    factorLevels = setNames(list(c("B", "A", "C")), transformed)
  )

  out <- results(
    fit,
    contrast = c(transformedAlias, "A", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$stat[[1]], 6)
  testthat::expect_equal(out$pvalue[[1]], 0.02)
  testthat::expect_equal(out$log2FoldChange[[2]], 2)
  testthat::expect_equal(attr(out, "resultName"), paste0(transformedAlias, "_A_vs_B"))
  testthat::expect_equal(
    attr(out, "comparison"),
    paste0("factor-level contrast: ", transformed, " A vs B; LRT full vs reduced")
  )
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(attr(out, "contrast")[[paste0(transformed, "_A_vs_B")]], 1)
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
