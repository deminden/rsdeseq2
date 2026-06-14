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
  filtered <- results(fit, name = "condition_B_vs_A")

  testthat::expect_s3_class(fit, "rsdeseq2PrimitiveWaldFit")
  testthat::expect_equal(out$log2FoldChange, c(2, 4))
  testthat::expect_equal(out$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(attr(out, "comparison"), "coefficient condition_B_vs_A")
  testthat::expect_equal(attr(filtered, "resultName"), "condition_B_vs_A")
  testthat::expect_equal(attr(filtered, "comparison"), "coefficient condition_B_vs_A")
  testthat::expect_equal(attr(filtered, "contrast")[["Intercept"]], 0)
  testthat::expect_equal(attr(filtered, "contrast")[["condition_B_vs_A"]], 1)
  testthat::expect_true(attr(filtered, "independentFiltering")$enabled)
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

testthat::test_that("R wrapper Wald test preserves transformed factor metadata from objects", {
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

  fit <- nbinomWaldTestRust(object = object)
  out <- results(
    fit,
    contrast = c(transformedAlias, "A", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_equal(attr(out, "resultName"), paste0(transformedAlias, "_A_vs_B"))
  testthat::expect_equal(
    attr(out, "comparison"),
    paste0("factor-level contrast: ", transformed, " A vs B")
  )
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(fit$factorLevels[[transformed]], c("B", "A", "C"))
})

testthat::test_that("R wrapper Wald test preserves cleaned character contrast reference metadata", {
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

  fit <- nbinomWaldTestRust(object = object)
  out <- results(
    fit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    reference = "A.cell",
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(out$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(attr(out, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(attr(out, "comparison"), "factor-level contrast: cell type C cell vs B cell")
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(fit$factorLevels[["cell type"]], c("A cell", "B cell", "C cell"))
})

testthat::test_that("R wrapper Wald test extracts factor reference aliases from fitted objects", {
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
    factorLevels = list("cell type" = c("B cell", "A cell", "C cell")),
    referenceLevels = c("cell.type" = "A.cell")
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object$betaCovariance <- cov

  fit <- nbinomWaldTestRust(object = object)
  out <- results(
    fit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(fit$factorReferences[["cell type"]], "A cell")
  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(attr(out, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(attr(out, "comparison"), "factor-level contrast: cell type C cell vs B cell")
  testthat::expect_equal(attr(out, "contrast")[["cell type_B cell_vs_A cell"]], -1)
  testthat::expect_equal(attr(out, "contrast")[["cell type_C cell_vs_A cell"]], 1)
})

testthat::test_that("R wrapper Wald test validates primitive inputs", {
  testthat::expect_error(nbinomWaldTestRust())
  testthat::expect_error(nbinomWaldTestRust(baseMean = 1, beta = matrix(1, nrow = 2)))
})
