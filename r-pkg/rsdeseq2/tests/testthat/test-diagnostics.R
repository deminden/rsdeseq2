testthat::test_that("R diagnostic metadata helper returns empty pre-GLM shape", {
  testthat::expect_identical(
    rsdeseq2DiagnosticSchemaRust(),
    c(
      "betaConv",
      "fullBetaConv",
      "reducedBetaConv",
      "betaIter",
      "reducedBetaIter",
      "deviance",
      "maxCooks"
    )
  )
  testthat::expect_identical(
    rsdeseq2DiagnosticSchemaRust(native = TRUE),
    rsdeseq2DiagnosticSchemaRust()
  )

  metadata <- deseq2McolsDiagnosticsRust(2, rowNames = c("g1", "g2"))

  testthat::expect_s3_class(metadata, "data.frame")
  testthat::expect_identical(dim(metadata), c(2L, 0L))
  testthat::expect_identical(rownames(metadata), c("g1", "g2"))
})

testthat::test_that("R diagnostic metadata helper exposes Wald mcols names", {
  metadata <- deseq2McolsDiagnosticsRust(
    2,
    test = "Wald",
    rowNames = c("g1", "g2"),
    betaConv = c(TRUE, FALSE),
    betaIter = c(3, 0),
    deviance = c(10.5, NA_real_),
    maxCooks = c(0.1, NA_real_)
  )

  testthat::expect_identical(
    names(metadata),
    c("betaConv", "betaIter", "deviance", "maxCooks")
  )
  testthat::expect_identical(metadata$betaConv, c(TRUE, FALSE))
  testthat::expect_identical(metadata$betaIter, c(3L, 0L))
  testthat::expect_equal(metadata$deviance, c(10.5, NA_real_))
  testthat::expect_identical(rownames(metadata), c("g1", "g2"))
})

testthat::test_that("R diagnostic metadata helper exposes LRT mcols names", {
  metadata <- deseq2McolsDiagnosticsRust(
    2,
    test = "LRT",
    fullBetaConv = c(FALSE, TRUE),
    reducedBetaConv = c(FALSE, TRUE),
    betaIter = c(0, 4),
    reducedBetaIter = c(0, 2),
    deviance = c(NA_real_, 20.25),
    maxCooks = c(NA_real_, 0.2)
  )

  testthat::expect_identical(
    names(metadata),
    c(
      "fullBetaConv",
      "reducedBetaConv",
      "betaIter",
      "reducedBetaIter",
      "deviance",
      "maxCooks"
    )
  )
  testthat::expect_identical(metadata$fullBetaConv, c(FALSE, TRUE))
  testthat::expect_identical(metadata$reducedBetaConv, c(FALSE, TRUE))
  testthat::expect_identical(metadata$betaIter, c(0L, 4L))
  testthat::expect_identical(metadata$reducedBetaIter, c(0L, 2L))
})

testthat::test_that("R diagnostic metadata helper validates primitive shapes", {
  testthat::expect_error(deseq2McolsDiagnosticsRust(-1))
  testthat::expect_error(deseq2McolsDiagnosticsRust(2, rowNames = "g1"))
  testthat::expect_error(deseq2McolsDiagnosticsRust(2, test = "Wald", betaConv = TRUE))
  testthat::expect_error(deseq2McolsDiagnosticsRust(2, betaIter = c(1, -1)))
  testthat::expect_error(deseq2McolsDiagnosticsRust(2, deviance = c("x", "y")))
})
