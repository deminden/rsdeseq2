testthat::test_that("R fitted mcols representations preserve complete Wald and LRT results", {
  mcols <- .rsdeseq2_mcols_fixture()
  lrt_mcols <- cbind(mcols, stat = c(4, 3), pvalue = c(0.01, 0.02))
  counts <- .rsdeseq2_count_fixture(c(0, 0, 10, 0, 10, 0, 0, 5))
  design <- matrix(c(1, 0, 1, 0, 1, 1, 1, 1), nrow = 4, byrow = TRUE,
                   dimnames = list(colnames(counts), c("Intercept", "condition_B_vs_A")))
  fields <- list(assays = list(counts = counts), sampleLevels = c("A", "A", "B", "B"),
                 modelMatrix = design)
  s4 <- function(mcols) .rsdeseq2_s4_fixture(
    "Rsdeseq2FittedMock",
    c(mcols = "data.frame", assays = "list", sampleLevels = "character", modelMatrix = "matrix"),
    c(list(mcols = mcols), fields)
  )
  # The matrix case deliberately lacks count dimnames/modelMatrix: gene labels
  # must come from mcols, and sampleLevels alone must resolve the contrast.
  matrix_object <- list(mcols = as.matrix(lrt_mcols), counts = unname(counts),
                        sampleLevels = c("A", "A", "B", "B"))
  cases <- list(
    list_fields = list(wald = list(mcols = mcols), lrt = c(list(mcols = lrt_mcols), fields)),
    matrix_mcols = list(
      wald = utils::modifyList(matrix_object, list(mcols = as.matrix(mcols))), lrt = matrix_object
    ),
    s4_slots = list(wald = s4(mcols), lrt = s4(lrt_mcols))
  )
  for (case in names(cases)) {
    objects <- cases[[case]]
    testthat::expect_equal(resultsNames(objects$wald), colnames(design), info = case)
    wald <- results(objects$wald, name = "condition_B_vs_A", independentFiltering = FALSE)
    lrt <- results(objects$lrt, contrast = c("condition", "B", "A"), independentFiltering = FALSE)
    .rsdeseq2_expect_result(wald, c(2, -2), c(0.5, 0.5), info = case)
    .rsdeseq2_expect_result(lrt, c(2, -2), c(0.5, 0.5), c(4, 3), c(0.01, 0.02), info = case)
    # A name selects the reported effect; supplied LRT statistics still govern
    # inference. In particular, matrix mcols with stat/pvalue must dispatch LRT.
    .rsdeseq2_expect_result(
      results(objects$lrt, name = "condition_B_vs_A", independentFiltering = FALSE),
      c(2, -2), c(0.5, 0.5), c(4, 3), c(0.01, 0.02), info = case
    )
    packaged <- nbinomWaldTestRust(object = objects$wald)
    testthat::expect_s3_class(packaged, "rsdeseq2PrimitiveWaldFit")
    .rsdeseq2_expect_result(
      results(packaged, name = "condition_B_vs_A", independentFiltering = FALSE),
      c(2, -2), c(0.5, 0.5), info = case
    )
  }
})

testthat::test_that("R model matrix access routes align samples and preserve numeric contrasts", {
  object <- .rsdeseq2_primitive_fixture(c(1, 2, 3, 1, 2, 3), variances = c(0, 0.25, 0.25))
  object$counts <- .rsdeseq2_count_fixture(c(10, 12, 0, 0, 0, 0, 10, 12, 30, 36, 50, 60))
  design <- cbind(Intercept = 1, condition_B_vs_A = rep(c(0, 1, 0), each = 2),
                  condition_C_vs_A = rep(c(0, 0, 1), each = 2))
  rownames(design) <- colnames(object$counts)
  shuffled <- design[c("s5", "s1", "s3", "s6", "s2", "s4"), , drop = FALSE]
  list_object <- c(object, list(designMatrix = shuffled))
  mcols <- cbind(.rsdeseq2_mcols_fixture(), condition_C_vs_A = c(3, 3), SE_condition_C_vs_A = 0.5)
  mcols$condition_B_vs_A <- c(2, 2)
  s4_object <- .rsdeseq2_s4_fixture(
    "Rsdeseq2FittedDesignAliasMock",
    c(mcols = "data.frame", betaCovariance = "array", counts = "matrix", fullModelMatrix = "matrix"),
    list(mcols = mcols, betaCovariance = object$betaCovariance, counts = object$counts,
         fullModelMatrix = shuffled)
  )
  attribute_object <- .rsdeseq2_primitive_fixture(
    c(1, 2, 5, 1, 4, 7), c("Intercept", "condition_B_vs_A", "batch_Y_vs_X"), c(0, 0.25, 0.25)
  )
  attribute_object$counts <- .rsdeseq2_count_fixture(c(5, 6, 0, 0, 5, 6, 2, 3))
  attr(attribute_object, "modelMatrix") <- matrix(
    c(1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0), nrow = 4, byrow = TRUE,
    dimnames = list(c("s3", "s4", "s1", "s2"), colnames(attribute_object$beta))
  )
  cases <- list(
    designMatrix_field = list(object = list_object, lfc = c(0, -1)),
    fullModelMatrix_slot = list(object = s4_object, lfc = c(0, -1)),
    modelMatrix_attribute = list(object = attribute_object, lfc = c(0, -3))
  )
  for (case in names(cases)) {
    input <- cases[[case]]
    fit <- nbinomWaldTestRust(object = input$object)
    testthat::expect_equal(rownames(fit$modelMatrix), colnames(fit$counts), info = case)
    # Exercise direct dispatch as well as packaging, with full independent answers.
    for (route in list(input$object, fit)) {
      out <- results(route, contrast = c(0, 1, -1), independentFiltering = FALSE)
      .rsdeseq2_expect_result(out, input$lfc, rep(sqrt(0.5), 2), info = case)
      testthat::expect_equal(attr(out, "contrastAllZero"), "numeric", info = case)
    }
  }
  badObject <- list_object
  rownames(badObject$designMatrix)[[1L]] <- "unknown"
  testthat::expect_error(
    results(badObject, contrast = c(0, 1, -1), independentFiltering = FALSE),
    "modelMatrix row names must contain all count sample names"
  )
})

testthat::test_that("R coefficient metadata aliases preserve selected coefficient numerics", {
  # B is not the default (last) coefficient, so ignoring name changes answers.
  object <- .rsdeseq2_primitive_fixture(c(1, 2, 3, 1, -2, 7))
  cases <- list(
    list_aliases = list(baseMean = object$baseMean, coefficients = object$beta,
                        coefNames = colnames(object$beta), coefCovariance = object$betaCovariance),
    s4_aliases = .rsdeseq2_s4_fixture(
      "Rsdeseq2FittedCoefAliasMock",
      c(baseMean = "numeric", coef = "matrix", resultNames = "character", betaCov = "array"),
      list(baseMean = object$baseMean, coef = object$beta, resultNames = colnames(object$beta),
           betaCov = object$betaCovariance)
    )
  )
  for (case in names(cases)) {
    testthat::expect_equal(resultsNames(cases[[case]]), colnames(object$beta), info = case)
    .rsdeseq2_expect_result(
      results(cases[[case]], name = "condition_B_vs_A", independentFiltering = FALSE),
      c(2, -2), c(0.5, 0.5), info = case
    )
  }
})

testthat::test_that("R factor level aliases preserve nonreference contrast numerics", {
  object <- .rsdeseq2_primitive_fixture(c(1, 2, 3, 1, 4, 7))
  cases <- list(
    factorLevelNames_field = c(object, list(factorLevelNames = list(condition = c("A", "B", "C")))),
    levelNames_slot = .rsdeseq2_s4_fixture(
      "Rsdeseq2FittedFactorLevelAliasMock",
      c(baseMean = "numeric", beta = "matrix", betaCovariance = "array", levelNames = "list"),
      c(object, list(levelNames = list(condition = c("A", "B", "C"))))
    )
  )
  for (case in names(cases)) {
    out <- results(cases[[case]], contrast = c("condition", "C", "B"), independentFiltering = FALSE)
    .rsdeseq2_expect_result(out, c(1, 3), rep(sqrt(1.25), 2), info = case)
    testthat::expect_equal(attr(out, "contrast"),
                           c(Intercept = 0, condition_B_vs_A = -1, condition_C_vs_A = 1), info = case)
  }
})
