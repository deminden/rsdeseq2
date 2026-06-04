testthat::test_that("R wrapper results are explicit future work", {
  testthat::expect_error(results(NULL))
})

testthat::test_that("R wrapper resolves DESeq2 results contrast character triplets", {
  names <- c("Intercept", "group_2_vs_1", "condition_2_vs_1", "condition_3_vs_1")

  direct <- resolveResultsContrastRust(c("condition", "2", "1"), names)
  reverse <- resolveResultsContrastRust(c("condition", "1", "3"), names)
  shared <- resolveResultsContrastRust(c("condition", "2", "3"), names)

  testthat::expect_s3_class(direct, "rsdeseq2ResultsContrast")
  testthat::expect_equal(unname(direct$numeric), c(0, 0, 1, 0))
  testthat::expect_equal(unname(reverse$numeric), c(0, 0, 0, -1))
  testthat::expect_equal(unname(shared$numeric), c(0, 0, 1, -1))
  testthat::expect_equal(shared$resultName, "condition_2_vs_3")
  testthat::expect_equal(shared$comparison, "factor-level contrast: condition 2 vs 3")
  testthat::expect_equal(shared$allZero$type, "character")
  testthat::expect_equal(shared$allZero$numerator, "2")
  testthat::expect_equal(shared$allZero$denominator, "3")
})

testthat::test_that("R wrapper resolves DESeq2 results contrast list and numeric forms", {
  names <- c("Intercept", "condition_B_vs_A", "batch_Y_vs_X")

  listContrast <- resolveResultsContrastRust(
    list("condition_B_vs_A", "batch_Y_vs_X"),
    names,
    listValues = c(0.5, -0.5)
  )
  positiveOnly <- resolveResultsContrastRust(list("condition_B_vs_A"), names)
  numericContrast <- resolveResultsContrastRust(c(0, 1, -1), names)

  testthat::expect_equal(unname(listContrast$numeric), c(0, 0.5, -0.5))
  testthat::expect_equal(listContrast$resultName, "contrast")
  testthat::expect_equal(
    listContrast$comparison,
    "coefficient list contrast: 0.5 condition_B_vs_A vs 0.5 batch_Y_vs_X"
  )
  testthat::expect_equal(listContrast$allZero$type, "numeric")
  testthat::expect_equal(unname(positiveOnly$numeric), c(0, 1, 0))
  testthat::expect_equal(
    positiveOnly$comparison,
    "coefficient list contrast: condition_B_vs_A effect"
  )
  testthat::expect_equal(unname(numericContrast$numeric), c(0, 1, -1))
  testthat::expect_equal(names(numericContrast$numeric), names)
  testthat::expect_equal(numericContrast$comparison, "primitive numeric contrast")
})

testthat::test_that("R wrapper exposes primitive results names for contrast resolution", {
  beta <- matrix(
    c(
      1, 2, 3,
      2, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  wald <- waldFitRust(c(10, 20), beta = beta, betaSE = matrix(1, nrow = 2, ncol = 3))
  lrt <- lrtFitRust(c(10, 20), beta = beta, betaSE = matrix(1, nrow = 2, ncol = 3), lrtStat = c(6, 8), lrtDf = 1)

  testthat::expect_equal(resultsNamesRust(wald), colnames(beta))
  testthat::expect_equal(resultsNamesRust(lrt), colnames(beta))
  testthat::expect_equal(resultsNamesRust(list(beta = beta)), colnames(beta))
  testthat::expect_equal(resultsNamesRust(beta), colnames(beta))
  testthat::expect_equal(resultsNames(wald), colnames(beta))
  testthat::expect_equal(resultsNames(list(beta = beta)), colnames(beta))
  resolved <- resolveResultsContrastRust(
    c("condition", "C", "B"),
    resultsNames(wald),
    reference = "A"
  )
  testthat::expect_equal(unname(resolved$numeric), c(0, -1, 1))
  testthat::expect_error(resultsNamesRust(list(beta = matrix(1, nrow = 1))))
  testthat::expect_error(resultsNamesRust(NULL))
})

testthat::test_that("R wrapper resolves R-cleaned results contrast aliases", {
  names <- c("(Intercept)", "if.", "condition.B.1")

  coefList <- resolveResultsContrastRust(list("if", "condition-B 1"), names)
  expanded <- resolveResultsContrastRust(c("condition", "B-1", "A-1"), c("conditionA.1", "conditionB.1"))
  sharedNoSeparator <- resolveResultsContrastRust(c("condition", "C", "B"), c("Intercept", "conditionB_vs_A", "conditionC_vs_A"))

  testthat::expect_equal(unname(coefList$numeric), c(0, 1, -1))
  testthat::expect_equal(unname(expanded$numeric), c(-1, 1))
  testthat::expect_equal(unname(sharedNoSeparator$numeric), c(0, -1, 1))
})

testthat::test_that("R wrapper rejects invalid results contrast shapes", {
  names <- c("Intercept", "condition_B_vs_A", "condition_C_vs_A")
  duplicated <- resolveResultsContrastRust(list(c("condition_B_vs_A", "condition_B_vs_A")), names)

  testthat::expect_error(resolveResultsContrastRust(c("condition", "B"), names))
  testthat::expect_error(resolveResultsContrastRust(c("condition", "B", "B"), names))
  testthat::expect_error(resolveResultsContrastRust(c("condition", "D", "A"), names))
  testthat::expect_error(resolveResultsContrastRust(c(0, 0, 0), names))
  testthat::expect_error(resolveResultsContrastRust(c(0, 1), names))
  testthat::expect_error(resolveResultsContrastRust(list(character()), names))
  testthat::expect_error(resolveResultsContrastRust(list("condition_B_vs_A", "condition_B_vs_A"), names))
  testthat::expect_equal(unname(duplicated$numeric), c(0, 1, 0))
  testthat::expect_error(resolveResultsContrastRust(list("condition_B_vs_A", "condition_C_vs_A"), names, listValues = c(0, -1)))
  testthat::expect_error(resolveResultsContrastRust(list("condition_B_vs_A", "condition_C_vs_A"), names, listValues = c(1, 1)))
  testthat::expect_error(
    resolveResultsContrastRust(
      c("condition", "C", "B"),
      c("Intercept", "condition_B_vs_A", "conditionB_vs_A", "condition_C_vs_A")
    ),
    "resolves ambiguously"
  )
})

testthat::test_that("R results computes primitive Wald coefficient and contrast tables", {
  beta <- matrix(
    c(
      1, 2, 5,
      2, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "batch_Y_vs_X"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  cov[, 2, 3] <- 0.1
  cov[, 3, 2] <- 0.1
  fit <- waldFitRust(
    baseMean = c(10, 20),
    beta = beta,
    betaCovariance = cov,
    betaSE = matrix(
      c(
        0.1, 0.5, 1,
        0.2, 0.5, 1
      ),
      nrow = 2,
      byrow = TRUE
    )
  )

  named <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)
  contrast <- results(
    fit,
    contrast = list("batch_Y_vs_X", "condition_B_vs_A"),
    independentFiltering = FALSE
  )
  numeric <- results(
    fit,
    contrast = c(0, -1, 1),
    independentFiltering = FALSE
  )

  testthat::expect_equal(named$log2FoldChange, c(2, 4))
  testthat::expect_equal(named$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(attr(named, "resultName"), "condition_B_vs_A")
  testthat::expect_equal(attr(named, "comparison"), "coefficient condition_B_vs_A")
  testthat::expect_equal(attr(named, "contrastAllZero"), "none")
  testthat::expect_equal(contrast$log2FoldChange, c(3, 3))
  testthat::expect_equal(contrast$lfcSE, rep(sqrt(1 + 0.25 - 2 * 0.1), 2))
  testthat::expect_equal(numeric$log2FoldChange, contrast$log2FoldChange)
  testthat::expect_equal(numeric$lfcSE, contrast$lfcSE)
})

testthat::test_that("R results accepts list-like primitive Wald and LRT inputs", {
  beta <- matrix(
    c(
      1, 2, 3,
      2, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  primitive <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov
  )

  wald <- results(
    primitive,
    contrast = c("condition", "C", "B"),
    reference = "A",
    independentFiltering = FALSE
  )
  lrt <- results(
    c(primitive, list(lrtStat = c(6, 8), lrtPvalue = c(0.02, 0.01))),
    contrast = c("condition", "C", "B"),
    reference = "A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(wald$log2FoldChange, c(1, 3))
  testthat::expect_equal(wald$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(lrt$log2FoldChange, c(1, 3))
  testthat::expect_equal(lrt$stat, c(6, 8))
  testthat::expect_equal(lrt$pvalue, c(0.02, 0.01))
  testthat::expect_equal(attr(wald, "resultName"), "condition_C_vs_B")
  testthat::expect_equal(attr(lrt, "resultName"), "condition_C_vs_B")
})

testthat::test_that("R results accepts DESeq2-shaped fitted list fields", {
  mcols <- data.frame(
    baseMean = c(10, 20),
    Intercept = c(1, 1),
    condition_B_vs_A = c(2, -2),
    SE_Intercept = c(0.1, 0.1),
    SE_condition_B_vs_A = c(0.5, 0.5),
    stat = c(4, 3),
    pvalue = c(0.01, 0.02),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )
  counts <- matrix(
    c(0, 0, 10, 0, 10, 0, 0, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  modelMatrix <- matrix(
    c(1, 0, 1, 0, 1, 1, 1, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("Intercept", "condition_B_vs_A"))
  )
  object <- list(
    mcols = mcols,
    assays = list(counts = counts),
    sampleLevels = c("A", "A", "B", "B"),
    modelMatrix = modelMatrix
  )

  wald <- results(
    list(mcols = mcols[, setdiff(colnames(mcols), c("stat", "pvalue")), drop = FALSE]),
    name = "condition_B_vs_A",
    independentFiltering = FALSE
  )
  lrt <- results(
    object,
    contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(wald$log2FoldChange, c(2, -2))
  testthat::expect_equal(wald$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(lrt$log2FoldChange[[2]], -2)
  testthat::expect_equal(lrt$stat, c(4, 3))
  testthat::expect_equal(lrt$pvalue, c(0.01, 0.02))
  testthat::expect_equal(rownames(lrt), c("g1", "g2"))
})

testthat::test_that("R results uses fitted factor levels as character contrast reference", {
  beta <- matrix(
    c(
      1, 2, 20, 3,
      1, 4, 40, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "conditionB_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(2, 4, 4))
  cov[, 2, 2] <- 0.25
  cov[, 4, 4] <- 1
  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    factorLevels = list(condition = c("A", "B", "C"))
  )

  out <- results(
    fit,
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )
  fromColData <- results(
    list(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      colData = data.frame(condition = factor(c("A", "B", "C"), levels = c("A", "B", "C")))
    ),
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange, c(1, 3))
  testthat::expect_equal(fromColData$log2FoldChange, out$log2FoldChange)
  testthat::expect_equal(attr(out, "contrast"), c(Intercept = 0, condition_B_vs_A = -1, conditionB_vs_A = 0, condition_C_vs_A = 1))
})

testthat::test_that("R results accepts S4-shaped fitted object slots", {
  if (!methods::isClass("Rsdeseq2FittedMock")) {
    methods::setClass(
      "Rsdeseq2FittedMock",
      slots = c(
        mcols = "data.frame",
        assays = "list",
        sampleLevels = "character",
        modelMatrix = "matrix"
      )
    )
  }
  mcols <- data.frame(
    baseMean = c(10, 20),
    Intercept = c(1, 1),
    condition_B_vs_A = c(2, -2),
    SE_Intercept = c(0.1, 0.1),
    SE_condition_B_vs_A = c(0.5, 0.5),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )
  counts <- matrix(
    c(0, 0, 10, 0, 10, 0, 0, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  modelMatrix <- matrix(
    c(1, 0, 1, 0, 1, 1, 1, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("Intercept", "condition_B_vs_A"))
  )
  waldObject <- methods::new(
    "Rsdeseq2FittedMock",
    mcols = mcols,
    assays = list(counts = counts),
    sampleLevels = c("A", "A", "B", "B"),
    modelMatrix = modelMatrix
  )
  lrtObject <- methods::new(
    "Rsdeseq2FittedMock",
    mcols = cbind(mcols, stat = c(4, 3), pvalue = c(0.01, 0.02)),
    assays = list(counts = counts),
    sampleLevels = c("A", "A", "B", "B"),
    modelMatrix = modelMatrix
  )

  wald <- results(waldObject, name = "condition_B_vs_A", independentFiltering = FALSE)
  lrt <- results(lrtObject, contrast = c("condition", "B", "A"), independentFiltering = FALSE)
  packaged <- nbinomWaldTestRust(object = waldObject)

  testthat::expect_equal(resultsNames(waldObject), c("Intercept", "condition_B_vs_A"))
  testthat::expect_equal(wald$log2FoldChange, c(2, -2))
  testthat::expect_equal(lrt$stat, c(4, 3))
  testthat::expect_equal(lrt$pvalue, c(0.01, 0.02))
  testthat::expect_s3_class(packaged, "rsdeseq2PrimitiveWaldFit")
})

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
})

testthat::test_that("R primitive Wald fit results validate required contrast diagnostics", {
  beta <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(NULL, c("Intercept", "condition_B_vs_A")))
  fit <- waldFitRust(c(1, 2), beta = beta, betaSE = matrix(c(0.1, 0.2, 0.1, 0.2), nrow = 2))

  contrastWins <- results(
    fit,
    name = "Intercept",
    contrast = c(0, 1),
    independentFiltering = FALSE
  )
  testthat::expect_equal(contrastWins$log2FoldChange, c(3, 4))
  testthat::expect_equal(attr(contrastWins, "resultName"), "contrast")
  testthat::expect_error(results(fit, contrast = c(1, -1), independentFiltering = FALSE))
  testthat::expect_error(waldFitRust(c(1), beta = beta))
})

testthat::test_that("R results validates character contrast sample levels", {
  beta <- matrix(
    c(
      0, 2,
      0, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A"))
  )
  counts <- matrix(
    c(
      0, 0,
      10, 12
    ),
    nrow = 2,
    byrow = TRUE
  )
  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaSE = matrix(c(0.1, 0.5, 0.1, 0.5), nrow = 2, byrow = TRUE),
    counts = counts,
    sampleLevels = c("A", "A")
  )

  testthat::expect_error(
    results(fit, contrast = c("condition", "B", "A"), independentFiltering = FALSE),
    "sampleLevels to contain numerator and denominator"
  )
})

testthat::test_that("R results supports thresholded Wald alternatives for primitive fits", {
  beta <- matrix(
    c(
      0, 2,
      0, -2,
      0, 0.2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("Intercept", "condition_B_vs_A"))
  )
  betaSE <- matrix(
    c(
      0.1, 0.5,
      0.1, 0.5,
      0.1, 0.5
    ),
    nrow = 3,
    byrow = TRUE
  )
  fit <- waldFitRust(c(10, 20, 30), beta = beta, betaSE = betaSE)

  greaterAbs <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greaterAbs",
    independentFiltering = FALSE
  )
  greaterAbs2014 <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greaterAbs2014",
    independentFiltering = FALSE
  )
  lessAbs <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "lessAbs",
    independentFiltering = FALSE
  )
  greater <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "greater",
    independentFiltering = FALSE
  )
  less <- results(
    fit,
    name = "condition_B_vs_A",
    lfcThreshold = 1,
    altHypothesis = "less",
    independentFiltering = FALSE
  )

  testthat::expect_equal(greaterAbs$stat[[1]], 4)
  testthat::expect_equal(greaterAbs$pvalue[[1]], 0.02275013293476686, tolerance = 1e-11)
  testthat::expect_equal(greaterAbs2014$stat[[1]], 2)
  testthat::expect_equal(greaterAbs2014$pvalue[[1]], 0.04550026389635842, tolerance = 1e-11)
  testthat::expect_equal(lessAbs$stat[[3]], 1.6)
  testthat::expect_equal(lessAbs$pvalue[[3]], 0.054799291699557995, tolerance = 1e-11)
  testthat::expect_equal(greater$stat[[1]], 2)
  testthat::expect_equal(greater$pvalue[[1]], 0.02275013194817921, tolerance = 1e-11)
  testthat::expect_equal(less$stat[[2]], -2)
  testthat::expect_equal(less$pvalue[[2]], 0.02275013194817921, tolerance = 1e-11)
  testthat::expect_equal(attr(greater, "lfcThreshold"), 1)
  testthat::expect_equal(attr(greater, "altHypothesis"), "greater")
})

testthat::test_that("R results supports thresholded contrast tests with t tails", {
  beta <- matrix(
    c(0, 2, 4),
    nrow = 1,
    dimnames = list("g1", c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(1, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  fit <- waldFitRust(c(g1 = 10), beta = beta, betaCovariance = cov)

  out <- results(
    fit,
    contrast = list("condition_C_vs_A", "condition_B_vs_A"),
    lfcThreshold = 1,
    altHypothesis = "greater",
    useT = TRUE,
    degreesOfFreedom = 10,
    independentFiltering = FALSE
  )

  testthat::expect_equal(out$log2FoldChange[[1]], 2)
  testthat::expect_equal(out$lfcSE[[1]], sqrt(0.5))
  testthat::expect_equal(out$stat[[1]], (2 - 1) / sqrt(0.5))
  testthat::expect_equal(
    out$pvalue[[1]],
    stats::pt((2 - 1) / sqrt(0.5), df = 10, lower.tail = FALSE),
    tolerance = 1e-15
  )
  testthat::expect_equal(attr(out, "degreesOfFreedom"), 10)
})

testthat::test_that("R results validates thresholded Wald options", {
  beta <- matrix(c(1, 2), nrow = 1, dimnames = list("g1", c("Intercept", "condition_B_vs_A")))
  betaSE <- matrix(c(0.1, 0.5), nrow = 1)
  fit <- waldFitRust(c(g1 = 10), beta = beta, betaSE = betaSE)

  testthat::expect_error(results(fit, lfcThreshold = -1))
  testthat::expect_error(results(fit, lfcThreshold = 0, altHypothesis = "lessAbs"))
  testthat::expect_error(results(fit, useT = "yes", degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = c(TRUE, FALSE), degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = NA, degreesOfFreedom = 10))
  testthat::expect_error(results(fit, useT = TRUE))
  testthat::expect_error(results(fit, useT = TRUE, degreesOfFreedom = 10, altHypothesis = "greaterAbsUPSHOT"))
})

testthat::test_that("R wrapper assembles primitive DESeq2-shaped result tables", {
  out <- resultsTableRust(
    baseMean = c(g1 = 10, g2 = 20, g3 = 30),
    log2FoldChange = c(1, -0.5, NA),
    lfcSE = c(0.2, 0.3, NA),
    stat = c(5, -1.5, NA),
    pvalue = c(0.001, 0.05, NA),
    dispersion = c(0.1, 0.2, 0.3),
    converged = c(TRUE, TRUE, FALSE)
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(
    colnames(out),
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged")
  )
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.001, 0.05, NA), method = "BH"))
})

testthat::test_that("R wrapper result table helper validates result vectors", {
  testthat::expect_error(resultsTableRust(c(1, NA)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1, 2)))
  testthat::expect_error(resultsTableRust(c(1, 2), padj = c(0.1, 0.2)))
  testthat::expect_error(resultsTableRust(c(1, 2), converged = c(TRUE)))

  out <- resultsTableRust(c(1, 2), rowNames = c("a", "b"))
  testthat::expect_equal(rownames(out), c("a", "b"))
  testthat::expect_true(all(is.na(out$pvalue)))
  testthat::expect_true(all(is.na(out$padj)))
})

testthat::test_that("R wrapper applies primitive independent filtering", {
  table <- resultsTableRust(
    baseMean = c(g1 = 0, g2 = 1, g3 = 2),
    stat = c(1, 1, 1),
    pvalue = c(0.01, 0.02, 0.03)
  )

  filtered <- applyIndependentFilteringRust(
    table,
    alpha = 0.5,
    theta = c(1, 1)
  )
  metadata <- attr(filtered, "independentFiltering")

  testthat::expect_equal(filtered$filtered, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(filtered$padj, c(NA, NA, 0.03))
  testthat::expect_true(metadata$enabled)
  testthat::expect_equal(metadata$selected, 1L)
  testthat::expect_equal(metadata$filterTheta, 1)
  testthat::expect_equal(metadata$filterThreshold, 2)
  testthat::expect_equal(metadata$numRej, c(1L, 1L))
})

testthat::test_that("R wrapper can disable primitive independent filtering", {
  table <- resultsTableRust(
    baseMean = c(0, 1, 2),
    pvalue = c(0.01, NA, 0.03)
  )

  filtered <- applyIndependentFilteringRust(table, enabled = FALSE)
  metadata <- attr(filtered, "independentFiltering")

  testthat::expect_equal(filtered$padj, stats::p.adjust(c(0.01, NA, 0.03), method = "BH"))
  testthat::expect_true(all(is.na(filtered$filtered)))
  testthat::expect_false(metadata$enabled)
  testthat::expect_equal(metadata$alpha, 0.1)
})

testthat::test_that("R wrapper independent filtering validates primitive inputs", {
  table <- resultsTableRust(c(1, 2), pvalue = c(0.1, 0.2))

  testthat::expect_error(applyIndependentFilteringRust(table[, c("baseMean", "pvalue")]))
  testthat::expect_error(applyIndependentFilteringRust(table, alpha = 1))
  testthat::expect_error(applyIndependentFilteringRust(table, theta = c(0.5)))
  testthat::expect_error(applyIndependentFilteringRust(table, theta = c(0, 2)))
})

testthat::test_that("R wrapper applies Cook's cutoff masking and BH adjustment", {
  out <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02, g3 = 0.5),
    maxCooks = c(0, 10, NA),
    cooksCutoff = 5
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02, g3 = 0.5),
    maxCooks = c(0, 10, NA),
    cooksCutoff = 5,
    native = TRUE
  )

  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(out$cooksOutlier, c(FALSE, TRUE, NA))
  testthat::expect_equal(out$pvalue, c(0.01, NA, 0.5))
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.01, NA, 0.5), method = "BH"))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper treats NaN Cook's values as missing", {
  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02, 0.03),
    maxCooks = c(NaN, NA_real_, 10),
    cooksCutoff = 5
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02, 0.03),
    maxCooks = c(NaN, NA_real_, 10),
    cooksCutoff = 5,
    native = TRUE
  )

  testthat::expect_equal(out$cooksOutlier, c(NA, NA, TRUE))
  testthat::expect_equal(out$pvalue, c(0.01, 0.02, NA))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper can disable Cook's cutoff masking", {
  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(100, 200),
    cooksCutoff = FALSE
  )

  testthat::expect_equal(out$pvalue, c(0.01, 0.02))
  testthat::expect_true(all(is.na(out$cooksOutlier)))
})

testthat::test_that("R wrapper low-count Cook's heuristic spares eligible rows", {
  counts <- matrix(c(1, 5, 6, 7, 9, 5, 6, 7), nrow = 2, byrow = TRUE)
  cooks <- matrix(c(10, 0.1, 0.2, 0.3, 10, 0.1, 0.2, 0.3), nrow = 2, byrow = TRUE)

  out <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(10, 10),
    cooksCutoff = 5,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = TRUE
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(0.01, 0.02),
    maxCooks = c(10, 10),
    cooksCutoff = 5,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = TRUE,
    native = TRUE
  )

  testthat::expect_equal(out$cooksOutlier, c(FALSE, TRUE))
  testthat::expect_equal(out$pvalue, c(0.01, NA))
  testthat::expect_equal(nativeOut, out)
})

testthat::test_that("R wrapper Cook's cutoff helper validates primitive inputs", {
  testthat::expect_error(applyCooksCutoffRust(c(-0.1), c(1), 5))
  testthat::expect_error(applyCooksCutoffRust(c(0.1), c(Inf), 5))
  testthat::expect_error(applyCooksCutoffRust(c(0.1), c(1), TRUE))
  testthat::expect_error(
    applyCooksCutoffRust(
      c(0.1),
      c(1),
      5,
      counts = matrix(1, nrow = 1, ncol = 2),
      cooks = matrix(1, nrow = 1, ncol = 1),
      lowCountHeuristic = TRUE
    )
  )
})
