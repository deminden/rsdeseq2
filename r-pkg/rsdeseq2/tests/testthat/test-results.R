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
  testthat::expect_equal(shared$resultName, "condition_X2_vs_X3")
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
  testthat::expect_equal(resultsNamesRust(list(resultNames = colnames(beta))), colnames(beta))
  testthat::expect_equal(resultsNamesRust(list(coefNames = colnames(beta))), colnames(beta))
  testthat::expect_equal(resultsNamesRust(list(coefNames = colnames(beta), coefficients = beta)), colnames(beta))
  testthat::expect_equal(resultsNamesRust(list(coefficients = beta)), colnames(beta))
  formulaBeta <- beta
  colnames(formulaBeta) <- c("Intercept", "cell type_as.double", "condition_C_vs_A")
  shuffledBeta <- formulaBeta[, c("condition_C_vs_A", "Intercept", "cell type_as.double")]
  colnames(shuffledBeta) <- c("condition_C_vs_A", "Intercept", "cell.type_as.double")
  shuffledFit <- waldFitRust(
    c(10, 20),
    beta = shuffledBeta,
    resultsNames = colnames(formulaBeta),
    betaSE = matrix(1, nrow = 2, ncol = 3)
  )
  testthat::expect_equal(colnames(shuffledFit$beta), colnames(formulaBeta))
  testthat::expect_equal(unname(shuffledFit$beta), unname(formulaBeta))
  resolved <- resolveResultsContrastRust(
    c("condition", "C", "B"),
    resultsNames(wald),
    reference = "A"
  )
  testthat::expect_equal(unname(resolved$numeric), c(0, -1, 1))
  testthat::expect_error(resultsNamesRust(list(beta = matrix(1, nrow = 1))))
  testthat::expect_error(
    resultsNamesRust(list(coefNames = c("Intercept", "condition_B_vs_A"), coefficients = beta)),
    "resultsNames must have one value per beta column"
  )
  testthat::expect_error(
    resultsNamesRust(c("Intercept", "condition_B_vs_A", "condition_B_vs_A")),
    "resultsNames must be unique"
  )
  testthat::expect_error(
    waldFitRust(
      c(10, 20),
      beta = matrix(
        c(1, 2, 3, 2, 4, 7),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_B_vs_A"))
      )
    ),
    "resultsNames must be unique"
  )
  duplicateBetaColumns <- beta
  colnames(duplicateBetaColumns) <- c("Intercept", "condition.B.vs.A", "condition.B.vs.A")
  testthat::expect_error(
    waldFitRust(
      c(10, 20),
      beta = duplicateBetaColumns,
      resultsNames = c("Intercept", "condition_B_vs_A", "condition_C_vs_A")
    ),
    "beta column names must be unique when supplied"
  )
  testthat::expect_error(resultsNamesRust(NULL))
})

testthat::test_that("R wrapper resultsNames aliases preserve formula-derived coefficient names", {
  beta <- matrix(
    c(
      1, 2,
      1, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_as.double"))
  )
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    coefficients = beta,
    coefNames = c("Intercept", "cell type_as.double"),
    betaSE = matrix(1, nrow = 2, ncol = 2)
  )

  out <- results(object, name = "cell.type_as.double", independentFiltering = FALSE)

  testthat::expect_equal(resultsNamesRust(object), c("Intercept", "cell type_as.double"))
  testthat::expect_equal(resultsNames(object), c("Intercept", "cell type_as.double"))
  testthat::expect_equal(out$log2FoldChange, c(2, 3))
  testthat::expect_equal(attr(out, "resultName"), "cell type_as.double")
  testthat::expect_equal(attr(out, "contrast")[["cell type_as.double"]], 1)
})

testthat::test_that("R wrapper resolves R-cleaned results contrast aliases", {
  names <- c("(Intercept)", "if.", "condition.B.1", "cell type_as.double")

  coefList <- resolveResultsContrastRust(list("if", "condition-B 1"), names)
  numericTransform <- resolveResultsContrastRust(list("cell.type_as.double"), names)
  expanded <- resolveResultsContrastRust(c("condition", "B-1", "A-1"), c("conditionA.1", "conditionB.1"))
  sharedNoSeparator <- resolveResultsContrastRust(c("condition", "C", "B"), c("Intercept", "conditionB_vs_A", "conditionC_vs_A"))
  exactWins <- resolveResultsContrastRust(list("cell.type_as.double"), c("Intercept", "cell.type_as.double", "cell type_as.double"))

  testthat::expect_equal(unname(coefList$numeric), c(0, 1, -1, 0))
  testthat::expect_equal(unname(numericTransform$numeric), c(0, 0, 0, 1))
  testthat::expect_equal(unname(exactWins$numeric), c(0, 1, 0))
  testthat::expect_equal(unname(expanded$numeric), c(-1, 1))
  testthat::expect_equal(unname(sharedNoSeparator$numeric), c(0, -1, 1))
  testthat::expect_error(
    resolveResultsContrastRust(list("cell.type_as.double"), c("Intercept", "cell type_as.double", "cell-type_as.double")),
    "results name alias cell.type_as.double is ambiguous"
  )
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

testthat::test_that("R object extraction coerces DataFrame-like mcols before column lookup", {
  mcols <- matrix(
    c(
      10, 1, 2, 0.1, 0.5, 4, 0.01,
      20, 1, -2, 0.1, 0.5, 3, 0.02
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("g1", "g2"),
      c("baseMean", "Intercept", "condition_B_vs_A", "SE_Intercept", "SE_condition_B_vs_A", "stat", "pvalue")
    )
  )
  object <- list(
    mcols = mcols,
    counts = matrix(
      c(0, 0, 10, 0, 10, 0, 0, 5),
      nrow = 2,
      byrow = TRUE
    ),
    sampleLevels = c("A", "A", "B", "B")
  )

  wald <- results(
    object,
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
  testthat::expect_equal(rownames(wald), c("g1", "g2"))
  testthat::expect_equal(lrt$stat, c(4, 3))
  testthat::expect_equal(lrt$pvalue, c(0.01, 0.02))
})

testthat::test_that("R object extraction uses DESeq2-shaped model matrix attributes", {
  beta <- matrix(
    c(
      1, 2, 5,
      1, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "batch_Y_vs_X"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  counts <- matrix(
    c(
      5, 6, 0, 0,
      5, 6, 2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  modelMatrix <- matrix(
    c(
      1, 1, 0,
      1, 0, 1,
      1, 0, 0,
      1, 0, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("s3", "s4", "s1", "s2"), colnames(beta))
  )
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    counts = counts
  )
  attr(object, "modelMatrix") <- modelMatrix

  fit <- nbinomWaldTestRust(object = object)
  out <- results(fit, contrast = c(0, 1, -1), independentFiltering = FALSE)

  testthat::expect_equal(rownames(fit$modelMatrix), colnames(counts))
  testthat::expect_equal(out$log2FoldChange, c(0, -3))
  testthat::expect_equal(out$stat[[1]], 0)
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
})

testthat::test_that("R results resolves R-cleaned mcols beta and SE aliases", {
  resultsNames <- c("(Intercept)", "condition-B vs A")
  mcols <- data.frame(
    baseMean = c(10, 20),
    "(Intercept)" = c(1, 1),
    "condition-B vs A" = c(2, -2),
    "SE_(Intercept)" = c(0.1, 0.1),
    "SE_condition-B vs A" = c(0.5, 0.5),
    row.names = c("g1", "g2")
  )
  object <- list(
    mcols = mcols,
    resultsNames = resultsNames
  )

  out <- results(
    object,
    name = "condition-B vs A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(resultsNames(object), resultsNames)
  testthat::expect_equal(out$log2FoldChange, c(2, -2))
  testthat::expect_equal(out$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(rownames(out), c("g1", "g2"))
})

testthat::test_that("R results rejects ambiguous cleaned mcols aliases", {
  object <- list(
    mcols = data.frame(
      baseMean = c(10, 20),
      "condition B" = c(2, -2),
      "condition-B" = c(3, -3),
      "SE_condition B" = c(0.5, 0.5),
      "SE_condition-B" = c(0.6, 0.6),
      check.names = FALSE
    ),
    resultsNames = "condition.B"
  )

  testthat::expect_error(
    results(object, name = "condition.B", independentFiltering = FALSE),
    "mcols column alias condition.B is ambiguous"
  )
})

testthat::test_that("R results accepts common object aliases for model matrices", {
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
      10, 12, 0, 0, 0, 0,
      10, 12, 30, 36, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
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
    byrow = TRUE,
    dimnames = list(colnames(counts), colnames(beta))
  )
  shuffledModelMatrix <- modelMatrix[c("s5", "s1", "s3", "s6", "s2", "s4"), , drop = FALSE]
  listObject <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    counts = counts,
    designMatrix = shuffledModelMatrix
  )

  listOut <- results(
    listObject,
    contrast = c(0, 1, -1),
    independentFiltering = FALSE
  )

  testthat::expect_equal(listOut$log2FoldChange[[1]], 0)
  testthat::expect_equal(listOut$pvalue[[1]], 1)
  testthat::expect_equal(attr(listOut, "contrastAllZero"), "numeric")

  if (!methods::isClass("Rsdeseq2FittedDesignAliasMock")) {
    methods::setClass(
      "Rsdeseq2FittedDesignAliasMock",
      slots = c(
        mcols = "data.frame",
        betaCovariance = "array",
        counts = "matrix",
        fullModelMatrix = "matrix"
      )
    )
  }
  mcols <- data.frame(
    baseMean = c(10, 20),
    Intercept = c(1, 1),
    condition_B_vs_A = c(2, 2),
    condition_C_vs_A = c(3, 3),
    SE_Intercept = c(0.1, 0.1),
    SE_condition_B_vs_A = c(0.5, 0.5),
    SE_condition_C_vs_A = c(0.5, 0.5),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )
  s4Object <- methods::new(
    "Rsdeseq2FittedDesignAliasMock",
    mcols = mcols,
    betaCovariance = cov,
    counts = counts,
    fullModelMatrix = shuffledModelMatrix
  )

  s4Out <- results(
    s4Object,
    contrast = c(0, 1, -1),
    independentFiltering = FALSE
  )

  testthat::expect_equal(s4Out$log2FoldChange[[1]], 0)
  testthat::expect_equal(s4Out$pvalue[[1]], 1)
  testthat::expect_equal(attr(s4Out, "contrastAllZero"), "numeric")

  badObject <- listObject
  rownames(badObject$designMatrix)[[1L]] <- "unknown"
  testthat::expect_error(
    results(badObject, contrast = c(0, 1, -1), independentFiltering = FALSE),
    "modelMatrix row names must contain all count sample names"
  )
})

testthat::test_that("R results accepts common object aliases for coefficient metadata", {
  beta <- matrix(
    c(
      1, 2,
      1, -2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A"))
  )
  cov <- array(0, dim = c(2, 2, 2))
  cov[, 2, 2] <- 0.25
  listObject <- list(
    baseMean = c(g1 = 10, g2 = 20),
    coefficients = beta,
    coefNames = colnames(beta),
    coefCovariance = cov
  )

  listOut <- results(
    listObject,
    name = "condition_B_vs_A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(resultsNames(listObject), colnames(beta))
  testthat::expect_equal(listOut$log2FoldChange, c(2, -2))
  testthat::expect_equal(listOut$lfcSE, c(0.5, 0.5))

  if (!methods::isClass("Rsdeseq2FittedCoefAliasMock")) {
    methods::setClass(
      "Rsdeseq2FittedCoefAliasMock",
      slots = c(
        baseMean = "numeric",
        coef = "matrix",
        resultNames = "character",
        betaCov = "array"
      )
    )
  }
  s4Object <- methods::new(
    "Rsdeseq2FittedCoefAliasMock",
    baseMean = c(g1 = 10, g2 = 20),
    coef = beta,
    resultNames = colnames(beta),
    betaCov = cov
  )

  s4Out <- results(
    s4Object,
    name = "condition_B_vs_A",
    independentFiltering = FALSE
  )

  testthat::expect_equal(resultsNames(s4Object), colnames(beta))
  testthat::expect_equal(s4Out$log2FoldChange, c(2, -2))
  testthat::expect_equal(s4Out$lfcSE, c(0.5, 0.5))
})

testthat::test_that("R results accepts common object aliases for factor levels", {
  beta <- matrix(
    c(
      1, 2, 3,
      1, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "condition_B_vs_A", "condition_C_vs_A"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  listObject <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    factorLevelNames = list(condition = c("A", "B", "C"))
  )

  listOut <- results(
    listObject,
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(listOut$log2FoldChange, c(1, 3))
  testthat::expect_equal(attr(listOut, "contrast")[["condition_B_vs_A"]], -1)
  testthat::expect_equal(attr(listOut, "contrast")[["condition_C_vs_A"]], 1)

  if (!methods::isClass("Rsdeseq2FittedFactorLevelAliasMock")) {
    methods::setClass(
      "Rsdeseq2FittedFactorLevelAliasMock",
      slots = c(
        baseMean = "numeric",
        beta = "matrix",
        betaCovariance = "array",
        levelNames = "list"
      )
    )
  }
  s4Object <- methods::new(
    "Rsdeseq2FittedFactorLevelAliasMock",
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    levelNames = list(condition = c("A", "B", "C"))
  )

  s4Out <- results(
    s4Object,
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(s4Out$log2FoldChange, c(1, 3))
  testthat::expect_equal(attr(s4Out, "contrast")[["condition_B_vs_A"]], -1)
  testthat::expect_equal(attr(s4Out, "contrast")[["condition_C_vs_A"]], 1)
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

  fromCharacterColData <- results(
    list(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta[, c("Intercept", "condition_B_vs_A", "condition_C_vs_A")],
      betaCovariance = cov[, c(1, 2, 4), c(1, 2, 4)],
      counts = matrix(
        c(
          10, 12, 20, 24, 50, 60,
          10, 12, 20, 24, 50, 60
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("g1", "g2"), paste0("s", 1:6))
      ),
      colData = data.frame(condition = c("C", "A", "B", "C", "A", "B"))
    ),
    contrast = c("condition", "A", "B"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(fromCharacterColData$log2FoldChange, c(-2, -4))
  testthat::expect_equal(attr(fromCharacterColData, "contrast")[["condition_B_vs_A"]], -1)

  relevelFit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    factorLevels = list(condition = c("B", "A", "C")),
    factorReferences = c(condition = "A")
  )
  relevelOut <- results(
    relevelFit,
    contrast = c("condition", "C", "B"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(relevelOut$log2FoldChange, c(1, 3))
  testthat::expect_equal(attr(relevelOut, "contrast")[["condition_B_vs_A"]], -1)
  testthat::expect_equal(attr(relevelOut, "contrast")[["condition_C_vs_A"]], 1)

  cleanedReferenceFit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = matrix(
      c(
        1, 2, 3,
        1, 4, 7
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
    ),
    betaCovariance = cov[, c(1, 2, 4), c(1, 2, 4)],
    factorLevels = list("cell type" = c("B cell", "A cell", "C cell")),
    factorReferences = c("cell.type" = "A.cell")
  )
  cleanedReferenceOut <- results(
    cleanedReferenceFit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    independentFiltering = FALSE
  )
  testthat::expect_equal(cleanedReferenceOut$log2FoldChange, c(1, 3))
  testthat::expect_equal(attr(cleanedReferenceOut, "resultName"), "cell.type_C.cell_vs_B.cell")
  testthat::expect_equal(attr(cleanedReferenceOut, "comparison"), "factor-level contrast: cell type C cell vs B cell")
  testthat::expect_equal(cleanedReferenceFit$factorReferences[["cell type"]], "A cell")
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      factorLevels = list("cell type" = c("A", "B"), "cell-type" = c("A", "B")),
      factorReferences = c("cell.type" = "A")
    ),
    "factorLevels names must resolve to unique R-cleaned aliases"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      factorLevels = list(condition = c("A", "B", "B")),
      factorReferences = c(condition = "A")
    ),
    "factorLevels entries must contain unique levels"
  )
  testthat::expect_error(
    lrtFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      lrtStat = c(g1 = 5, g2 = 7),
      lrtPvalue = c(g1 = 0.03, g2 = 0.02),
      factorLevels = list(condition = c("A", "B", "B")),
      factorReferences = c(condition = "A")
    ),
    "factorLevels entries must contain unique levels"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      factorLevels = list("cell type" = c("A cell", "A-cell", "B cell")),
      factorReferences = c("cell.type" = "B.cell")
    ),
    "factorLevels entries must resolve to unique R-cleaned level aliases"
  )
  testthat::expect_error(
    lrtFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      lrtStat = c(g1 = 5, g2 = 7),
      lrtPvalue = c(g1 = 0.03, g2 = 0.02),
      factorLevels = list("cell type" = c("A cell", "A-cell", "B cell")),
      factorReferences = c("cell.type" = "B.cell")
    ),
    "factorLevels entries must resolve to unique R-cleaned level aliases"
  )
  testthat::expect_error(
    waldFitRust(
      baseMean = c(g1 = 10, g2 = 20),
      beta = beta,
      betaCovariance = cov,
      factorLevels = list("cell type" = c("A", "B")),
      factorReferences = c("cell type" = "A", "cell.type" = "A")
    ),
    "factorReferences aliases must resolve to unique factor names"
  )
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

testthat::test_that("R results accepts S4 assay containers with bracket access", {
  if (!methods::isClass("Rsdeseq2AssaysMock")) {
    methods::setClass("Rsdeseq2AssaysMock", slots = c(data = "list"))
  }
  methods::setMethod(
    "[[",
    signature(x = "Rsdeseq2AssaysMock", i = "ANY", j = "missing"),
    function(x, i, j, ..., exact = TRUE) {
      x@data[[i]]
    }
  )
  if (!methods::isClass("Rsdeseq2FittedNestedAssaysMock")) {
    methods::setClass(
      "Rsdeseq2FittedNestedAssaysMock",
      slots = c(
        mcols = "data.frame",
        assays = "ANY",
        colData = "data.frame"
      )
    )
  }
  mcols <- data.frame(
    baseMean = c(10, 20),
    Intercept = c(1, 1),
    condition_B_vs_A = c(2, -2),
    condition_C_vs_A = c(3, 1),
    SE_Intercept = c(0.1, 0.1),
    SE_condition_B_vs_A = c(0.5, 0.5),
    SE_condition_C_vs_A = c(0.5, 0.5),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )
  counts <- matrix(
    c(
      0, 0, 0, 0, 50, 60,
      10, 12, 20, 24, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  object <- methods::new(
    "Rsdeseq2FittedNestedAssaysMock",
    mcols = mcols,
    assays = methods::new("Rsdeseq2AssaysMock", data = list(counts = counts)),
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
  testthat::expect_equal(out$log2FoldChange[[2]], -2)
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
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

testthat::test_that("R wrapper object list contrasts accept cleaned formula coefficient aliases", {
  beta <- matrix(
    c(
      1, 2, 3,
      1, 4, 7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 1
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(6, 8),
    lrtPvalue = c(0.02, 0.01)
  )

  wald <- results(
    nbinomWaldTestRust(object = object),
    contrast = list("cell.type_C.cell_vs_A.cell", "cell.type_B.cell_vs_A.cell"),
    independentFiltering = FALSE
  )
  lrt <- results(
    nbinomLRTRust(object = object),
    contrast = list("cell.type_C.cell_vs_A.cell", "cell.type_B.cell_vs_A.cell"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(wald$log2FoldChange, c(1, 3))
  testthat::expect_equal(wald$lfcSE, rep(sqrt(1.25), 2))
  testthat::expect_equal(attr(wald, "resultName"), "contrast")
  testthat::expect_equal(
    attr(wald, "comparison"),
    "coefficient list contrast: cell.type_C.cell_vs_A.cell vs cell.type_B.cell_vs_A.cell"
  )
  testthat::expect_equal(attr(wald, "contrast")[["cell type_C cell_vs_A cell"]], 1)
  testthat::expect_equal(attr(wald, "contrast")[["cell type_B cell_vs_A cell"]], -1)
  testthat::expect_equal(lrt$log2FoldChange, c(1, 3))
  testthat::expect_equal(lrt$stat, c(6, 8))
  testthat::expect_equal(lrt$pvalue, c(0.02, 0.01))
  testthat::expect_equal(
    attr(lrt, "comparison"),
    "coefficient list contrast: cell.type_C.cell_vs_A.cell vs cell.type_B.cell_vs_A.cell; LRT full vs reduced"
  )
  testthat::expect_equal(attr(lrt, "contrast")[["cell type_C cell_vs_A cell"]], 1)
  testthat::expect_equal(attr(lrt, "contrast")[["cell type_B cell_vs_A cell"]], -1)
})

testthat::test_that("R wrapper object results accept cleaned formula numeric coefficient aliases", {
  beta <- matrix(
    c(
      1, 2.5,
      1, 4.5
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_as.double"))
  )
  cov <- array(0, dim = c(2, 2, 2))
  cov[, 2, 2] <- 0.25
  object <- list(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(5, 7),
    lrtPvalue = c(0.03, 0.02)
  )

  wald <- results(
    nbinomWaldTestRust(object = object),
    name = "cell.type_as.double",
    independentFiltering = FALSE
  )
  lrt <- results(
    nbinomLRTRust(object = object),
    name = "cell.type_as.double",
    independentFiltering = FALSE
  )

  testthat::expect_equal(wald$log2FoldChange, c(2.5, 4.5))
  testthat::expect_equal(wald$lfcSE, c(0.5, 0.5))
  testthat::expect_equal(attr(wald, "resultName"), "cell type_as.double")
  testthat::expect_equal(attr(wald, "comparison"), "coefficient cell type_as.double")
  testthat::expect_equal(attr(wald, "contrast")[["cell type_as.double"]], 1)
  testthat::expect_equal(lrt$log2FoldChange, c(2.5, 4.5))
  testthat::expect_equal(lrt$stat, c(5, 7))
  testthat::expect_equal(lrt$pvalue, c(0.03, 0.02))
  testthat::expect_equal(attr(lrt, "comparison"), "coefficient cell type_as.double; LRT full vs reduced")
  testthat::expect_equal(attr(lrt, "contrast")[["cell type_as.double"]], 1)
})

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

testthat::test_that("R primitive constructors canonicalize factor references through colData levels", {
  beta <- matrix(
    c(
      1, 2, 3,
      1, 2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
  )
  cov <- array(0, dim = c(2, 3, 3))
  cov[, 2, 2] <- 0.25
  cov[, 3, 3] <- 0.25
  counts <- matrix(
    c(
      0, 0, 50, 60, 0, 0,
      10, 12, 30, 36, 50, 60
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))
  )
  colData <- data.frame(check.names = FALSE, row.names = colnames(counts))
  colData[["cell type"]] <- factor(
    c("A cell", "A cell", "B cell", "B cell", "C cell", "C cell"),
    levels = c("B cell", "A cell", "C cell")
  )

  waldFit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    counts = counts,
    colData = colData,
    factorReferences = c("cell.type" = "A.cell")
  )
  lrtFit <- lrtFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
    lrtStat = c(g1 = 5, g2 = 7),
    lrtPvalue = c(g1 = 0.03, g2 = 0.02),
    counts = counts,
    colData = colData,
    factorReferences = c("cell.type" = "A.cell")
  )
  wald <- results(
    waldFit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    independentFiltering = FALSE
  )
  lrt <- results(
    lrtFit,
    contrast = c("cell.type", "C.cell", "B.cell"),
    independentFiltering = FALSE
  )

  testthat::expect_equal(waldFit$factorLevels[["cell type"]], c("B cell", "A cell", "C cell"))
  testthat::expect_equal(waldFit$factorReferences[["cell type"]], "A cell")
  testthat::expect_equal(lrtFit$factorReferences[["cell type"]], "A cell")
  testthat::expect_equal(attr(wald, "contrast")[["cell type_C cell_vs_A cell"]], 1)
  testthat::expect_equal(attr(wald, "contrast")[["cell type_B cell_vs_A cell"]], -1)
  testthat::expect_equal(attr(wald, "comparison"), "factor-level contrast: cell type C cell vs B cell")
  testthat::expect_equal(attr(lrt, "comparison"), "factor-level contrast: cell type C cell vs B cell; LRT full vs reduced")
  testthat::expect_equal(attr(wald, "contrastAllZero"), "character")
  testthat::expect_equal(attr(lrt, "contrastAllZero"), "character")
})

testthat::test_that("R results resolves ordered formula-transform character contrasts", {
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
  fit <- waldFitRust(
    baseMean = c(g1 = 10, g2 = 20),
    beta = beta,
    betaCovariance = cov,
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
  testthat::expect_equal(out$pvalue[[1]], 1)
  testthat::expect_true(out$pvalue[[2]] < 1)
  testthat::expect_equal(attr(out, "resultName"), paste0(transformedAlias, "_A_vs_B"))
  testthat::expect_equal(
    attr(out, "comparison"),
    paste0("factor-level contrast: ", transformed, " A vs B")
  )
  testthat::expect_equal(attr(out, "contrastAllZero"), "character")
  testthat::expect_equal(attr(out, "contrast")[[paste0(transformed, "_A_vs_B")]], 1)
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
    "sampleLevels or colData to contain numerator level B"
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
    log2FoldChange = c(g3 = NA, g1 = 1, g2 = -0.5),
    lfcSE = c(g2 = 0.3, g3 = NA, g1 = 0.2),
    stat = c(g1 = 5, g3 = NA, g2 = -1.5),
    pvalue = c(g2 = 0.05, g1 = 0.001, g3 = NA),
    dispersion = c(g3 = 0.3, g2 = 0.2, g1 = 0.1),
    converged = c(g2 = TRUE, g1 = TRUE, g3 = FALSE),
    contrast = c(Intercept = 0, condition_B_vs_A = 1)
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(rownames(out), c("g1", "g2", "g3"))
  testthat::expect_equal(
    colnames(out),
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged")
  )
  testthat::expect_equal(out$log2FoldChange, c(1, -0.5, NA))
  testthat::expect_equal(out$lfcSE, c(0.2, 0.3, NA))
  testthat::expect_equal(out$stat, c(5, -1.5, NA))
  testthat::expect_equal(out$dispersion, c(0.1, 0.2, 0.3))
  testthat::expect_equal(out$converged, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(attr(out, "contrast"), c(Intercept = 0, condition_B_vs_A = 1))
  testthat::expect_equal(out$padj, stats::p.adjust(c(0.001, 0.05, NA), method = "BH"))
})

testthat::test_that("R wrapper result table helper validates result vectors", {
  testthat::expect_error(resultsTableRust(c(1, NA)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1)))
  testthat::expect_error(resultsTableRust(c(1, 2), pvalue = c(0.1, 2)))
  testthat::expect_error(resultsTableRust(c(1, 2), padj = c(0.1, 0.2)))
  testthat::expect_error(resultsTableRust(c(1, 2), converged = c(TRUE)))
  testthat::expect_error(resultsTableRust(c(1, 2), contrast = numeric()))
  testthat::expect_error(resultsTableRust(c(1, 2), contrast = c(0, NA)))
  testthat::expect_error(
    resultsTableRust(
      c(g1 = 1, g2 = 2),
      pvalue = c(g1 = 0.1, gx = 0.2)
    ),
    "pvalue names must contain all beta row names"
  )

  out <- resultsTableRust(c(1, 2), rowNames = c("a", "b"))
  testthat::expect_equal(rownames(out), c("a", "b"))
  testthat::expect_true(all(is.na(out$pvalue)))
  testthat::expect_true(all(is.na(out$padj)))
})

testthat::test_that("R wrapper applies primitive independent filtering", {
  table <- resultsTableRust(
    baseMean = c(g1 = 0, g2 = 1, g3 = 2),
    stat = c(1, 1, 1),
    pvalue = c(0.01, 0.02, 0.03),
    contrast = c(Intercept = 0, condition_B_vs_A = 1)
  )
  attr(table, "resultName") <- "condition_B_vs_A"
  attr(table, "comparison") <- "coefficient condition_B_vs_A"
  attr(table, "contrastAllZero") <- "none"

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
  testthat::expect_equal(attr(filtered, "resultName"), "condition_B_vs_A")
  testthat::expect_equal(attr(filtered, "comparison"), "coefficient condition_B_vs_A")
  testthat::expect_equal(attr(filtered, "contrast"), c(Intercept = 0, condition_B_vs_A = 1))
  testthat::expect_equal(attr(filtered, "contrastAllZero"), "none")
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
  counts <- matrix(
    c(1, 5, 6, 7, 9, 5, 6, 7),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  cooks <- matrix(
    c(10, 0.1, 0.2, 0.3, 10, 0.1, 0.2, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )

  out <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02),
    maxCooks = c(g2 = 10, g1 = 10),
    cooksCutoff = 5,
    counts = counts[c("g2", "g1"), c("s3", "s1", "s4", "s2")],
    cooks = cooks[c("g2", "g1"), c("s2", "s3", "s1", "s4")],
    lowCountHeuristic = TRUE
  )
  nativeOut <- applyCooksCutoffRust(
    pvalue = c(g1 = 0.01, g2 = 0.02),
    maxCooks = c(g2 = 10, g1 = 10),
    cooksCutoff = 5,
    counts = counts[c("g2", "g1"), c("s3", "s1", "s4", "s2")],
    cooks = cooks[c("g2", "g1"), c("s2", "s3", "s1", "s4")],
    lowCountHeuristic = TRUE,
    native = TRUE
  )

  testthat::expect_equal(rownames(out), c("g1", "g2"))
  testthat::expect_equal(out$cooksOutlier, c(FALSE, TRUE))
  testthat::expect_equal(out$pvalue, c(0.01, NA))
  testthat::expect_equal(nativeOut, out)

  duplicateCooks <- cooks
  colnames(duplicateCooks) <- c("s1", "s2", "s2", "s4")
  testthat::expect_error(
    applyCooksCutoffRust(
      pvalue = c(g1 = 0.01, g2 = 0.02),
      maxCooks = c(10, 10),
      cooksCutoff = 5,
      counts = counts,
      cooks = duplicateCooks,
      lowCountHeuristic = TRUE
    ),
    "cooks column names must be unique"
  )
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
