# Character, numeric, list, and named result selectors.

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
