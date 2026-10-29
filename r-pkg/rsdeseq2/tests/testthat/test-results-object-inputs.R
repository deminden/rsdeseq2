# Result queries through matrix, list, and dataset inputs.

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
