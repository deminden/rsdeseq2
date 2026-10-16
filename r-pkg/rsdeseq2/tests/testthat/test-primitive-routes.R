testthat::test_that("R primitive argument packaging preserves Wald and chi-square LRT answers", {
  object <- .rsdeseq2_primitive_fixture(
    c(1, 2, 2, 4), c("Intercept", "condition_B_vs_A"), c(0.01, 0.25)
  )
  object$betaCovariance <- NULL
  object$betaSE <- matrix(c(0.1, 0.5, 0.1, 0.5), nrow = 2, byrow = TRUE)
  for (method in c("Wald", "LRT")) {
    args <- object
    if (method == "LRT") args <- c(args, list(lrtStat = c(6, 8), lrtDf = 1))
    constructor <- if (method == "Wald") nbinomWaldTestRust else nbinomLRTRust
    fit <- do.call(constructor, args)
    out <- results(fit, name = "condition_B_vs_A", independentFiltering = FALSE)
    if (method == "Wald") {
      .rsdeseq2_expect_result(out, c(2, 4), c(0.5, 0.5), info = method)
    } else {
      .rsdeseq2_expect_result(out, c(2, 4), c(0.5, 0.5), c(6, 8),
        stats::pchisq(c(6, 8), df = 1, lower.tail = FALSE), info = method)
    }
    testthat::expect_s3_class(fit,
      if (method == "Wald") "rsdeseq2PrimitiveWaldFit" else "rsdeseq2PrimitiveLrtFit")
    comparison <- paste0("coefficient condition_B_vs_A",
      if (method == "LRT") "; LRT full vs reduced" else "")
    testthat::expect_equal(attr(out, "comparison"), comparison, info = method)
    filtered <- results(fit, name = "condition_B_vs_A")
    testthat::expect_equal(attr(filtered, "resultName"), "condition_B_vs_A", info = method)
    testthat::expect_equal(attr(filtered, "comparison"), comparison, info = method)
    testthat::expect_equal(attr(filtered, "contrast"), c(Intercept = 0, condition_B_vs_A = 1), info = method)
    testthat::expect_true(attr(filtered, "independentFiltering")$enabled, info = method)
  }
})

testthat::test_that("R primitive list dispatch and packaging preserve contrast answers and precedence", {
  object <- .rsdeseq2_primitive_fixture()
  for (method in c("Wald", "LRT")) {
    input <- object
    if (method == "LRT") input <- c(input, list(lrtStat = c(6, 8), lrtPvalue = c(0.02, 0.01)))
    constructor <- if (method == "Wald") nbinomWaldTestRust else nbinomLRTRust
    # Raw dispatch and constructor extraction must both work.
    for (route in c("direct", "packaged")) {
      fit <- if (route == "direct") input else constructor(object = input)
      for (with_name in c(FALSE, TRUE)) {
        label <- paste(method, route, with_name)
        out <- results(fit, contrast = c("condition", "C", "B"), reference = "A",
          name = if (with_name) "condition_B_vs_A" else NULL, independentFiltering = FALSE)
        if (method == "Wald") {
          .rsdeseq2_expect_result(out, c(1, 3), rep(sqrt(1.25), 2), info = label)
        } else {
          .rsdeseq2_expect_result(out, c(1, 3), rep(sqrt(1.25), 2), c(6, 8), c(0.02, 0.01), info = label)
        }
        testthat::expect_equal(attr(out, "resultName"), "condition_C_vs_B", info = label)
        testthat::expect_equal(attr(out, "contrastAllZero"), "character", info = label)
        testthat::expect_equal(attr(out, "contrast"),
          c(Intercept = 0, condition_B_vs_A = -1, condition_C_vs_A = 1), info = label)
      }
    }
  }
})

testthat::test_that("R cleaned formula selectors retain original coefficient names and method metadata", {
  factor_object <- .rsdeseq2_primitive_fixture(c(1, 2, 3, 1, 4, 7),
    c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell"))
  numeric_object <- .rsdeseq2_primitive_fixture(c(1, 2.5, 1, 4.5),
    c("Intercept", "cell type_as.double"), c(0, 0.25))
  cases <- list(
    factor_list = list(object = factor_object,
      selector = list(contrast = list("cell.type_C.cell_vs_A.cell", "cell.type_B.cell_vs_A.cell")),
      lfc = c(1, 3), se = rep(sqrt(1.25), 2), stat = c(6, 8), pvalue = c(0.02, 0.01),
      result_name = "contrast", weights = c(0, -1, 1),
      comparison = "coefficient list contrast: cell.type_C.cell_vs_A.cell vs cell.type_B.cell_vs_A.cell"),
    numeric_name = list(object = numeric_object, selector = list(name = "cell.type_as.double"),
      lfc = c(2.5, 4.5), se = c(0.5, 0.5), stat = c(5, 7), pvalue = c(0.03, 0.02),
      result_name = "cell type_as.double", weights = c(0, 1), comparison = "coefficient cell type_as.double")
  )
  for (case in names(cases)) for (method in c("Wald", "LRT")) {
    input <- cases[[case]]
    # Retain extraction of a Wald fit from an object that also has LRT fields.
    object <- c(input$object, list(lrtStat = input$stat, lrtPvalue = input$pvalue))
    constructor <- if (method == "Wald") nbinomWaldTestRust else nbinomLRTRust
    fit <- constructor(object = object)
    out <- do.call(results, c(list(object = fit, independentFiltering = FALSE), input$selector))
    label <- paste(case, method)
    if (method == "Wald") {
      .rsdeseq2_expect_result(out, input$lfc, input$se, info = label)
    } else {
      .rsdeseq2_expect_result(out, input$lfc, input$se, input$stat, input$pvalue, info = label)
    }
    testthat::expect_equal(attr(out, "resultName"), input$result_name, info = label)
    testthat::expect_equal(attr(out, "comparison"), paste0(input$comparison,
      if (method == "LRT") "; LRT full vs reduced" else ""), info = label)
    testthat::expect_equal(attr(out, "contrast"),
      stats::setNames(input$weights, colnames(input$object$beta)), info = label)
  }
})
