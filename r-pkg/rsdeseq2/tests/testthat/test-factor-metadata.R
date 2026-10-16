testthat::test_that("R factor metadata routes preserve Wald tails and supplied LRT tests", {
  transformed <- 'ordered(condition, levels = c("B", "A"))'
  ordered <- .rsdeseq2_primitive_fixture(
    c(1, 2, 1, 2), c("Intercept", paste0(transformed, "_A_vs_B")), c(0, 0.25)
  )
  counts <- .rsdeseq2_count_fixture(c(0, 0, 0, 0, 50, 60, 10, 12, 20, 24, 50, 60))
  ordered$counts <- counts
  ordered$colData <- data.frame(row.names = colnames(counts))
  ordered$colData[[transformed]] <- factor(c("A", "A", "B", "B", "C", "C"), levels = c("B", "A", "C"))
  cleaned <- .rsdeseq2_primitive_fixture(
    coefficients = c("Intercept", "cell type_B cell_vs_A cell", "cell type_C cell_vs_A cell")
  )
  cleaned$counts <- counts
  cleaned$colData <- data.frame(row.names = colnames(counts))
  cleaned$colData[["cell type"]] <- factor(
    c("A cell", "A cell", "B cell", "B cell", "C cell", "C cell"), levels = c("A cell", "B cell", "C cell")
  )
  references <- .rsdeseq2_primitive_fixture(coefficients = colnames(cleaned$beta))
  references$factorLevels <- list("cell type" = c("B cell", "A cell", "C cell"))
  cases <- list(
    transformed_object = list(object = ordered, route = "object", factor = transformed,
      levels = c("B", "A", "C"), contrast = c(make.names(transformed), "A", "B"),
      lfc = c(0, 2), se = c(0.5, 0.5), weights = c(0, 1),
      result_name = paste0(make.names(transformed), "_A_vs_B"),
      comparison = paste0("factor-level contrast: ", transformed, " A vs B")),
    cleaned_reference = list(object = cleaned, route = "object", factor = "cell type",
      levels = c("A cell", "B cell", "C cell"), contrast = c("cell.type", "C.cell", "B.cell"),
      reference = "A.cell", lfc = c(1, 3), se = rep(sqrt(1.25), 2), weights = c(0, -1, 1),
      result_name = "cell.type_C.cell_vs_B.cell", comparison = "factor-level contrast: cell type C cell vs B cell"),
    reference_aliases = list(object = references, route = "object", factor = "cell type",
      levels = c("B cell", "A cell", "C cell"), contrast = c("cell.type", "C.cell", "B.cell"),
      lfc = c(1, 3), se = rep(sqrt(1.25), 2), weights = c(0, -1, 1),
      result_name = "cell.type_C.cell_vs_B.cell", comparison = "factor-level contrast: cell type C cell vs B cell")
  )
  # Direct constructors and object packaging are distinct public routes.
  cases$transformed_constructor <- cases$transformed_object
  cases$transformed_constructor$route <- "constructor"
  for (case in names(cases)) {
    input <- cases[[case]]
    for (method in c("Wald", "LRT")) {
      label <- paste(case, method)
      object <- input$object
      lfc <- input$lfc
      if (case == "reference_aliases") {
        # Retain both alias spellings and the original method-specific inputs.
        alias <- if (method == "Wald") "referenceLevels" else "factorReferenceLevels"
        object[[alias]] <- c("cell.type" = "A.cell")
        if (method == "LRT") {
          object$beta[2, ] <- c(1, 2, 3)
          lfc <- c(1, 1)
        }
      }
      if (method == "LRT") {
        object$lrtStat <- c(6, 8)
        object$lrtPvalue <- c(0.02, 0.01)
      }
      if (input$route == "constructor") {
        object$factorLevels <- setNames(list(input$levels), input$factor)
        constructor <- if (method == "Wald") waldFitRust else lrtFitRust
        fit <- do.call(constructor, object)
      } else {
        constructor <- if (method == "Wald") nbinomWaldTestRust else nbinomLRTRust
        fit <- constructor(object = object)
      }
      out <- results(fit, contrast = input$contrast, reference = input$reference,
                     independentFiltering = FALSE)
      if (method == "Wald") {
        .rsdeseq2_expect_result(out, lfc, input$se, info = label)
      } else {
        .rsdeseq2_expect_result(out, lfc, input$se, c(6, 8), c(0.02, 0.01), info = label)
      }
      testthat::expect_s3_class(fit, if (method == "Wald") "rsdeseq2PrimitiveWaldFit" else "rsdeseq2PrimitiveLrtFit")
      testthat::expect_equal(fit$factorLevels[[input$factor]], input$levels, info = label)
      testthat::expect_equal(attr(out, "resultName"), input$result_name, info = label)
      comparison <- paste0(input$comparison, if (method == "LRT") "; LRT full vs reduced" else "")
      testthat::expect_equal(attr(out, "comparison"), comparison, info = label)
      testthat::expect_equal(
        attr(out, "contrast"), stats::setNames(input$weights, colnames(object$beta)), info = label
      )
      if (case == "reference_aliases") {
        testthat::expect_equal(fit$factorReferences[["cell type"]], "A cell", info = label)
      } else {
        testthat::expect_equal(attr(out, "contrastAllZero"), "character", info = label)
      }
    }
  }
})
