test_that("saved locfit and translated custom curves retain independent frozen VST precision", {
  skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  means <- exp(seq(log(.3), log(5000), length.out = 60L))
  dispersions <- .03 + 1.5 / means + .015 * (1 + sin(log(means)))
  saved <- DESeq2:::localDispersionFit(means, dispersions, 1e-8)
  attr(saved, "fitType") <- "local"
  path <- tempfile(fileext = ".rds")
  saveRDS(saved, path)
  saved <- readRDS(path)
  unlink(path)
  custom <- local({
    asymptote <- .15
    function(mu) { asymptote + .3 / sqrt(pmax(mu, 1e-8)) }
  })
  attr(custom, "fitType") <- "local"
  floor <- function(mu) rep(1e-8, length(mu))
  attr(floor, "fitType") <- "local"
  for (trend in list(saved, custom, floor)) for (matrix_factors in c(FALSE, TRUE)) {
    dds <- DESeq2::DESeqDataSetFromMatrix(counts,
      data.frame(condition = factor(rep(c("A", "B"), each = 4L)), row.names = colnames(counts)),
      ~ condition)
    if (matrix_factors) {
      DESeq2::normalizationFactors(dds) <- outer(seq(.8, 1.2, length.out = nrow(counts)), seq(.7, 1.4, length.out = ncol(counts)))
    } else DESeq2::sizeFactors(dds) <- seq(.7, 1.4, length.out = ncol(counts))
    DESeq2::dispersionFunction(dds) <- trend
    expected <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
    actual <- rsdeseq2::varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
    expect_equal(SummarizedExperiment::assay(actual), SummarizedExperiment::assay(expected), tolerance = 2e-6)
  }
})

test_that("frozen function translation rejects executable bindings before execution", {
  touched <- FALSE
  functions <- list(
    function(mu) { touched <<- TRUE; rep(.2, length(mu)) },
    function(mu) stop("must never execute"),
    local({ hidden <- function(x) { touched <<- TRUE; x }; function(mu) hidden(mu) })
  )
  for (fun in functions) expect_error(rsdeseq2:::.rsdeseq2_frozen_local_model(fun), "unsupported|one expression")
  env <- new.env(parent = baseenv())
  makeActiveBinding("constant", function() { touched <<- TRUE; .2 }, env)
  active <- function(mu) rep(constant, length(mu))
  environment(active) <- env
  expect_error(rsdeseq2:::.rsdeseq2_frozen_local_model(active), "active")
  rm("constant", envir = env)
  delayedAssign("constant", { touched <<- TRUE; .2 }, assign.env = env)
  expect_error(rsdeseq2:::.rsdeseq2_frozen_local_model(active), "delayed")
  expect_false(touched)
})

test_that("parallel extrema retain original numeric missingness and recycling", {
  operations <- list(pmax = base::pmax, pmin = base::pmin)
  for (name in names(operations)) for (remove in c(FALSE, TRUE)) {
    x <- c(NA_real_, NaN, -Inf, Inf, -1, 0, 1)
    y <- rev(x)
    call <- list(type = "call", function_name = name, arguments = list(
      list(name = NULL, expression = list(type = "symbol", name = "x")),
      list(name = NULL, expression = list(type = "symbol", name = "y")),
      list(name = "na.rm", expression = list(type = "literal", value = remove))))
    actual <- rsdeseq2:::.rsdeseq2_formula_request(list(operation = "evaluate",
      expressions = list(call), variables = list(x = x, y = y), collation = "root"))[[1L]]
    expected <- operations[[name]](x, y, na.rm = remove)
    expect_identical(.row_reference_bits(actual), .row_reference_bits(expected))
  }
})

test_that("saved-fit imports reject altered predictors and transforms without execution", {
  skip_if_not_installed("DESeq2")
  means <- exp(seq(log(.3),log(5000),length.out=60L))
  curve <- DESeq2:::localDispersionFit(means,.03+1.5/means,1e-8)
  original <- get("fit",environment(curve),inherits=FALSE)
  env <- new.env(parent=environment(curve))
  environment(curve) <- env
  env$fit <- original
  touched <- FALSE
  altered <- original
  altered$trans <- function(x) { touched <<- TRUE; x }
  env$fit <- altered
  expect_error(rsdeseq2:::.rsdeseq2_frozen_local_model(curve),"identity transform")
  env$fit <- original
  makeActiveBinding("predict",function() { touched <<- TRUE; stats::predict },env)
  expect_error(rsdeseq2:::.rsdeseq2_frozen_local_model(curve),"masked saved local")
  expect_false(touched)
})
