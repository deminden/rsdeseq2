testthat::test_that("custom likelihood import records syntax and scalar captures without calls", {
  scale <- .7
  shift <- .2
  likelihood <- list(name = "custom_logistic",
    lcdfFUN = function(x) -log1p(exp(-(x - shift) / scale)),
    lpdfFUN = function(x) -(x - shift) / scale - 2 * log1p(exp(-(x - shift) / scale)) - log(scale))
  parsed <- rsdeseq2:::.rsdeseq2_restricted_likelihood(likelihood)
  testthat::expect_identical(parsed$name, "custom_logistic")
  testthat::expect_identical(parsed$log_cdf$arguments, "x")
  testthat::expect_identical(parsed$log_cdf$constants, list(shift = shift, scale = scale))
  testthat::expect_identical(parsed$log_cdf$expression$type, "call")
  testthat::expect_null(parsed$first)
  testthat::expect_null(parsed$second)
  likelihood$etruncFUN <- function(lower, upper) (lower + upper) / 2
  parsed <- rsdeseq2:::.rsdeseq2_restricted_likelihood(likelihood)
  testthat::expect_identical(parsed$first$arguments, c("lower", "upper"))
})

testthat::test_that("custom likelihood import rejects executable bindings and unused callbacks", {
  touched <- FALSE
  danger <- function(x) { touched <<- TRUE; stop("must not execute") }
  inspect <- function(fun) rsdeseq2:::.rsdeseq2_restricted_likelihood_function(fun, 1L, "custom likelihood")
  for (fun in list(function(x) danger(x), function(x) ifelse(TRUE, x, danger(x)),
                   function(x) { x + 1; x * 2 }, function(x) sum(x))) {
    testthat::expect_error(inspect(fun), "unsupported|one expression")
  }
  env <- new.env(parent = baseenv())
  fun <- function(x) x + constant
  environment(fun) <- env
  makeActiveBinding("constant", function() { touched <<- TRUE; .2 }, env)
  testthat::expect_error(inspect(fun), "active")
  rm("constant", envir = env)
  delayedAssign("constant", { touched <<- TRUE; .2 }, assign.env = env)
  testthat::expect_error(inspect(fun), "delayed")
  rm("constant", envir = env)
  env$constant <- c(1, 2)
  testthat::expect_error(inspect(fun), "scalars")
  env$exp <- function(x) { touched <<- TRUE; x }
  masked <- function(x) exp(x)
  environment(masked) <- env
  testthat::expect_error(inspect(masked), "custom or unsupported")
  testthat::expect_false(touched)
})
