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
testthat::test_that("already forced closure parameters are read without repeating their expressions", {
  testthat::skip_if(getRversion() < "4.6.0")
  state <- new.env(parent = emptyenv())
  state$touched <- FALSE
  make_function <- function(constant) {
    force(constant)
    function(x) x + constant
  }
  fun <- make_function({ state$touched <- TRUE; .25 })
  testthat::expect_true(state$touched)
  state$touched <- FALSE
  parsed <- rsdeseq2:::.rsdeseq2_restricted_likelihood_function(fun, 1L, "custom likelihood")
  testthat::expect_identical(parsed$constants, list(constant = .25))
  testthat::expect_false(state$touched)
  testthat::expect_error(.Call("rsdeseq2_stored_binding", NA_character_, environment(fun), "test", PACKAGE="rsdeseq2"), "requires")
})

testthat::test_that("canonical likelihood descriptors retain parameter captures and reject masked bodies", {
  testthat::skip_if_not_installed("ashr")
  translate <- rsdeseq2:::.rsdeseq2_native_ash_likelihood
  testthat::expect_identical(translate(ashr::lik_normal()), list(type = "normal"))
  testthat::skip_if(getRversion() < "4.6.0")
  testthat::expect_identical(translate(ashr::lik_t(5)), list(type = "t", df = 5))
  normal <- ashr::lik_normal()
  normal$etruncFUN <- function(a, b) stop("must not execute")
  testthat::expect_error(translate(normal), "unsupported")
})

testthat::test_that("regression likelihood syntax supports vector densities without widening formulas", {
  translate <- rsdeseq2:::.rsdeseq2_apeglm_likelihood
  sigma <- .7
  parsed <- translate(function(y,x,beta,param,offset)
    stats::dnorm(y, mean=as.numeric(x %*% beta)+offset, sd=sigma, log=TRUE))
  testthat::expect_identical(parsed$type, "restricted")
  testthat::expect_identical(parsed$model$arguments, c("y","x","beta","param","offset"))
  testthat::expect_identical(parsed$model$constants, list(sigma=sigma))
  testthat::expect_error(rsdeseq2:::.rsdeseq2_validate_formula_calls(~stats::dnorm(x)), "unsupported")
  testthat::expect_error(rsdeseq2:::.rsdeseq2_validate_formula_calls(~x %*% y), "unsupported")
  if (requireNamespace("apeglm", quietly=TRUE)) {
    testthat::expect_identical(translate(apeglm::logLikNB), list(type="negative_binomial"))
    # Identical multi-statement text outside the original namespace is still
    # a custom function and must satisfy the single-expression restriction.
    fake <- apeglm::logLikNB
    environment(fake) <- baseenv()
    testthat::expect_error(translate(fake), "one expression")
  }
})

testthat::test_that("regression imports reject callbacks, masked densities and executable captures", {
  translate <- rsdeseq2:::.rsdeseq2_apeglm_likelihood
  touched <- FALSE
  danger <- function(...) {touched <<- TRUE; stop("must not execute")}
  testthat::expect_error(translate(function(y,x,beta,param,offset)
    ifelse(TRUE, stats::dnorm(y, log=TRUE), danger(y))), "unsupported")
  env <- new.env(parent=baseenv())
  env$dnorm <- danger
  fun <- function(y,x,beta,param,offset)dnorm(y, log=TRUE)
  environment(fun) <- env
  testthat::expect_error(translate(fun), "custom or unsupported")
  fun <- function(y,x,beta,param,offset)stats::dnorm(y, sd=sigma, log=TRUE)
  environment(fun) <- env
  makeActiveBinding("sigma", function() {touched <<- TRUE; .7}, env)
  testthat::expect_error(translate(fun), "active")
  rm("sigma", envir=env)
  delayedAssign("sigma", {touched <<- TRUE; .7}, assign.env=env)
  testthat::expect_error(translate(fun), "delayed")
  testthat::expect_false(touched)
})
