#!/usr/bin/env Rscript
# Original ashr fits whose mixsqp subproblems have numerically singular
# Hessians: every standard error is large next to the default normal grid, so
# the columns of the likelihood matrix are nearly identical. The fitted prior
# is then a single vertex, and a wrong search direction ends at another one.
# Base R, ashr and jsonlite only.
#
# Usage: Rscript --vanilla generate_ashr_singular_direction_reference.R FRESH_OUTPUT_JSON
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2',
  as.character(packageVersion('ashr')) == '2.2.63',
  as.character(packageVersion('mixsqp')) == '0.3.54', !file.exists(args[[1L]]))
bits <- function(x) lapply(as.numeric(x), function(v)
  paste0('bits:', paste(format(writeBin(v, raw(), size = 8, endian = 'big')), collapse = '')))
expected <- function(x) lapply(as.numeric(x), function(v) sprintf('%.17g', v))
inputs <- list(
  golden_panel_eight_genes = list(
    b = c(0.8938459130692038, -0.30807406620808614, 0.5188649759446969,
      0.6867711762458827, -1.0225709949451243, 0.7285431855079626,
      0.78147041493471, 0.4468318106955514),
    s = c(4.047762889349763, 4.133025421302511, 2.9889064674996098,
      2.839866842552506, 3.4874172938544956, 2.566068424736632,
      4.09496696651059, 2.579278191601985)))
set.seed(20261125)
for (i in 1:3) {
  n <- c(6L, 12L, 25L)[[i]]
  inputs[[paste0('wide_errors_', n)]] <- list(b = rnorm(n, 0, 0.8), s = runif(n, 2, 5))
}
cases <- lapply(names(inputs), function(name) {
  b <- inputs[[name]]$b; s <- inputs[[name]]$s
  warnings <- character()
  fit <- withCallingHandlers(
    ashr::ash(b, s, mixcompdist = 'normal', method = 'shrink', outputlevel = 3),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    })
  g <- ashr::get_fitted_g(fit)
  list(name = name, estimates = bits(b), standard_errors = bits(s),
    log_likelihood = sprintf('%.17g', fit$loglik),
    fitted_weights = expected(g$pi), grid_sds = expected(g$sd),
    posterior_mean = expected(fit$result$PosteriorMean),
    posterior_sd = expected(fit$result$PosteriorSD),
    warnings = as.list(warnings))
})
jsonlite::write_json(list(source = list(R = as.character(getRversion()),
  ashr = as.character(packageVersion('ashr')),
  mixsqp = as.character(packageVersion('mixsqp')),
  call = 'ash(b, s, mixcompdist = "normal", method = "shrink")',
  scope = 'mixsqp fits with numerically singular reduced Hessians'),
  cases = cases), args[[1L]], auto_unbox = TRUE, pretty = TRUE, null = 'null')
