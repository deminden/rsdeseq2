#!/usr/bin/env Rscript
# Original ashr only: preserve recycling warnings and output-dependent order.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion()) == "4.6.2",
          as.character(packageVersion("ashr")) == "2.2.63", !file.exists(args[[1L]]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), ".."))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
bits <- function(v) lapply(as.numeric(v), function(x) paste0("bits:", hex(x)))
expected <- function(v) if (is.null(v)) NULL else lapply(v, function(x)
  if (is.na(x) && !is.nan(x)) "NA" else sprintf("%.17g", x))
x <- c(-2, -.4, 0, .8, 2, 4)
s <- c(.3, .5, 1, 1.2, .7, .4)
g <- ashr::unimix(c(.2, .3, .5), c(0, -.4, -3), c(0, .4, 3))
outputs <- list(zero=0, one=1, two=2, three=3, four=4, five=5,
  mean="PosteriorMean", sd="PosteriorSD", negative="NegativeProb", positive="PositiveProb",
  lfsr="lfsr", svalue="svalue", lfdr="lfdr", qvalue="qvalue", likelihood="loglik",
  ratio="logLR", flash="flash_data", prior="fitted_g", empty=character(),
  duplicate=c("PosteriorMean", "PosteriorMean", "fitted_g"),
  selected=c("PosteriorSD", "svalue", "qvalue", "fitted_g"))
settings <- lapply(names(outputs), function(name)
  list(name=paste0("nonmultiple_", name), output=outputs[[name]]))
settings <- c(settings, list(
  list(name="nonmultiple_fitted", output=3, fitted=TRUE),
  list(name="nonmultiple_missing", output=3, missing=TRUE),
  list(name="nonmultiple_long_weights", output=3, long=TRUE),
  list(name="multiple_control", output=3, multiple=TRUE),
  list(name="equal_total_control", output=3, equal=TRUE)))
cases <- lapply(settings, function(p) {
  y <- x
  if (isTRUE(p$missing)) y[[2L]] <- NA_real_
  weights <- if (isTRUE(p$long)) c(.4, .3, .2, .1) else c(.7, .3)
  scales <- cbind(seq(.3, 1, length.out=6L), seq(.5, 2, length.out=6L), seq(.7, 1.5, length.out=6L))
  if (isTRUE(p$multiple)) scales <- cbind(scales, seq(1, 3, length.out=6L))
  if (isTRUE(p$equal)) {
    weights <- c(.1, .2, .1, .15, .25, .2)
    scales <- matrix(c(.3, .5, .7, 1, 1.5, 2), nrow=2L)
  }
  lik <- ashr::lik_normalmix(weights, scales)
  options <- list(mixcompdist="uniform", pointmass=TRUE, prior="uniform",
    mixsd=c(.4, 3), lik=lik, outputlevel=p$output)
  fixed <- !isTRUE(p$fitted)
  if (fixed) { options$g <- g; options$fixg <- TRUE } else options$mode <- 0
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(do.call(ashr::ash, c(list(betahat=y, sebetahat=s), options)),
    warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error=identity)
  supplied <- if (!fixed) NULL else list(component_class="unimix", weights=bits(g$pi), means=NULL,
    sds=NULL, lower=bits(g$a), upper=bits(g$b))
  location <- list(name=p$name, builtin_likelihood=list(type="normalmix", weights=bits(weights),
    scales=.formula_descriptor(scales)), estimates=bits(y), standard_errors=bits(s), family="uniform",
    point_mass=TRUE, nullweight=NULL, scales=bits(c(.4, 3)), mode=bits(0)[[1L]], range=bits(c(-Inf, Inf)),
    alpha=bits(0)[[1L]], prune_threshold=bits(1e-10)[[1L]], observation_weights=NULL,
    fixed_prior=fixed, supplied_prior=supplied, df=NULL, thresholds=list(), optimizer="mixSQP", control=list())
  list(name=p$name, estimates=bits(y), standard_errors=bits(s), output=.formula_descriptor(p$output),
    threshold=bits(0)[[1L]], query_svalues=TRUE, location=location,
    expected_error=if (inherits(fit, "error")) conditionMessage(fit) else NULL,
    posterior_mean=if (inherits(fit, "error")) NULL else expected(fit$result$PosteriorMean),
    posterior_sd=if (inherits(fit, "error")) NULL else expected(fit$result$PosteriorSD),
    svalues=if (inherits(fit, "error")) NULL else expected(fit$result$svalue),
    fitted_prior=if (inherits(fit, "error")) NULL else !is.null(fit$fitted_g), warnings=as.list(warnings))
})
# Reuse only primitive inputs/syntax from the earlier independent fixture.
# Every expected output and warning below is freshly evaluated by original ashr.
previous <- jsonlite::fromJSON(file.path(root, "crates/rsdeseq2/tests/data/custom_likelihood_r462.json"), simplifyVector=FALSE)
laplace <- previous$cases[[which(vapply(previous$cases, function(c) c$name == "laplace", logical(1)))]]
input_fields <- c("custom_likelihood", "estimates", "standard_errors", "family", "point_mass", "nullweight", "scales", "mode",
                  "range", "alpha", "prune_threshold", "observation_weights", "fixed_prior", "supplied_prior", "df")
restricted <- list(name="custom_laplace",
  lcdfFUN=function(x) ifelse(x < 0, log(.5) + x/.8, log1p(-.5*exp(-x/.8))),
  lpdfFUN=function(x) -abs(x)/.8-log(2*.8))
laplace_x <- c(-4, -1, 0, .3, 2, 5)
laplace_s <- c(.3, .7, 1, 2, .4, 1)
for (name in c("zero", "one", "two", "four", "five", "mean", "sd", "likelihood")) {
  output <- outputs[[name]]
  warnings <- character()
  fit <- withCallingHandlers(ashr::ash(laplace_x, laplace_s, g=g, fixg=TRUE, lik=restricted, outputlevel=output),
    warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  location <- laplace[input_fields]
  location$name <- paste0("restricted_", name)
  location$thresholds <- list()
  cases[[length(cases)+1L]] <- list(name=location$name, estimates=bits(laplace_x), standard_errors=bits(laplace_s),
    output=.formula_descriptor(output), threshold=bits(0)[[1L]], query_svalues=TRUE, location=location,
    expected_error=NULL, posterior_mean=expected(fit$result$PosteriorMean), posterior_sd=expected(fit$result$PosteriorSD),
    svalues=expected(fit$result$svalue), fitted_prior=!is.null(fit$fitted_g), warnings=as.list(warnings))
}
# Original likelihood-only mode searches also repeat the selected optimum.
.warning_mode_capture <- new.env()
suppressMessages(trace("ash.estmode", where=asNamespace("ashr"), print=FALSE,
  tracer=quote({.warning_mode_capture$bounds <- c(modemin,modemax)}),
  exit=quote({.warning_mode_capture$mode <- returnValue()})))
suppressMessages(trace("ash.workhorse", where=asNamespace("ashr"), print=FALSE, exit=quote({
  if (identical(outputlevel,"loglik")) .warning_mode_capture$trials <- c(.warning_mode_capture$trials,list(c(mode,returnValue()$loglik)))
})))
mode_cases <- lapply(c(FALSE,TRUE),function(restricted_mode) {
  location <- if(restricted_mode) laplace[input_fields] else cases[[1L]]$location
  location$name <- if(restricted_mode) "restricted_mode" else "recycled_mode"
  location["nullweight"] <- list(NULL)
  location$fixed_prior <- FALSE
  location$supplied_prior <- NULL
  location$thresholds <- list()
  location$scales <- bits(c(.4,3))
  y <- if(restricted_mode)laplace_x else x
  error <- if(restricted_mode)laplace_s else s
  lik <- if(restricted_mode)restricted else ashr::lik_normalmix(c(.7,.3),
    cbind(seq(.3,1,length.out=6),seq(.5,2,length.out=6),seq(.7,1.5,length.out=6)))
  .warning_mode_capture$trials <- list()
  warnings <- character()
  fit <- withCallingHandlers(ashr::ash(y,error,mode="estimate",mixcompdist="uniform",pointmass=TRUE,
      prior="uniform",mixsd=c(.4,3),lik=lik,outputlevel=3),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  c(location,list(mode_input=.formula_descriptor("estimate"),expected_error=NULL,native_expected_error=NULL,
    expected_mode=expected(.warning_mode_capture$mode)[[1L]], search_bounds=bits(.warning_mode_capture$bounds),
    trials=lapply(.warning_mode_capture$trials,expected),posterior_mean=expected(fit$result$PosteriorMean),
    posterior_sd=expected(fit$result$PosteriorSD),svalues=expected(fit$result$svalue),
    log_likelihood=expected(fit$loglik)[[1L]],warnings=as.list(warnings)))
})
suppressMessages(untrace("ash.workhorse",where=asNamespace("ashr")))
suppressMessages(untrace("ash.estmode",where=asNamespace("ashr")))
jsonlite::write_json(list(source=list(R=as.character(getRversion()), ashr=as.character(packageVersion("ashr"))),
  cases=cases, mode_cases=mode_cases), args[[1L]], auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
