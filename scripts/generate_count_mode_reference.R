#!/usr/bin/env Rscript
# Original-package mode bounds/traces, independent of candidate results.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2', as.character(packageVersion('ashr')) == '2.2.63', !file.exists(args[[1L]]))
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~', ' ', script, fixed=TRUE)), '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits <- function(v) lapply(as.numeric(v), function(x) paste0('bits:', hex(x)))
expected <- function(v) lapply(as.numeric(v), function(x) if (is.na(x) && !is.nan(x)) 'NA' else sprintf('%.17g', x))
decode <- function(v) vapply(v, function(h) readBin(as.raw(strtoi(substring(sub('bits:', '', h), seq(1,15,2), seq(2,16,2)),16L)), 'double', n=1L, size=8L, endian='big'), numeric(1))
reference <- jsonlite::fromJSON(file.path(root, 'crates/rsdeseq2/tests/data/count_likelihood_r462.json'), simplifyVector=FALSE)$cases
.count_mode_capture <- new.env()
suppressMessages(trace('ash.estmode', where=asNamespace('ashr'), print=FALSE,
  tracer=quote({.count_mode_capture$bounds <- c(modemin, modemax)}), exit=quote({.count_mode_capture$mode <- returnValue()})))
suppressMessages(trace('ash.workhorse', where=asNamespace('ashr'), print=FALSE, exit=quote({
  if (identical(outputlevel, 'loglik')) .count_mode_capture$trials <- c(.count_mode_capture$trials, list(c(mode, returnValue()$loglik)))
})))
cases <- list()
for (name in c('pois_identity_rows', 'pois_log_rows', 'binom_identity_rows', 'binom_logit_rows')) {
  for (variant in c('explicit_grid', 'automatic_grid', 'expanded_seed')) {
    c <- reference[[which(vapply(reference, function(c) identical(c$name, paste0(name, '_fitted')), logical(1)))]]
    m <- c$builtin_likelihood
    lik <- if (m$type == 'pois') ashr::lik_pois(decode(m$y), decode(m$scale), m$link) else ashr::lik_binom(decode(m$y), decode(m$trials), m$link)
    mode <- if (variant == 'expanded_seed') c(-2, 2) else 'estimate'
    scales <- if (variant == 'automatic_grid') NULL else decode(c$scales)
    .count_mode_capture$bounds <- NULL; .count_mode_capture$mode <- NULL; .count_mode_capture$trials <- list()
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(ashr::ash(decode(c$estimates), decode(c$standard_errors), mixcompdist='uniform',
      pointmass=FALSE, prior='uniform', mixsd=scales, mode=mode, lik=lik, outputlevel=3),
      warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }), error=identity)
    c$name <- paste(name, variant, sep='_'); c$mode_input <- .formula_descriptor(mode)
    c['scales'] <- list(if (is.null(scales)) NULL else bits(scales))
    c$expected_error <- if (inherits(fit, 'error')) conditionMessage(fit) else NULL
    c$expected_mode <- if (is.null(.count_mode_capture$mode)) NULL else expected(.count_mode_capture$mode)[[1L]]
    c$search_bounds <- bits(.count_mode_capture$bounds); c$trials <- lapply(.count_mode_capture$trials, expected)
    c$warnings <- as.list(warnings)
    if (!inherits(fit, 'error')) {
      c$posterior_mean <- expected(fit$result$PosteriorMean); c$posterior_sd <- expected(fit$result$PosteriorSD)
      c$svalues <- expected(fit$result$svalue); c$log_likelihood <- expected(fit$loglik)[[1L]]
    }
    cases[[length(cases)+1L]] <- c
  }
}
suppressMessages(untrace('ash.workhorse', where=asNamespace('ashr')))
suppressMessages(untrace('ash.estmode', where=asNamespace('ashr')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()), ashr=as.character(packageVersion('ashr'))), cases=cases),
  args[[1L]], auto_unbox=TRUE, pretty=TRUE, null='null', digits=NA)
