#!/usr/bin/env Rscript
# Original ashr only. Accept immutable binary input JSON, never native outputs.
args <- commandArgs(TRUE)
if(length(args)!=2L)stop('usage: INPUT_JSON FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',
  as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[2L]]))
input <- jsonlite::fromJSON(args[[1L]],simplifyVector=FALSE)
c <- if(is.null(input$location_cases)) input else input$location_cases[[1L]]
decode <- function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')},numeric(1L))
expected <- function(x)as.list(sprintf('%.17g',as.numeric(x)))
fit <- ashr::ash(betahat=decode(c$estimates),sebetahat=decode(c$standard_errors),
  mixcompdist='normal',method='shrink')
record <- list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),
  scope='positive-coordinate convergence in automatic 21-component normal prior'),
  estimates=c$estimates,standard_errors=c$standard_errors,
  fitted_weights=expected(fit$fitted_g$pi),fitted_sds=expected(fit$fitted_g$sd),
  posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD))
jsonlite::write_json(record,args[[2L]],auto_unbox=TRUE,pretty=TRUE)
