#!/usr/bin/env Rscript
# Original-package query policy, independent of the native implementation.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2',
  as.character(packageVersion('ashr')) == '2.2.63', !file.exists(args[[1L]]))
bits <- function(x) lapply(as.numeric(x), function(v)
  paste0('bits:', paste(format(writeBin(v, raw(), size=8, endian='big')), collapse='')))
expected <- function(x) lapply(as.numeric(x), function(v)
  if (is.na(v) && !is.nan(v)) 'NA' else sprintf('%.17g',v))
g <- ashr::normalmix(c(.1,.2,.3,.4), c(0,-.4,.8,2), c(0,.3,2,1))
u <- ashr::unimix(c(.2,.3,.5), c(0,-1,-4), c(0,1,4))
x <- c(-2,-.3,0,.3,2,NA_real_,NaN,.6,-.6)
s <- c(.2,.7,1,.7,.2,1,1,0,Inf)
profiles <- list(
  no_svalues=list(), sign_svalues=list(svalue=TRUE), threshold=list(threshold=.4),
  threshold_missing=list(x=x,s=s,threshold=.4),
  alpha_fit_normal_query=list(alpha=.5,threshold=.4),
  forwarded_fit=list(fit_x=c(-1,.2,2),fit_s=c(.6,.8,.4),threshold=.4),
  forwarded_sign=list(fit_x=c(-1,.2,2),fit_s=c(.6,.8,.4),svalue=TRUE),
  t_fit_normal_query=list(g=u,df=5,threshold=.4),
  unnormalized=list(g=ashr::normalmix(c(1,2,3,4),g$mean,g$sd),x=x,s=s,threshold=.4),
  undefined_uniform=list(g=u,x=c(-1e6,0,1e6),s=rep(.001,3),threshold=.4))
cases <- lapply(names(profiles), function(name) {
  p <- profiles[[name]]
  original_x <- if(is.null(p$x)) c(-2,-.3,0,.3,2,.6) else p$x
  original_s <- if(is.null(p[['s']])) c(.2,.7,1,.7,.2,.8) else p[['s']]
  fit_x <- if(is.null(p$fit_x)) original_x else p$fit_x
  fit_s <- if(is.null(p$fit_s)) original_s else p$fit_s
  prior <- if(is.null(p$g)) g else p$g
  alpha <- if(is.null(p$alpha)) 0 else p$alpha
  threshold <- if(is.null(p$threshold)) 0 else p$threshold
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr,warning=function(w){
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
  fit <- capture(ashr::ash(fit_x,fit_s,g=prior,fixg=TRUE,alpha=alpha,
    df=if(is.null(p$df)) Inf else p$df, mixcompdist=if(is.null(p$df)) 'normal' else 'uniform'))
  # Assigning to the original result population follows data.frame recycling.
  result <- data.frame(original=original_x)
  result$mean <- fit$result$PosteriorMean
  result$sd <- fit$result$PosteriorSD
  sv <- NULL
  if(threshold > 0) {
    data <- ashr::set_data(original_x,original_s)
    positive <- capture(ashr::cdf_post(fit$fitted_g,threshold,data))
    negative <- capture(ashr::cdf_post(fit$fitted_g,-threshold,data))
    rate <- ifelse(result$mean > 0,positive,1-negative)
    sv <- cumsum(sort(rate,na.last=TRUE))/seq_along(rate)
    sv <- sv[rank(rate,ties.method='first',na.last=TRUE)]
  } else if(isTRUE(p$svalue)) {
    result$svalue <- fit$result$svalue; sv <- result$svalue
  }
  list(name=name,estimates=bits(fit_x),standard_errors=bits(fit_s),
    query_estimates=bits(original_x),query_standard_errors=bits(original_s),
    prior=list(component_class=class(prior)[[1L]],weights=bits(prior$pi),
      means=if(is.null(prior$mean)) NULL else bits(prior$mean),
      sds=if(is.null(prior$sd)) NULL else bits(prior$sd),
      lower=if(is.null(prior$a)) NULL else bits(prior$a),
      upper=if(is.null(prior$b)) NULL else bits(prior$b)),
    alpha=bits(alpha)[[1L]],df=if(is.null(p$df)) NULL else bits(p$df)[[1L]],
    threshold=bits(threshold)[[1L]],svalue=isTRUE(p$svalue),
    posterior_mean=expected(result$mean),posterior_sd=expected(result$sd),
    svalues=if(is.null(sv)) NULL else expected(sv),warnings=as.list(warnings))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),
  ashr=as.character(packageVersion('ashr')),scope='original wrapper query and recycling policy'),
  cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
