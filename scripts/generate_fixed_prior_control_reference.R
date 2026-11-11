#!/usr/bin/env Rscript
# Original ashr only: fixed priors do not execute optimizer control handling.
args <- commandArgs(TRUE)
if(length(args)!=1L) stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion())=="4.6.2",as.character(packageVersion("ashr"))=="2.2.63",!file.exists(args[[1L]]))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
bits <- function(v) lapply(as.numeric(v),function(x)paste0("bits:",hex(x)))
expected <- function(v) if(is.null(v))NULL else lapply(v,function(x)if(is.na(x)&&!is.nan(x))"NA"else sprintf("%.17g",x))
x <- c(-2,-.4,0,.8,2,4);se <- c(.3,.5,1,1.2,.7,.4)
g <- ashr::normalmix(c(.5,.3,.2),rep(0,3),c(0,.5,2))
methods <- c("mixSQP","mixIP","mixEM","w_mixEM","mixVBEM","cxxMixSquarem")
settings <- list()
for(method in methods) {
  bad <- if(method=="mixSQP")list(maxiter.sqp=-1)else if(method=="mixIP")list(rtol=-1)else list(maxiter=-1)
  for(control in list(list(unused="ignored"),"not a list",bad))
    settings[[length(settings)+1L]] <- list(method=method,control=control,fixed=TRUE)
  settings[[length(settings)+1L]] <- list(method=method,control=list(unused="invalid when fitted"),fixed=FALSE)
}
settings <- c(settings,list(list(method="not_an_optimizer",control=list(),fixed=TRUE),
                           list(method="mixSQP",control="unused",fixed=TRUE,missing_g=TRUE)))
cases <- lapply(seq_along(settings),function(i) {
  p <- settings[[i]];warnings <- character();messages <- character()
  arguments <- list(betahat=x,sebetahat=se,mixcompdist="normal",method="shrink",fixg=p$fixed,
                    optmethod=p$method,control=p$control,outputlevel=3)
  if(!isTRUE(p$missing_g))arguments$g <- g
  fit <- tryCatch(withCallingHandlers(do.call(ashr::ash,arguments),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},
    message=function(m){messages <<- c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),error=identity)
  prior <- if(isTRUE(p$missing_g))NULL else list(component_class="normalmix",weights=bits(g$pi),means=bits(g$mean),
                                                sds=bits(g$sd),lower=NULL,upper=NULL)
  name <- paste0(p$method,"_",i)
  location <- list(name=name,estimates=bits(x),standard_errors=bits(se),thresholds=list(),family="normal",point_mass=FALSE,
    nullweight=NULL,scales=NULL,mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(0)[[1L]],prune_threshold=bits(1e-10)[[1L]],
    observation_weights=NULL,fixed_prior=p$fixed,supplied_prior=prior,df=NULL,optimizer=p$method,control=p$control)
  list(name=name,estimates=bits(x),standard_errors=bits(se),output=.formula_descriptor(3),threshold=bits(0)[[1L]],query_svalues=TRUE,
    location=location,expected_error=if(inherits(fit,"error"))conditionMessage(fit)else NULL,
    posterior_mean=if(inherits(fit,"error"))NULL else expected(fit$result$PosteriorMean),
    posterior_sd=if(inherits(fit,"error"))NULL else expected(fit$result$PosteriorSD),
    svalues=if(inherits(fit,"error"))NULL else expected(fit$result$svalue),
    fitted_prior=if(inherits(fit,"error"))NULL else !is.null(fit$fitted_g),warnings=as.list(warnings),original_messages=as.list(messages))
})
jsonlite::write_json(list(scope="fixed_prior_controls",source=list(R=as.character(getRversion()),ashr=as.character(packageVersion("ashr"))),cases=cases),
                    args[[1L]],auto_unbox=TRUE,pretty=TRUE,null="null",digits=NA)
