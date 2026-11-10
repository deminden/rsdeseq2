#!/usr/bin/env Rscript
# Independent complete normal/t likelihood batches; never load rsdeseq2.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_location_shrinkage_reference.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
bits <- function(x)lapply(as.numeric(x),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
expected <- function(x)lapply(as.numeric(x),function(v)if(is.na(v)&&!is.nan(v))'NA' else sprintf('%.17g',v))
components <- function(g) list(component_class=class(g)[[1L]],weights=bits(g$pi),means=if(is.null(g$mean))NULL else bits(g$mean),sds=if(is.null(g$sd))NULL else bits(g$sd),lower=if(is.null(g$a))NULL else bits(g$a),upper=if(is.null(g$b))NULL else bits(g$b))
i <- seq_len(96);x <- 1.4*sin(i*.73)+.6*cos(i*.11);s <- .2+(i%%11)/13
normal <- ashr::normalmix(c(.1,.3,.6),c(0,-.4,.8),c(0,.3,2))
uniform <- ashr::unimix(c(.1,.3,.6),c(0,-.5,-4),c(0,.5,4))
truncated <- ashr::tnormalmix(c(.1,.3,.6),rep(0,3),c(0,.4,2),c(-Inf,-Inf,0),c(0,0,Inf))
missing_x <- c(-2,-.3,0,.3,2,NA_real_,NaN,.6,-.6)
missing_s <- c(.2,.7,1,.7,.2,1,1,0,Inf)
profiles <- list(
  automatic=list(x=x,s=s,arguments=list(mixcompdist='normal',method='shrink')),
  automatic_null=list(x=x,s=s,arguments=list(mixcompdist='normal',pointmass=TRUE,prior='nullbiased')),
  scalar_se=list(x=x[1:24],s=.7,arguments=list(mixcompdist='normal',method='shrink')),
  alpha_half=list(x=x,s=s,arguments=list(mixcompdist='normal',method='shrink',alpha=.5)),
  alpha_one=list(x=x[1:24],s=s[1:24],arguments=list(g=normal,fixg=TRUE,alpha=1)),
  weighted=list(x=x,s=s,arguments=list(mixcompdist='normal',method='shrink',weights=as.numeric(i%%7))),
  fixed_normal_missing=list(x=missing_x,s=missing_s,arguments=list(g=normal,fixg=TRUE)),
  fixed_normal_missing_non_coprime=list(x=missing_x,s=missing_s,arguments=list(g=ashr::normalmix(c(.1,.2,.3,.4),c(0,-.4,.8,2),c(0,.3,2,1)),fixg=TRUE)),
  fixed_uniform_missing=list(x=missing_x,s=missing_s,arguments=list(g=uniform,fixg=TRUE)),
  fixed_halfnormal_missing=list(x=missing_x,s=missing_s,arguments=list(g=truncated,fixg=TRUE)),
  fixed_pruned=list(x=missing_x,s=missing_s,arguments=list(g=normal,fixg=TRUE,pi_thresh=.2)),
  supplied_refit=list(x=x,s=s,arguments=list(g=normal)),
  fixed_unnormalized=list(x=missing_x,s=missing_s,arguments=list(g=ashr::normalmix(c(1,2,3),c(0,-.4,.8),c(0,.3,2)),fixg=TRUE)),
  uniform_fit=list(x=x[1:24],s=s[1:24],arguments=list(mixcompdist='uniform',mixsd=c(.2,1,4),pointmass=TRUE,prior='nullbiased')),
  halfuniform_fit=list(x=x[1:24],s=s[1:24],arguments=list(mixcompdist='halfuniform',mixsd=c(.2,1,4),pointmass=TRUE,prior='nullbiased')),
  halfnormal_fit=list(x=x[1:24],s=s[1:24],arguments=list(mixcompdist='halfnormal',mixsd=c(.2,1,4),pointmass=TRUE,prior='nullbiased')),
  positive_uniform=list(x=abs(x[1:24]),s=s[1:24],arguments=list(mixcompdist='+uniform',mixsd=c(.2,1,4),pointmass=TRUE,prior='nullbiased')),
  negative_uniform=list(x=-abs(x[1:24]),s=s[1:24],arguments=list(mixcompdist='-uniform',mixsd=c(.2,1,4),pointmass=TRUE,prior='nullbiased')),
  student_uniform=list(x=x[1:24],s=s[1:24],arguments=list(g=uniform,fixg=TRUE,df=5)),
  all_excluded=list(x=c(-1,0,1),s=c(Inf,0,Inf),arguments=list(g=ashr::normalmix(1,0,1),fixg=TRUE)))
cutoffs <- c(-1,0,1)
cases <- lapply(names(profiles),function(name){
 p <- profiles[[name]];a <- p$arguments;warnings <- character()
 record <- list(name=name,estimates=bits(p$x),standard_errors=bits(p$s),thresholds=bits(cutoffs),
  family=if(is.null(a$mixcompdist))'normal'else a$mixcompdist,
  point_mass=if(identical(a$method,'shrink'))FALSE else if(is.null(a$pointmass))TRUE else a$pointmass,
  nullweight=if(!identical(a$method,'shrink')&&(is.null(a$prior)||a$prior=='nullbiased'))bits(10)[[1L]]else NULL,
  scales=if(is.null(a$mixsd))NULL else bits(a$mixsd),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),
  alpha=bits(if(is.null(a$alpha))0 else a$alpha)[[1L]],
  prune_threshold=bits(if(is.null(a$pi_thresh))1e-10 else a$pi_thresh)[[1L]],
  observation_weights=if(is.null(a$weights))NULL else bits(a$weights),
  fixed_prior=isTRUE(a$fixg),supplied_prior=if(is.null(a$g))NULL else components(a$g),
  df=if(is.null(a$df))NULL else bits(a$df)[[1L]])
 result <- tryCatch({
 fit <- withCallingHandlers(do.call(ashr::ash,c(list(betahat=p$x,sebetahat=p$s,outputlevel=3),a)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 data <- ashr::set_data(p$x,p$s,alpha=if(is.null(a$alpha))0 else a$alpha,lik=if(is.null(a$df))NULL else ashr::lik_t(a$df))
 cdf <- withCallingHandlers(t(vapply(cutoffs,function(c)ashr::cdf_post(fit$fitted_g,c,data),numeric(length(p$x)))),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 list(
  fitted_prior=components(fit$fitted_g),
  posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD),
  negative_probability=expected(fit$result$NegativeProb),positive_probability=expected(fit$result$PositiveProb),
  local_false_discovery_rate=expected(fit$result$lfdr),local_false_sign_rate=expected(fit$result$lfsr),
  svalues=expected(fit$result$svalue),qvalues=expected(fit$result$qvalue),
  cdf=expected(t(cdf)),log_likelihood=expected(fit$loglik)[[1L]],warnings=as.list(warnings))
 },error=function(e)list(expected_error=conditionMessage(e),warnings=as.list(warnings)))
 c(record,result)
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),
  input_encoding='f64_bits',scope='complete original-package fit and summary batches'),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
