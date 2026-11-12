#!/usr/bin/env Rscript
# Original ashr only. Trace its own search bounds, trial locations and likelihoods.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
decode<-function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(sub('bits:','',h),seq(1,15,2),seq(2,16,2)),16L)),'double',n=1L,size=8,endian='big'),numeric(1))
original<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/ashr_optimizer_r462.json'),simplifyVector=FALSE)$cases[[1L]]
x<-decode(original$estimates)+.8;s<-decode(original$standard_errors)
custom<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/custom_likelihood_r462.json'),simplifyVector=FALSE)$cases[[1L]]$custom_likelihood
.mode_reference_capture<-new.env()
suppressMessages(trace('ash.estmode',where=asNamespace('ashr'),print=FALSE,
  tracer=quote({.mode_reference_capture$bounds<-c(modemin,modemax)}),exit=quote({.mode_reference_capture$mode<-returnValue()})))
suppressMessages(trace('ash.workhorse',where=asNamespace('ashr'),print=FALSE,exit=quote({
  if(identical(outputlevel,'loglik')) .mode_reference_capture$trials<-c(.mode_reference_capture$trials,list(c(mode,returnValue()$loglik)))
})))
profiles<-list(normal_auto=list(),normal_scales=list(scales=c(.2,.7,2)),uniform=list(family='uniform'),
  halfuniform=list(family='halfuniform'),halfnormal=list(family='halfnormal'),positive=list(family='+uniform',x=abs(x)+.3),
  negative=list(family='-uniform',x=-abs(x)-.3),student=list(family='uniform',df=4),alpha=list(alpha=.5),
  em=list(optimizer='mixEM'),weighted=list(optimizer='w_mixEM',weights=as.numeric(seq_along(x)%%7)),
  variational=list(optimizer='mixVBEM',penalty='unit'),compiled=list(optimizer='cxxMixSquarem'),
  vector_lower=list(mode=c(-8,8)),vector_upper=list(mode=c(8,-8)),vector_inner=list(mode=c(.2,.4)),
  vector_missing_first=list(mode=c(NA_real_,8)),text_pair=list(mode=c('ignored','pair')),
  fdr=list(point_mass=TRUE,nullweight=10),symmetric=list(x=c(-1,0,1),s=rep(.5,3)),
  shifted=list(x=x+20),constant=list(x=rep(1,length(x))),missing=list(x=replace(x,2,NA_real_)),
  infinite=list(x=replace(x,2,Inf)),negative_se=list(s=replace(s,2,-1)),
  custom_logistic=list(family='uniform',lik=list(name='custom_logistic',lcdfFUN=function(x)-log1p(exp(-x)),lpdfFUN=function(x)-x-2*log1p(exp(-x)))))
cases<-list()
for(name in names(profiles)){
 p<-profiles[[name]];value<-function(n,d)if(is.null(p[[n]]))d else p[[n]]
 y<-value('x',x);se<-value('s',s);family<-value('family','normal');mode<-value('mode','estimate');optimizer<-value('optimizer','mixSQP')
 scales<-value('scales',if(name=='normal_auto')NULL else c(.2,.7,2));pointmass<-value('point_mass',FALSE)
 parameters<-list(mixcompdist=family,optmethod=optimizer,mode=mode,pointmass=pointmass,prior=value('penalty',if(pointmass)'nullbiased'else'uniform'),
   mixsd=scales,alpha=value('alpha',0),outputlevel=3)
 if(!is.null(p$df))parameters$df<-p$df
 if(!is.null(p$weights))parameters$weights<-p$weights
 if(!is.null(p$lik))parameters$lik<-p$lik
 warnings<-character();.mode_reference_capture$trials<-list();.mode_reference_capture$mode<-NULL;.mode_reference_capture$bounds<-NULL
 fit<-tryCatch(withCallingHandlers(do.call(ashr::ash,c(list(betahat=y,sebetahat=se),parameters)),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 cases[[length(cases)+1L]]<-list(name=name,estimates=bits(y),standard_errors=bits(se),family=family,point_mass=pointmass,
   nullweight=if(pointmass)bits(10)[[1L]]else NULL,penalty=if(optimizer=='mixVBEM')'unit'else NULL,
   scales=bits(scales),mode=bits(0)[[1L]],mode_input=.formula_descriptor(mode),range=bits(c(-Inf,Inf)),
   alpha=bits(value('alpha',0))[[1L]],prune_threshold=bits(1e-10)[[1L]],observation_weights=bits(p$weights),
   fixed_prior=FALSE,supplied_prior=NULL,df=if(is.null(p$df))NULL else bits(p$df)[[1L]],thresholds=list(),optimizer=optimizer,control=list(),
   custom_likelihood=if(is.null(p$lik))NULL else custom,
   expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,
   native_expected_error=if(name=='negative_se')'standard errors must be nonnegative'else NULL,
   expected_mode=if(is.null(.mode_reference_capture$mode))NULL else expected(.mode_reference_capture$mode)[[1L]],
   search_bounds=bits(.mode_reference_capture$bounds),trials=lapply(.mode_reference_capture$trials,expected),
   posterior_mean=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorMean),
   posterior_sd=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorSD),
   svalues=if(inherits(fit,'error'))NULL else expected(fit$result$svalue),
   log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$loglik)[[1L]],warnings=as.list(warnings))
 # NULL storage has meaning: absence differs from an explicitly empty scale grid.
 if(is.null(scales))cases[[length(cases)]]['scales']<-list(NULL)
 if(is.null(p$weights))cases[[length(cases)]]['observation_weights']<-list(NULL)
}
suppressMessages(untrace('ash.workhorse',where=asNamespace('ashr')))
suppressMessages(untrace('ash.estmode',where=asNamespace('ashr')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr'))),cases=cases),
 args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
