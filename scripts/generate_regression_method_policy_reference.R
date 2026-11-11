#!/usr/bin/env Rscript
# Original-package calculations only; reuse frozen independent inputs, never candidate outputs.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
decode<-function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(sub('bits:','',h),seq(1,15,2),seq(2,16,2)),16L)),'double',n=1,size=8,endian='big'),numeric(1))
mat<-function(m)if(is.null(m))NULL else matrix(decode(m$values),nrow=m$rows,ncol=m$columns)
descriptor<-function(m)list(rows=nrow(m),columns=ncol(m),values=bits(m))
read_cases<-function(name)jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data',name),simplifyVector=FALSE)$cases
regression<-read_cases('regression_posterior_r462.json');beta<-read_cases('beta_binomial_posterior_r462.json')
normal<-function(y,x,beta,param,offset)stats::dnorm(y,mean=x%*%beta+offset,sd=param,log=TRUE)
.method_policy_calls<-new.env()
suppressMessages(trace('optimBetabin',where=asNamespace('apeglm'),print=FALSE,
  tracer=quote({.method_policy_calls$beta<-.method_policy_calls$beta+1L})))
methods<-c('general','nbinomR','betabinR','betabinCR','betabinC','betabinC*')
profiles<-unlist(lapply(methods,function(method)list(
 list(name=paste0(method,'_empty_bounds'),method=method,optimizer='L-BFGS-B',bounds=numeric()),
 list(name=paste0(method,'_unknown_optimizer'),method=method,optimizer='not_an_optimizer',bounds=numeric()))),recursive=FALSE)
profiles<-c(profiles,list(
 list(name='general_nonfinite_bounds',method='general',optimizer='L-BFGS-B',bounds=c(NA_real_,NaN)),
 list(name='nbinomR_short_bounds',method='nbinomR',optimizer='L-BFGS-B',bounds=-1),
 list(name='betabinR_extra_bounds',method='betabinR',optimizer='L-BFGS-B',bounds=c(-10,10,0)),
 list(name='betabinR_ignored_bounds',method='betabinR',optimizer='BFGS',bounds=c(Inf,-Inf,NaN)),
 list(name='betabinCR_unknown_tight_bounds',method='betabinCR',optimizer='not_an_optimizer',bounds=c(-.01,.01)),
 list(name='betabinCR_refinement_unknown',method='betabinCR',optimizer='not_an_optimizer',bounds=numeric(),refinement=TRUE),
 list(name='betabinCR_refinement_empty_bounds',method='betabinCR',optimizer='L-BFGS-B',bounds=numeric(),refinement=TRUE),
 list(name='general_prefix_ignores_bounds',method='general',optimizer='L-B',bounds=c(-.001,.001)),
 list(name='nbinomR_prefix_ignores_bounds',method='nbinomR',optimizer='L-B',bounds=c(-.001,.001))))
cases<-lapply(profiles,function(p){
 c<-if(startsWith(p$method,'betabin'))beta[[1L]]else regression[[if(p$method=='general')1L else 7L]]
 c$name<-p$name;c$method<-p$method;c$optimizer<-p$optimizer;c$bounds<-bits(p$bounds)
 c$expected_error<-c$native_expected_error<-NULL;c$rng_state<-c$expected_random_state<-NULL
 y<-mat(c$observations);x<-mat(c$design);param<-mat(c$parameters)
 if(startsWith(p$method,'betabin')&&!isTRUE(p$refinement)){
  # Balanced interior counts, selected solely using original-package prefits.
  y<-matrix(rep(2:5,2),nrow=1);param<-cbind(8,matrix(10,nrow=1,ncol=8))
  c$observations<-descriptor(y);c$parameters<-descriptor(param)
 }
 pc<-list(no.shrink=as.numeric(unlist(c$nuisance))+1L,prior.mean=decode(c$prior_mean),
  prior.scale=decode(c$scale),prior.df=decode(c$prior_df),prior.no.shrink.mean=decode(c$nuisance_mean),prior.no.shrink.scale=decode(c$nuisance_scale))
 state<-function(v)lapply(v,function(x)paste(format(writeBin(as.integer(x),raw(),size=4L,endian='big')),collapse=''))
 if(p$method=='betabinC*'){RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42L);c$rng_state<-state(.Random.seed)}
 warnings<-character();.method_policy_calls$beta<-0L
 fit<-tryCatch(withCallingHandlers(apeglm::apeglm(Y=y,x=x,log.lik=normal,param=param,coef=c$coefficient+1L,
  weights=mat(c$weights),offset=if(p$method=='general')matrix(0,nrow(y),ncol(y))else mat(c$offsets),
  prior.control=pc,method=p$method,optim.method=p$optimizer,bounds=p$bounds,threshold=decode(c$threshold)),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 if(p$method=='betabinC*')c$expected_random_state<-state(.Random.seed)
 c$original_beta_optimizer_calls<-.method_policy_calls$beta
 if(p$method%in%c('betabinC','betabinC*','betabinCR')){
  if(isTRUE(p$refinement))stopifnot(.method_policy_calls$beta>0L)else stopifnot(.method_policy_calls$beta==0L)
 }
 if(inherits(fit,'error')){
  c$expected_error<-conditionMessage(fit);c$native_expected_error<-'unsupported posterior optimizer'
 }else{
  c$mode<-expected(fit$map);c$sd<-expected(fit$sd);c$fsr<-expected(fit$fsr);c$sign_svalues<-expected(fit$svalue)
  c$fsos<-expected(fit$thresh);c$threshold_svalues<-expected(apeglm::svalue(fit$thresh))
  c$fitted_scale<-expected(fit$prior.control$prior.scale)[[1L]];c$prior_variance<-NULL
 }
 c$warnings<-as.list(warnings);c
})
suppressMessages(untrace('optimBetabin',where=asNamespace('apeglm')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm'))),cases=cases),
 args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
