#!/usr/bin/env Rscript
# Independent original-ashr count kernels and complete count-likelihood workflows.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
models<-c('pois_identity','pois_identity_rows','pois_log','pois_log_rows','pois_log_zero','pois_large','binom_identity','binom_identity_rows','binom_logit','binom_logit_rows','binom_logit_zero','binom_large')
make<-function(name,n){
 y<-if(grepl('rows',name))rep(c(0,1,3,7,2),length.out=n)else if(grepl('zero',name))0 else if(grepl('large',name))800 else 2
 if(startsWith(name,'pois')){
   scale<-if(grepl('rows',name))seq(.5,2,length.out=n)else 1.3
   link<-if(grepl('log',name))'log'else 'identity'
   return(list(lik=ashr::lik_pois(y,scale,link),model=list(type='pois',y=bits(y),scale=bits(scale),link=link)))
 }
 trials<-if(grepl('large',name))1000 else if(grepl('rows',name))rep(c(8,1,5,9,2),length.out=n)else 8
 link<-if(grepl('logit',name))'logit'else 'identity'
 list(lik=ashr::lik_binom(y,trials,link),model=list(type='binom',y=bits(y),trials=bits(trials),link=link))
}
q<-c(-20,-2,-1,-.3,0,.3,1,2,20);a<-c(-Inf,-4,-1,-.5,0,0,.5,2,10);b<-c(-3,-2,0,.5,0,1,1,4,Inf)
kernels<-lapply(models,function(name){
 p<-make(name,length(q));lik<-ashr:::add_etruncFUN(p$lik);warnings<-character()
 value<-tryCatch(withCallingHandlers(list(log_cdf=expected(lik$lcdfFUN(q)),log_density=expected(lik$lpdfFUN(q)),
   first=expected(lik$etruncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),second=expected(lik$e2truncFUN(matrix(a,ncol=1),matrix(b,ncol=1)))),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(expected_error=conditionMessage(e)))
 c(list(name=name,builtin_likelihood=p$model,arguments=bits(q),lower=bits(a),upper=bits(b)),value,list(warnings=as.list(warnings)))
})
x<-rep(0,6);s<-rep(1,6)
settings<-c(lapply(models,function(n)list(name=paste0(n,'_fixed'),model=n,fixed=TRUE)),
 lapply(c('pois_identity_rows','pois_log_rows','binom_identity_rows','binom_logit_rows'),function(n)list(name=paste0(n,'_fitted'),model=n)),
 list(list(name='pois_missing',model='pois_identity_rows',fixed=TRUE,missing=TRUE),list(name='binom_missing',model='binom_identity_rows',fixed=TRUE,missing=TRUE),
 list(name='pois_global_orientation',model='pois_identity',fixed=TRUE,global=TRUE),list(name='binom_global_orientation',model='binom_identity',fixed=TRUE,global=TRUE),list(name='pois_fdr',model='pois_identity_rows',fdr=TRUE),list(name='binom_fdr',model='binom_logit_rows',fdr=TRUE),list(name='pois_bad_prior',model='pois_identity',family='normal'),list(name='binom_bad_prior',model='binom_identity',family='normal')))
cases<-lapply(settings,function(p){
 datax<-x;if(isTRUE(p$missing))datax[[2L]]<-NA_real_
 if(isTRUE(p$global))datax<-c(0,.7,.3,.9,.2,.8)
 se<-s
 model<-make(p$model,length(x));fixed<-isTRUE(p$fixed);family<-if(is.null(p$family))'uniform'else p$family
 bounded<-identical(model$model$link,'identity');atom<-isTRUE(p$fdr)
 lower<-if(bounded) c(0,0)else c(-1,-4);upper<-if(bounded)if(model$model$type=='binom')c(.3,1)else c(1,5)else c(1,4)
 g<-if(isTRUE(p$global))ashr::unimix(1,0,1)else ashr::unimix(c(.4,.6),lower,upper)
 warnings<-character();parameters<-list(mixcompdist=family,pointmass=atom,prior=if(isTRUE(p$fdr))'nullbiased'else'uniform',mixsd=c(.4,3),lik=model$lik,outputlevel=3)
 if(fixed){parameters$g<-g;parameters$fixg<-TRUE}else parameters$mode<-0
 fit<-tryCatch(withCallingHandlers(do.call(ashr::ash,c(list(betahat=datax,sebetahat=se),parameters)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 cutoffs<-c(-1,0,1)
 cdf<-if(inherits(fit,'error'))NULL else vapply(cutoffs,function(q)ashr::cdf_post(fit$fitted_g,q,ashr::set_data(datax,se,lik=ashr:::add_etruncFUN(model$lik))),numeric(length(x)))
 supplied<-if(!fixed)NULL else list(component_class='unimix',weights=bits(g$pi),means=NULL,sds=NULL,lower=bits(g$a),upper=bits(g$b))
 list(name=p$name,builtin_likelihood=model$model,estimates=bits(datax),standard_errors=bits(se),family=family,point_mass=atom,nullweight=if(isTRUE(p$fdr))bits(1)[[1L]]else NULL,
   scales=bits(c(.4,3)),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(0)[[1L]],prune_threshold=bits(1e-10)[[1L]],
   observation_weights=NULL,fixed_prior=fixed,supplied_prior=supplied,df=NULL,thresholds=bits(cutoffs),optimizer='mixSQP',control=list(),
   expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,
   posterior_mean=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorMean),posterior_sd=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorSD),
   svalues=if(inherits(fit,'error'))NULL else expected(fit$result$svalue),cdf=if(is.null(cdf))NULL else expected(t(cdf)),log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$loglik)[[1L]],warnings=as.list(warnings))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr'))),kernel_cases=kernels,cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
