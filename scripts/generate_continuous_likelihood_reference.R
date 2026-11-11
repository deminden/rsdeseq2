#!/usr/bin/env Rscript
# Independent original-ashr kernels and full continuous-likelihood workflows.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
models<-c('normalmix','normalmix_rows','normalmix_unnormalized','normalmix_zero','logF','logF_rows','logF_small','logF_large','logF_first_infinite','logF_second_infinite','logF_both_infinite')
make<-function(name,n){
 if(startsWith(name,'normalmix')){
   w<-if(name=='normalmix_unnormalized')c(1,2)else c(.7,.3)
   s<-if(name=='normalmix_rows')cbind(seq(.3,1.2,length.out=n),seq(1.3,3,length.out=n))else if(name=='normalmix_zero')c(0,1)else c(.5,2)
   return(list(lik=ashr::lik_normalmix(w,s),model=list(type='normalmix',weights=bits(w),scales=.formula_descriptor(s))))
 }
 d1<-switch(name,logF_rows=seq(2,9,length.out=n),logF_small=.5,logF_large=1000,logF_first_infinite=Inf,logF_both_infinite=Inf,5)
 d2<-switch(name,logF_rows=seq(10,3,length.out=n),logF_small=.75,logF_large=1500,logF_second_infinite=Inf,logF_both_infinite=Inf,8)
 list(lik=ashr::lik_logF(d1,d2),model=list(type='logF',df1=bits(d1),df2=bits(d2)))
}
q<-c(-50,-4,-.1,0,.1,4,50,9.5,38);a<-c(-Inf,-5,-.5,0,.2,2,10,9.5,9.53);b<-c(-4,-3,.7,0,2,5,Inf,38,40)
kernels<-lapply(models,function(name){
 p<-make(name,length(q));lik<-ashr:::add_etruncFUN(p$lik);warnings<-character()
 value<-tryCatch(withCallingHandlers(list(log_cdf=expected(lik$lcdfFUN(q)),log_density=expected(lik$lpdfFUN(q)),
   first=expected(lik$etruncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),second=expected(lik$e2truncFUN(matrix(a,ncol=1),matrix(b,ncol=1)))),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(expected_error=conditionMessage(e)))
 c(list(name=name,builtin_likelihood=p$model,arguments=bits(q),lower=bits(a),upper=bits(b)),value,list(warnings=as.list(warnings)))
})
x<-c(-2,-.4,0,.8,2,4);s<-c(.3,.5,1,1.2,.7,.4)
settings<-c(lapply(models,function(n)list(name=paste0(n,'_fixed'),model=n,fixed=TRUE)),
 list(list(name='normalmix_fitted',model='normalmix'),list(name='normalmix_rows_fitted',model='normalmix_rows'),
 list(name='logF_fitted',model='logF'),list(name='logF_rows_fitted',model='logF_rows'),
 list(name='normalmix_missing',model='normalmix_rows',fixed=TRUE,missing=TRUE),list(name='normalmix_fitted_missing',model='normalmix_rows',missing=TRUE),list(name='logF_missing',model='logF_rows',fixed=TRUE,missing=TRUE),
 list(name='normalmix_bad_prior',model='normalmix',family='normal'),list(name='logF_bad_prior',model='logF',family='normal'),
 list(name='logF_tail_cancellation',model='logF_first_infinite',fixed=TRUE,tail=TRUE)))
cases<-lapply(settings,function(p){
 datax<-x;if(isTRUE(p$missing))datax[[2L]]<-NA_real_
 se<-s;if(isTRUE(p$tail)){datax[1:3]<-c(3.2,3.3,3.4);se[1:3]<-.14}
 model<-make(p$model,length(x));fixed<-isTRUE(p$fixed);family<-if(is.null(p$family))'uniform'else p$family
 atom<-p$model!='normalmix_zero';g<-if(atom)ashr::unimix(c(.2,.3,.5),c(0,-.4,-3),c(0,.4,3))else ashr::unimix(c(.4,.6),c(-.4,-3),c(.4,3))
 warnings<-character();parameters<-list(mixcompdist=family,pointmass=atom,prior='uniform',mixsd=c(.4,3),lik=model$lik,outputlevel=3)
 if(fixed){parameters$g<-g;parameters$fixg<-TRUE}
 fit<-tryCatch(withCallingHandlers(do.call(ashr::ash,c(list(betahat=datax,sebetahat=se),parameters)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 cutoffs<-c(-1,0,1)
 cdf<-if(inherits(fit,'error'))NULL else vapply(cutoffs,function(q)ashr::cdf_post(fit$fitted_g,q,ashr::set_data(datax,se,lik=ashr:::add_etruncFUN(model$lik))),numeric(length(x)))
 supplied<-if(!fixed)NULL else list(component_class='unimix',weights=bits(g$pi),means=NULL,sds=NULL,lower=bits(g$a),upper=bits(g$b))
 list(name=p$name,builtin_likelihood=model$model,estimates=bits(datax),standard_errors=bits(se),family=family,point_mass=atom,nullweight=NULL,
   scales=bits(c(.4,3)),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(0)[[1L]],prune_threshold=bits(1e-10)[[1L]],
   observation_weights=NULL,fixed_prior=fixed,supplied_prior=supplied,df=NULL,thresholds=bits(cutoffs),optimizer='mixSQP',control=list(),
   expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,
   posterior_mean=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorMean),posterior_sd=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorSD),
   svalues=if(inherits(fit,'error'))NULL else expected(fit$result$svalue),cdf=if(is.null(cdf))NULL else expected(t(cdf)),log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$loglik)[[1L]],warnings=as.list(warnings))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr'))),kernel_cases=kernels,cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
