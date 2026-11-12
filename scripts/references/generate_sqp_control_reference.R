#!/usr/bin/env Rscript
# Original ashr/mixsqp only; capture full fits and primitive likelihood inputs.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
decode<-function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(sub('bits:','',h),seq(1,15,2),seq(2,16,2)),16L)),'double',n=1L,size=8,endian='big'),numeric(1))
original<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/ashr_optimizer_r462.json'),simplifyVector=FALSE)$cases[[1L]]
x<-decode(original$estimates);s<-decode(original$standard_errors)
.sqp_reference_capture<-new.env()
suppressMessages(trace('estimate_mixprop',where=asNamespace('ashr'),print=FALSE,tracer=quote({c<-get('.sqp_reference_capture',.GlobalEnv);c$data<-data;c$g<-g;c$prior<-prior})))
suppressMessages(trace('tsvd',where=asNamespace('mixsqp'),print=FALSE,exit=quote({c<-get('.sqp_reference_capture',.GlobalEnv);v<-returnValue();c$rank<-if(is.null(v))NULL else ncol(v$U)})))
profiles<-list(default=list(),active_tight=list(convtol.activeset=1e-12),active_loose=list(convtol.activeset=1e-4),
 direction_tight=list(zero.threshold.searchdir=1e-14),direction_loose=list(zero.threshold.searchdir=1e-3),
 armijo=list(suffdecr.linesearch=.2),reduce_half=list(stepsizereduce=.5),reduce_slow=list(stepsizereduce=.9),
 step_large=list(minstepsize=.1),step_full=list(minstepsize=1),ridge_small=list(identity.contrib.increase=2),ridge_large=list(identity.contrib.increase=100),
 normalize_off=list(normalize.rows=FALSE),svd_off=list(tol.svd=0),svd_tight=list(tol.svd=1e-8),svd_medium=list(tol.svd=.001),svd_coarse=list(tol.svd=.05),
 active_auto=list(maxiter.activeset=NULL),combined=list(convtol.activeset=1e-12,zero.threshold.searchdir=1e-12,suffdecr.linesearch=.1,stepsizereduce=.5,minstepsize=1e-6),
 short_sqp=list(maxiter.sqp=1,numiter.em=0),short_active=list(maxiter.activeset=1,numiter.em=0),
 invalid_active=list(convtol.activeset=-1),invalid_direction=list(zero.threshold.searchdir=-1),invalid_reduction=list(stepsizereduce=1),
 invalid_step=list(minstepsize=0),invalid_ridge=list(identity.contrib.increase=0),invalid_svd=list(tol.svd=-1))
settings<-lapply(names(profiles),function(n)list(name=n,control=profiles[[n]],family='normal',pointmass=FALSE))
settings<-c(settings,list(list(name='asymmetric_line',control=list(suffdecr.linesearch=.2,stepsizereduce=.5),family='halfuniform',pointmass=TRUE),
 list(name='prior_svd',control=list(tol.svd=.05),family='normal',pointmass=TRUE),list(name='prior_auto',control=list(maxiter.activeset=NULL),family='normal',pointmass=TRUE),
 list(name='wide_svd',control=list(tol.svd=.01),family='normal',pointmass=FALSE,scales=exp(seq(log(.05),log(6),length.out=20))),
 list(name='wide_prior_svd',control=list(tol.svd=.01),family='normal',pointmass=TRUE,scales=exp(seq(log(.05),log(6),length.out=20))),
 list(name='wide_prior_truncated',control=list(tol.svd=.1),family='normal',pointmass=TRUE,scales=exp(seq(log(.05),log(6),length.out=20)))))
cases<-list();direct<-list()
exact_svd_fit<-function(expression){
 ns<-asNamespace('mixsqp');old<-get('tsvd',ns)
 on.exit({unlockBinding('tsvd',ns);assign('tsvd',old,ns);lockBinding('tsvd',ns)},add=TRUE)
 # Independent LAPACK spectral oracle. Keep original rank-discovery/fallback
 # bounds; replace randomized partial vectors only in this labeled reference.
 exact<-function(X,tol){
   v<-base::svd(X);r<-min(dim(X));probe<-2L
   repeat{if(probe>floor(r/2))return(NULL);if(v$d[[probe]]<tol)break;probe<-min(2L*probe,r)}
   keep<-seq_len(max(2L,sum(v$d>tol)));d<-sqrt(v$d[keep])
   list(U=sweep(v$u[,keep,drop=FALSE],2L,d,'*'),V=sweep(v$v[,keep,drop=FALSE],2L,d,'*'))
 }
 unlockBinding('tsvd',ns);assign('tsvd',exact,ns);lockBinding('tsvd',ns)
 force(expression)
}
for(p in settings){
 set.seed(6001);warnings<-character();.sqp_reference_capture$rank<-NULL
 scales<-if(is.null(p$scales))c(.1,.2,.5,1,2,4)else p$scales
 fit<-tryCatch(withCallingHandlers(ashr::ash(x,s,mixcompdist=p$family,pointmass=p$pointmass,prior=if(p$pointmass)'nullbiased'else'uniform',
   mixsd=scales,optmethod='mixSQP',control=p$control,outputlevel=3),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 c<-.sqp_reference_capture;ll<-t(get(paste0('log_comp_dens_conv.',class(c$g)[[1L]]),asNamespace('ashr'))(c$g,c$data))
 direct[[length(direct)+1L]]<-list(name=p$name,optimizer='mixSQP',observations=nrow(ll),components=ncol(ll),log_likelihoods=bits(t(ll)),
   observation_weights=NULL,dirichlet=bits(c$prior),initial_weights=bits(c$g$pi),control=p$control,
   expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,expected_weights=if(inherits(fit,'error'))NULL else expected(fit$fitted_g$pi),
   expected_log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$loglik)[[1L]],
   expected_penalized_log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$fit_details$penloglik)[[1L]],
   converged=if(inherits(fit,'error'))NULL else fit$fit_details$optreturn$converged,source_svd_rank=c$rank)
 cases[[length(cases)+1L]]<-list(name=p$name,estimates=bits(x),standard_errors=bits(s),family=p$family,point_mass=p$pointmass,
   nullweight=if(p$pointmass)bits(10)[[1L]]else NULL,scales=bits(scales),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(0)[[1L]],
   prune_threshold=bits(1e-10)[[1L]],observation_weights=NULL,fixed_prior=FALSE,supplied_prior=NULL,df=NULL,thresholds=list(),optimizer='mixSQP',control=p$control,
   expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,posterior_mean=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorMean),
   posterior_sd=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorSD),svalues=if(inherits(fit,'error'))NULL else expected(fit$result$svalue),
   log_likelihood=if(inherits(fit,'error'))NULL else expected(fit$loglik)[[1L]],warnings=as.list(warnings))
 if(!inherits(fit,'error')&&!is.null(c$rank)&&c$rank<ncol(ll)){
   deterministic<-exact_svd_fit(suppressWarnings(ashr::ash(x,s,mixcompdist=p$family,pointmass=p$pointmass,
     prior=if(p$pointmass)'nullbiased'else'uniform',mixsd=scales,optmethod='mixSQP',control=p$control,outputlevel=3)))
   corrected<-list(provenance='Original mixsqp with an independent deterministic LAPACK SVD',
     expected_weights=expected(deterministic$fitted_g$pi),expected_log_likelihood=expected(deterministic$loglik)[[1L]],
     expected_penalized_log_likelihood=expected(deterministic$fit_details$penloglik)[[1L]],converged=deterministic$fit_details$optreturn$converged,
     posterior_mean=expected(deterministic$result$PosteriorMean),posterior_sd=expected(deterministic$result$PosteriorSD),svalues=expected(deterministic$result$svalue),log_likelihood=expected(deterministic$loglik)[[1L]])
   direct[[length(direct)]]$deterministic_reference<-corrected;cases[[length(cases)]]$deterministic_reference<-corrected
 }
}
suppressMessages(untrace('estimate_mixprop',where=asNamespace('ashr')));suppressMessages(untrace('tsvd',where=asNamespace('mixsqp')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),mixsqp=as.character(packageVersion('mixsqp')),
 irlba=as.character(packageVersion('irlba')),optimizer_family='SQP'),cases=cases,em_fit_cases=direct),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
