#!/usr/bin/env Rscript
# Original ashr only. Capture primitive optimizer inputs, then invoke the
# original fit; candidate code and candidate outputs never enter this process.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
suppressPackageStartupMessages(library(ashr))
bits<-function(v)lapply(as.numeric(v),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
expected<-function(v)lapply(as.numeric(v),function(v)if(is.na(v)&&!is.nan(v))'NA'else sprintf('%.17g',v))
components<-function(g)list(component_class=class(g)[[1L]],weights=bits(g$pi),means=if(is.null(g$mean))NULL else bits(g$mean),sds=if(is.null(g$sd))NULL else bits(g$sd),lower=if(is.null(g$a))NULL else bits(g$a),upper=if(is.null(g$b))NULL else bits(g$b))
log_components<-function(g,data)get(paste0('log_comp_dens_conv.',class(g)[[1L]]),envir=asNamespace('ashr'))(g,data)
.ashr_reference_capture<-new.env()
suppressMessages(trace('estimate_mixprop',where=asNamespace('ashr'),print=FALSE,tracer=quote({
 c<-get('.ashr_reference_capture',envir=.GlobalEnv);c$data<-data;c$g<-g;c$prior<-prior
})))
on.exit(suppressMessages(untrace('estimate_mixprop',where=asNamespace('ashr'))))
i<-seq_len(48);x<-1.4*sin(i*.73)+.6*cos(i*.11);s<-.2+(i%%11)/13
profiles<-list(default=list(),method_one=list(control=list(method=1)),method_two=list(control=list(method=2)),
 budget_one=list(control=list(maxiter=1)),budget_two=list(control=list(maxiter=2)),budget_three=list(control=list(maxiter=3)),
 budget_four=list(control=list(maxiter=4)),budget_ten=list(control=list(maxiter=10)),tight=list(control=list(tol=1e-10)),
 bounded_steps=list(control=list(step.min0=.5,step.max0=2,mstep=2)),unmonitored=list(control=list(objfn.inc=Inf)),
 monotone=list(control=list(objfn.inc=0)),no_square=list(control=list(square=FALSE)),
 asymmetric=list(family='halfuniform'),half_normal=list(family='halfnormal'),
 uniform=list(family='uniform'),positive=list(family='+uniform',x=abs(x)),
 alpha=list(alpha=.5),null_optimal=list(x=rep(0,48)),
 early_first=list(control=list(tol=1)),early_second=list(control=list(tol=.1)),
 unit_generated=list(prior='unit'),zero_weights=list(weights=rep(0,48)),
 supplied=list(g=ashr::normalmix(c(.7,.1,.2),rep(0,3),c(0,.4,2))),
 supplied_zeros=list(g=ashr::normalmix(c(0,0,1),rep(0,3),c(0,.4,2))))
methods<-c('mixEM','w_mixEM','mixVBEM','cxxMixSquarem')
cases<-list();direct<-list()
for(name in names(profiles))for(method in methods){
 if(name=='unit_generated'&&method!='mixVBEM')next
 if(name=='zero_weights'&&method!='w_mixEM')next
 p<-profiles[[name]];getp<-function(key,default)if(is.null(p[[key]]))default else p[[key]]
 y<-getp('x',x);se<-s;family<-getp('family','normal');control<-getp('control',list())
 parameters<-list(mixcompdist=family,optmethod=method,pointmass=TRUE,prior=getp('prior','nullbiased'),mixsd=c(.2,.7,2),control=control,alpha=getp('alpha',0))
 if(method=='w_mixEM')parameters$weights<-getp('weights',as.numeric(i%%7))
 if(!is.null(p$g))parameters$g<-p$g
 warnings<-character()
 result<-tryCatch({
  fit<-withCallingHandlers(do.call(ashr::ash,c(list(betahat=y,sebetahat=se,outputlevel=3),parameters)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  ll<-t(log_components(.ashr_reference_capture$g,.ashr_reference_capture$data));ll<-ll[!ashr:::get_exclusions(.ashr_reference_capture$data),,drop=FALSE]
  pi<-fit$fitted_g$pi
  native_pi<-pi;fresh_original_objective<-NULL;corrected_weights<-NULL
  if(method=='cxxMixSquarem'&&!is.null(fit$fit_details$optreturn$B)) {
    mat<-fit$fit_details$matrix_lik;pr<-.ashr_reference_capture$prior
    fresh_original_objective<-expected(ashr:::negpenloglik(pi,mat,pr))[[1L]]
    p1<-ashr:::fixpoint(pi,mat,pr);p2<-ashr:::fixpoint(p1,mat,pr)
    tol<-ashr:::set_control_squarem(control,nrow(mat))$tol
    if(isTRUE(fit$fit_details$optreturn$converged)&&sqrt(sum((p1-pi)^2))>=tol&&sqrt(sum((p2-p1)^2))<tol) {
      native_pi<-p2;corrected_weights<-expected(p2)
    }
  }
  # likelihood input bits, initial state and prior are captured from original
  # preparation. This panel separates optimizer differences from components.
  direct[[length(direct)+1L]]<-list(name=paste(name,method,sep='_'),optimizer=method,observations=nrow(ll),components=ncol(ll),
    log_likelihoods=bits(t(ll)),observation_weights=if(is.null(parameters$weights))NULL else bits(parameters$weights),dirichlet=bits(.ashr_reference_capture$prior),initial_weights=bits(.ashr_reference_capture$g$pi),control=control,
    expected_weights=expected(pi),corrected_weights=corrected_weights,fresh_original_objective=fresh_original_objective,converged=fit$fit_details$optreturn$converged,iterations=fit$fit_details$optreturn$niter,optimizer_objective=if(is.null(fit$fit_details$optreturn$B))NULL else expected(fit$fit_details$optreturn$B)[[1L]],native_objective_reference=if(is.null(fit$fit_details$optreturn$B))NULL else if(method=='cxxMixSquarem')expected(ashr:::negpenloglik(native_pi,fit$fit_details$matrix_lik,.ashr_reference_capture$prior))[[1L]]else expected(fit$fit_details$optreturn$B)[[1L]],posterior_concentrations=if(is.null(fit$fit_details$optreturn$post))NULL else expected(fit$fit_details$optreturn$post))
  list(fitted_prior=components(fit$fitted_g),posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD),svalues=expected(fit$result$svalue),local_false_sign_rate=expected(fit$result$lfsr),log_likelihood=expected(fit$loglik)[[1L]],converged=fit$fit_details$optreturn$converged)
 },error=function(e)list(expected_error=conditionMessage(e)))
 cases[[length(cases)+1L]]<-c(list(name=paste(name,method,sep='_'),estimates=bits(y),standard_errors=bits(se),optimizer=method,control=control,
   family=family,point_mass=TRUE,nullweight=if(name=='unit_generated')NULL else bits(10)[[1L]],penalty=if(name=='unit_generated')'unit'else NULL,scales=bits(c(.2,.7,2)),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(getp('alpha',0))[[1L]],prune_threshold=bits(1e-10)[[1L]],observation_weights=if(is.null(parameters$weights))NULL else bits(parameters$weights),fixed_prior=FALSE,supplied_prior=if(is.null(p$g))NULL else components(p$g),thresholds=list()),result,list(warnings=as.list(warnings)))
}
# Explicit variational priors, including source unit=1/K and positive parameters
# below one. No convex optimizer is permitted to use these priors.
for(prior in list('unit',c(.2,.4,.6))) {
 g<-ashr::normalmix(c(.6,.2,.2),rep(0,3),c(0,.4,2));w<-character()
 fit<-withCallingHandlers(ashr::ash(x,s,g=g,prior=prior,optmethod='mixVBEM',outputlevel=3),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 ll<-t(log_components(g,.ashr_reference_capture$data))
 name<-if(is.character(prior))'unit_variational'else'fractional_variational'
 direct[[length(direct)+1L]]<-list(name=name,optimizer='mixVBEM',observations=nrow(ll),components=ncol(ll),log_likelihoods=bits(t(ll)),observation_weights=NULL,dirichlet=bits(.ashr_reference_capture$prior),initial_weights=bits(g$pi),control=list(),expected_weights=expected(fit$fitted_g$pi),converged=fit$fit_details$optreturn$converged,iterations=fit$fit_details$optreturn$niter,optimizer_objective=if(is.null(fit$fit_details$optreturn$B))NULL else expected(fit$fit_details$optreturn$B)[[1L]],native_objective_reference=expected(fit$fit_details$optreturn$B)[[1L]],posterior_concentrations=expected(fit$fit_details$optreturn$post))
 cases[[length(cases)+1L]]<-list(name=name,estimates=bits(x),standard_errors=bits(s),optimizer='mixVBEM',control=list(),family='normal',point_mass=TRUE,nullweight=NULL,penalty=if(is.character(prior))'unit'else bits(prior),scales=bits(c(.2,.7,2)),mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),alpha=bits(0)[[1L]],prune_threshold=bits(1e-10)[[1L]],observation_weights=NULL,fixed_prior=FALSE,supplied_prior=components(g),thresholds=list(),fitted_prior=components(fit$fitted_g),posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD),svalues=expected(fit$result$svalue),local_false_sign_rate=expected(fit$result$lfsr),log_likelihood=expected(fit$loglik)[[1L]],converged=fit$fit_details$optreturn$converged,warnings=as.list(w))
}
# Original compiled second-probe convergence returns the initial state even
# when a categorical update reaches a different stationary point exactly.
L<-matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);a<-c(2,3);start<-c(.5,.5)
f<-suppressWarnings(ashr:::cxxMixSquarem(L,a,start,ashr:::set_control_squarem(list(),4)))
truth<-c(3+a[[1L]]-1,1+a[[2L]]-1);truth<-truth/sum(truth)
stopifnot(identical(as.numeric(f$pihat),start),max(abs(truth-start))>.05)
direct[[length(direct)+1L]]<-list(name='categorical_compiled_stop',optimizer='cxxMixSquarem',observations=4L,components=2L,log_likelihoods=bits(t(log(L))),observation_weights=NULL,dirichlet=bits(a),initial_weights=bits(start),control=list(),expected_weights=expected(f$pihat),fresh_original_objective=expected(ashr:::negpenloglik(f$pihat,L,a))[[1L]],corrected_weights=expected(truth),converged=f$converged,iterations=f$niter,posterior_concentrations=NULL,optimizer_objective=expected(f$B)[[1L]],native_objective_reference=expected(ashr:::negpenloglik(truth,L,a))[[1L]],correction='commit the second stationary update instead of returning unchanged starting weights')
suppressMessages(untrace('estimate_mixprop',where=asNamespace('ashr')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),SQUAREM=as.character(packageVersion('SQUAREM'))),cases=cases,em_fit_cases=direct),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
