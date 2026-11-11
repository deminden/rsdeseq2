#!/usr/bin/env Rscript
# Original R packages and frozen original inputs only; never load candidate Rust.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
suppressPackageStartupMessages(library(ashr))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
base<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/ashr_optimizer_r462.json'),simplifyVector=FALSE)
bits<-function(v)lapply(as.numeric(v),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
expected<-function(v)lapply(as.numeric(v),function(v)if(is.na(v)&&!is.nan(v))'NA'else sprintf('%.17g',v))
decode<-function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(sub('bits:','',h),seq(1,15,2),seq(2,16,2)),16L)),'double',n=1L,size=8L,endian='big'),numeric(1))
components<-function(g)list(component_class=class(g)[[1L]],weights=bits(g$pi),means=if(is.null(g$mean))NULL else bits(g$mean),sds=if(is.null(g$sd))NULL else bits(g$sd),lower=if(is.null(g$a))NULL else bits(g$a),upper=if(is.null(g$b))NULL else bits(g$b))
profiles<-list(rre2=list(K=2),mpe2=list(K=2,method='mpe'),rre2_unsquared=list(K=2,square=FALSE),mpe2_unsquared=list(K=2,method='mpe',square=FALSE),
 rre3=list(K=3),mpe3=list(K=3,method='mpe'),mpe3_unsquared=list(K=3,method='mpe',square=FALSE),
 budget_zero=list(K=2,maxiter=0),budget_one=list(K=2,maxiter=1),budget_three=list(K=2,maxiter=3),budget_seven=list(K=2,maxiter=7),
 tight=list(K=2,method='mpe',tol=1e-10),monotone=list(K=2,method='mpe',objfn.inc=0))
cases<-list();direct<-list()
for(name in names(profiles))for(method in c('mixEM','w_mixEM','mixVBEM','cxxMixSquarem')) {
 input<-Filter(function(v)v$name==paste0('default_',method),base$cases)[[1L]]
 primitive<-Filter(function(v)v$name==paste0('default_',method),base$em_fit_cases)[[1L]]
 control<-profiles[[name]];label<-paste(name,method,sep='_');input$name<-primitive$name<-label
 input$control<-primitive$control<-control
 x<-decode(input$estimates);s<-decode(input$standard_errors)
 arguments<-list(betahat=x,sebetahat=s,mixcompdist='normal',optmethod=method,pointmass=TRUE,prior='nullbiased',mixsd=c(.2,.7,2),control=control,outputlevel=3)
 if(method=='w_mixEM')arguments$weights<-decode(input$observation_weights)
 warnings<-character()
 result<-tryCatch({
  f<-withCallingHandlers(do.call(ashr::ash,arguments),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  list(fitted_prior=components(f$fitted_g),posterior_mean=expected(f$result$PosteriorMean),posterior_sd=expected(f$result$PosteriorSD),svalues=expected(f$result$svalue),local_false_sign_rate=expected(f$result$lfsr),log_likelihood=expected(f$loglik)[[1L]],converged=f$fit_details$optreturn$converged)
 },error=function(e)list(expected_error=conditionMessage(e)))
 # Retain only original numerical outputs for this new control request.
 input<-input[c('name','estimates','standard_errors','optimizer','control','family','point_mass','nullweight','scales','mode','range','alpha','prune_threshold','observation_weights','fixed_prior','supplied_prior','thresholds')]
 cases[[length(cases)+1L]]<-c(input,result,list(warnings=as.list(warnings)))
 ll<-matrix(decode(primitive$log_likelihoods),nrow=primitive$observations,byrow=TRUE);lik<-exp(ll-apply(ll,1,max))
 opt<-tryCatch({
  parameters<-list(matrix_lik=lik,prior=decode(primitive$dirichlet),pi_init=decode(primitive$initial_weights),control=control)
  if(method=='mixVBEM')parameters$pi_init<-NULL
  if(method=='w_mixEM')parameters$weights<-decode(primitive$observation_weights)
  if(method=='cxxMixSquarem')parameters$control<-ashr:::set_control_squarem(control,nrow(lik))
  r<-suppressWarnings(do.call(get(method,asNamespace('ashr')),parameters))
  list(expected_weights=expected(r$pihat),converged=r$converged,iterations=r$niter,posterior_concentrations=if(is.null(r$post))NULL else expected(r$post))
 },error=function(e)list(expected_error=conditionMessage(e)))
 primitive<-primitive[c('name','optimizer','observations','components','log_likelihoods','observation_weights','dirichlet','initial_weights','control')]
 direct[[length(direct)+1L]]<-c(primitive,opt)
}
# A genuine independent optimizer reference for the same unregularized convex
# objective, not an assertion that original Mosek-backed mixIP was executed.
for(weighted in c(FALSE,TRUE)) {
 label<-if(weighted)'interior_weighted'else'interior_unweighted'
 method<-if(weighted)'w_mixEM'else'mixEM'
 input<-Filter(function(v)v$name==paste0('default_',method),base$cases)[[1L]]
 primitive<-Filter(function(v)v$name==paste0('default_',method),base$em_fit_cases)[[1L]]
 x<-decode(input$estimates);s<-decode(input$standard_errors)
 parameters<-list(betahat=x,sebetahat=s,mixcompdist='normal',optmethod='mixSQP',pointmass=TRUE,prior='nullbiased',mixsd=c(.2,.7,2),outputlevel=3,control=list(eps=0,convtol.sqp=1e-12,maxiter.sqp=1000,verbose=FALSE))
 if(weighted)parameters$weights<-decode(input$observation_weights)
 f<-suppressWarnings(do.call(ashr::ash,parameters))
 input<-input[c('estimates','standard_errors','family','point_mass','nullweight','scales','mode','range','alpha','prune_threshold','observation_weights','fixed_prior','supplied_prior','thresholds')]
 cases[[length(cases)+1L]]<-c(list(name=label,optimizer='mixIP',control=list(rtol=1e-6),reference_optimizer='original mixSQP with eps=0'),input,list(fitted_prior=components(f$fitted_g),posterior_mean=expected(f$result$PosteriorMean),posterior_sd=expected(f$result$PosteriorSD),svalues=expected(f$result$svalue),local_false_sign_rate=expected(f$result$lfsr),log_likelihood=expected(f$loglik)[[1L]],converged=f$fit_details$optreturn$converged,warnings=list()))
 primitive<-primitive[c('observations','components','log_likelihoods','observation_weights','dirichlet','initial_weights')]
 direct[[length(direct)+1L]]<-c(list(name=label,optimizer='mixIP',control=list(rtol=1e-6),reference_optimizer='original mixSQP with eps=0'),primitive,list(expected_weights=expected(f$fitted_g$pi)))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),SQUAREM=as.character(packageVersion('SQUAREM')),REBayes_available=requireNamespace('REBayes',quietly=TRUE),original_mixIP_executed=FALSE),cases=cases,em_fit_cases=direct),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
