#!/usr/bin/env Rscript
# Original ashr only; independent reference for the native convex weight fitter.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_mixture_fit_reference.R OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63')
numbers<-function(x)as.list(sprintf('%.17g',as.numeric(x)))
set.seed(1763)
x<-c(rnorm(128,0,.15),rnorm(64,0,1.3),rnorm(64,0,4))
s<-rep(c(.15,.3,.6,1),length.out=length(x))
base<-list(means=rep(0,7),sds=c(0,.03,.15,.5,1.5,4,10))
profiles<-list(
  null=list(x=seq(-.1,.1,length.out=128),s=rep(.7,128)),
  mixed=list(x=x,s=s),
  weighted=list(x=x,s=s,weights=rep(c(0,.2,1,2),length.out=length(x))),
  nullbiased=list(x=x,s=s,prior=c(10,rep(1,6))),
  shifted=list(x=x+rep(c(-1,0,2,3),length.out=length(x)),s=s,
    means=c(-2,-.5,0,1,3),sds=c(.1,.2,0,.8,3)),
  duplicate=list(x=x,s=s,means=rep(0,6),sds=c(0,.15,.5,.5,2,5)),
  near_collinear=list(x=x,s=s,means=rep(0,6),sds=c(0,.15,.5,.500001,2,5)),
  active_duplicate=list(x=x,s=s,means=rep(0,7),sds=c(0,.15,.15,.5,1.5,4,10)))
cases<-lapply(names(profiles),function(name){
 p<-utils::modifyList(base,profiles[[name]])
 k<-length(p$sds);g<-ashr::normalmix(rep(1/k,k),p$means,p$sds)
 prior<-if(is.null(p$prior))rep(1,k)else p$prior
 weights<-if(is.null(p$weights))rep(1,length(p$x))else p$weights
 f<-ashr::ash(p$x,p$s,g=g,prior=prior,optmethod='mixSQP',weights=weights,outputlevel=3)
 data<-ashr::set_data(p$x,p$s)
 log_likelihoods<-ashr:::log_comp_dens_conv.normalmix(g,data)
 list(name=name,observations=length(p$x),components=k,means=numbers(p$means),sds=numbers(p$sds),
   estimates=numbers(p$x),standard_errors=numbers(p$s),dirichlet=numbers(prior),observation_weights=numbers(weights),
   log_likelihoods=numbers(log_likelihoods),fitted_weights=numbers(f$fitted_g$pi),
   weighted_log_likelihood=numbers(sum(weights*ashr:::calc_vloglik(f$fitted_g,data))),
   posterior_mean=numbers(ashr::postmean(f$fitted_g,data)),posterior_sd=numbers(ashr::postsd(f$fitted_g,data)),
   converged=f$fit_details$optreturn$converged)
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),optimizer='mixSQP default controls'),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
