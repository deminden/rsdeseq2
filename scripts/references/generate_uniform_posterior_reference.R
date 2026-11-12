#!/usr/bin/env Rscript
# Independent original ashr references. Never load rsdeseq2.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_uniform_posterior_reference.R FRESH_OUTPUT_JSON')
if(file.exists(args[[1L]]))stop('Use a fresh path; references are never overwritten')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63')
numbers<-function(x)as.list(sprintf('%.17g',as.numeric(x)))
bits<-function(x)lapply(unname(as.numeric(x)),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
profiles<-list(
 broad=list(w=1,a=-100,b=100),
 single_interval=list(w=1,a=-1,b=1),
 point_zero=list(w=1,a=0,b=0),
 point_shifted=list(w=1,a=1.25,b=1.25),
 symmetric=list(w=c(.1,.3,.4,.2),a=c(0,-.1,-2,-100),b=c(0,.1,2,100)),
 halfuniform=list(w=c(.05,.2,.3,.3,.15),a=c(0,-.2,-3,0,0),b=c(0,0,0,.2,100)),
 shifted=list(w=c(.1,.3,.4,.2),a=c(-1.5,-2,1,-100),b=c(-1.5,-.2,3,100)),
 zero_weights=list(w=c(0,.4,0,.6),a=c(0,-.5,4,-100),b=c(0,.5,5,100)),
 narrow=list(w=c(.2,.3,.5),a=c(-1e-6,.5,-100),b=c(1e-6,.500002,100)),
 ultra_narrow=list(w=c(.4,.6),a=c(-1e-12,-100),b=c(1e-12,100)))
beta<-c(-40,-12,-6,-2,-.5,-1e-8,0,1e-8,.5,2,6,12,40)
se<-c(.01,.15,.7,2,8)
x<-rep(beta,each=length(se));s<-rep(se,times=length(beta))
cutoffs<-c(-Inf,-40,-2,-.5,0,.5,1.25,2,40,Inf)
cases<-lapply(names(profiles),function(name){
 p<-profiles[[name]];g<-ashr::unimix(p$w,p$a,p$b);data<-ashr::set_data(x,s)
 list(name=name,weights=bits(p$w),lower=bits(p$a),upper=bits(p$b),
  estimates=bits(x),standard_errors=bits(s),thresholds=bits(cutoffs),
  log_likelihoods=numbers(ashr:::log_comp_dens_conv.unimix(g,data)),
  posterior_weights=numbers(ashr:::comp_postprob(g,data)),posterior_mean=numbers(ashr::postmean(g,data)),posterior_sd=numbers(ashr::postsd(g,data)),
  cdf=numbers(t(vapply(cutoffs,function(c)ashr::cdf_post(g,c,data),numeric(length(x))))))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),input_encoding='f64_bits',platform=R.version$platform,
 dnorm_probe=numbers(stats::dnorm(37/.15,log=TRUE)),unfused_density_probe=numbers(-(.9189385332046728+.5*(37/.15)*(37/.15)))),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
