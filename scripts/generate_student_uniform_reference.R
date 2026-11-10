#!/usr/bin/env Rscript
# Independent original-package t-likelihood references. Never load rsdeseq2.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_student_uniform_reference.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
bits <- function(x)lapply(as.numeric(x),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
numbers <- function(x)as.list(sprintf('%.17g',as.numeric(x)))
profiles <- list(symmetric=list(w=c(.1,.4,.5),a=c(0,-.3,-4),b=c(0,.3,4)),
  halfuniform=list(w=c(.1,.3,.6),a=c(0,-3,0),b=c(0,0,3)),
  shifted=list(w=c(.3,.7),a=c(-2,.3),b=c(-.5,4)),
  point=list(w=1,a=1.25,b=1.25),
  narrow=list(w=c(.2,.8),a=c(-1e-6,-10),b=c(1e-6,10)))
x <- rep(c(-15,-4,-.5,0,.5,4,15),each=3L);s <- rep(c(.1,1,4),times=7L)
thresholds <- c(-Inf,-2,0,.3,1.25,2,Inf)
cases <- list()
for(df in c(.5,1,2,3,10,100))for(name in names(profiles)){
 p <- profiles[[name]];g <- ashr::unimix(p$w,p$a,p$b)
 data <- ashr::set_data(x,s,ashr::lik_t(df))
 warnings <- character()
 result <- withCallingHandlers(list(log_likelihoods=numbers(ashr:::log_comp_dens_conv.unimix(g,data)),
   posterior_weights=numbers(ashr:::comp_postprob(g,data)),posterior_mean=numbers(ashr::postmean(g,data)),
   posterior_sd=numbers(ashr::postsd(g,data)),
   cdf=numbers(t(vapply(thresholds,function(c)ashr::cdf_post(g,c,data),numeric(length(x)))))),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 cases[[length(cases)+1L]] <- c(list(name=paste(name,df,sep='_df'),df=bits(df)[[1L]],weights=bits(p$w),
  lower=bits(p$a),upper=bits(p$b),estimates=bits(x),standard_errors=bits(s),thresholds=bits(thresholds)),result,
  list(reference_warnings=as.list(unique(warnings))))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),
  etrunct=as.character(packageVersion('etrunct')),input_encoding='f64_bits'),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
