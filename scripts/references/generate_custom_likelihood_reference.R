#!/usr/bin/env Rscript
# Original ashr calculations with user-defined arithmetic functions only.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_custom_likelihood_reference.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
hex <- function(x)paste(format(writeBin(as.double(x),raw(),size=8,endian='big')),collapse='')
bits <- function(x)lapply(as.numeric(x),function(v)paste0('bits:',hex(v)))
expected <- function(x)lapply(as.numeric(x),function(v)if(is.na(v)&&!is.nan(v))'NA'else sprintf('%.17g',v))
value <- function(x)list(type=if(is.logical(x))'logical'else'number',integer=is.integer(x),values=if(is.logical(x))as.list(x)else lapply(as.numeric(x),hex),dimensions=NULL,names=NULL,dimnames=list(NULL,NULL),attributes=list())
ast <- function(x){
 if(is.symbol(x))return(list(type='symbol',name=as.character(x)))
 if(!is.call(x))return(list(type='literal',value=value(x)))
 args <- as.list(x)[-1L];labels <- names(args)
 list(type='call',function_name=as.character(x[[1L]]),arguments=lapply(seq_along(args),function(i)list(name=if(is.null(labels)||labels[[i]]=='')NULL else labels[[i]],expression=ast(args[[i]]))))
}
function_descriptor <- function(f,constants=list()){
 body <- body(f)
 if(is.call(body)&&identical(body[[1L]],as.name('{')))body <- body[[2L]]
 list(arguments=as.list(names(formals(f))),expression=ast(body),constants=lapply(constants,value))
}
# Captured scalar shifts/scales demonstrate closure-data import without calling
# closures in the candidate wrapper. The original reference does call them.
make_logistic <- function(scale,shift) {
 list(name='custom_logistic',lcdfFUN=function(x)-log1p(exp(-(x-shift)/scale)),
      lpdfFUN=function(x)-(x-shift)/scale-2*log1p(exp(-(x-shift)/scale))-log(scale))
}
make_laplace <- function(scale) {
 list(name='custom_laplace',lcdfFUN=function(x)ifelse(x<0,log(.5)+x/scale,log1p(-.5*exp(-x/scale))),
      lpdfFUN=function(x)-abs(x)/scale-log(2*scale))
}
profiles <- list(logistic=list(lik=make_logistic(1,0),constants=list(scale=1,shift=0)),
 shifted_logistic=list(lik=make_logistic(.7,.2),constants=list(scale=.7,shift=.2)),
 laplace=list(lik=make_laplace(.8),constants=list(scale=.8)))
x <- c(-4,-1,0,.3,2,5);s <- c(.3,.7,1,2,.4,1)
g <- ashr::unimix(c(.2,.3,.5),c(0,-.4,-3),c(0,.4,3));cutoffs <- c(-1,0,1)
cases <- lapply(names(profiles),function(name){
 p <- profiles[[name]];lik <- p$lik
 descriptor <- list(name=lik$name,log_cdf=function_descriptor(lik$lcdfFUN,p$constants),log_density=function_descriptor(lik$lpdfFUN,p$constants),first=NULL,second=NULL)
 warnings <- character()
 fit <- withCallingHandlers(ashr::ash(x,s,g=g,fixg=TRUE,lik=lik,outputlevel=3),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 data <- ashr::set_data(x,s,lik=lik)
 cdf <- withCallingHandlers(vapply(cutoffs,function(c)ashr::cdf_post(g,c,data),numeric(length(x))),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 list(name=name,custom_likelihood=descriptor,estimates=bits(x),standard_errors=bits(s),thresholds=bits(cutoffs),
  family='uniform',point_mass=TRUE,nullweight=bits(10)[[1L]],scales=NULL,mode=bits(0)[[1L]],range=bits(c(-Inf,Inf)),
  alpha=bits(0)[[1L]],prune_threshold=bits(1e-10)[[1L]],observation_weights=NULL,fixed_prior=TRUE,df=NULL,
  supplied_prior=list(component_class='unimix',weights=bits(g$pi),means=NULL,sds=NULL,lower=bits(g$a),upper=bits(g$b)),
  posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD),negative_probability=expected(fit$result$NegativeProb),positive_probability=expected(fit$result$PositiveProb),
  local_false_discovery_rate=expected(fit$result$lfdr),local_false_sign_rate=expected(fit$result$lfsr),svalues=expected(fit$result$svalue),qvalues=expected(fit$result$qvalue),cdf=expected(cdf),log_likelihood=expected(fit$loglik)[[1L]],warnings=as.list(warnings))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),input_encoding='f64_bits',cdf_encoding='cutoff_major'),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
