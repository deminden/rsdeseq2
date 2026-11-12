#!/usr/bin/env Rscript
# Independent original likelihood kernels, including both sides of large means.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
cases<-list();z<-c(-12,-8,-3,0,3,8,12)
for(size in c(1e3,1e5,1e6,1e8,1e10))for(kind in c('pois','pois_log','binom','binom_logit','logF','logF_first_infinite','logF_second_infinite')){
  if(startsWith(kind,'pois')){
    link<-if(kind=='pois')'identity'else'log';q<-1+z/sqrt(size)
    a<-rep(-1-2/sqrt(size),length(z));b<-rep(-1+2/sqrt(size),length(z))
    if(link=='log'){q<--log(q);a[]<--log(1+2/sqrt(size));b[]<--log(1-2/sqrt(size))}
    lik<-ashr::lik_pois(size,size,link);model<-list(type='pois',y=bits(size),scale=bits(size),link=link)
  }else if(startsWith(kind,'binom')){
    link<-if(kind=='binom')'identity'else'logit';delta<-1/sqrt(8*size);q<-.5+z*delta
    a<-rep(-.5-2*delta,length(z));b<-rep(-.5+2*delta,length(z))
    if(link=='logit'){q<-qlogis(q);a[]<-qlogis(.5-2*delta);b[]<-qlogis(.5+2*delta)}
    lik<-ashr::lik_binom(size,2*size,link);model<-list(type='binom',y=bits(size),trials=bits(2*size),link=link)
  }else{
    d1<-if(kind=='logF_first_infinite')Inf else 2*size
    d2<-if(kind=='logF_second_infinite')Inf else 2*size
    delta<-sqrt(2/d1+2/d2);q<-z*delta;a<-rep(-2*delta,length(z));b<-rep(2*delta,length(z))
    lik<-ashr:::add_etruncFUN(ashr::lik_logF(d1,d2));model<-list(type='logF',df1=bits(d1),df2=bits(d2))
  }
  warnings<-character()
  values<-withCallingHandlers(list(log_cdf=expected(lik$lcdfFUN(q)),log_density=expected(lik$lpdfFUN(q)),
    first=expected(lik$etruncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),second=expected(lik$e2truncFUN(matrix(a,ncol=1),matrix(b,ncol=1)))),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  cases[[length(cases)+1L]]<-c(list(name=paste(kind,size,sep='_'),builtin_likelihood=model,arguments=bits(q),lower=bits(a),upper=bits(b)),values,list(warnings=as.list(warnings)))
}
# Highly unequal beta shapes exercise normalization without a central-limit profile.
for(y in c(0,3,1000))for(link in c('identity','logit')){
  n<-1e10;p<-(y+1)/(n+2);q<-p*c(.01,.1,.5,1,2,5,20)
  a<-rep(-2*p,length(q));b<-rep(-p/2,length(q))
  if(link=='logit'){q<-qlogis(q);a[]<-qlogis(p/2);b[]<-qlogis(2*p)}
  lik<-ashr::lik_binom(y,n,link)
  cases[[length(cases)+1L]]<-list(name=paste('unequal',y,link,sep='_'),builtin_likelihood=list(type='binom',y=bits(y),trials=bits(n),link=link),arguments=bits(q),lower=bits(a),upper=bits(b),
    log_cdf=expected(lik$lcdfFUN(q)),log_density=expected(lik$lpdfFUN(q)),first=expected(lik$etruncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),second=expected(lik$e2truncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),warnings=list())
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr'))),kernel_cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
