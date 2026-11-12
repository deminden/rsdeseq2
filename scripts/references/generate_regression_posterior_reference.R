#!/usr/bin/env Rscript
# Original apeglm/statistical expressions only; candidate code is never loaded.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
hex<-function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian='big')),collapse='')
bits<-function(v)lapply(as.numeric(v),function(v)paste0('bits:',hex(v)))
expected<-function(v)lapply(as.numeric(v),function(v)if(is.na(v)&&!is.nan(v))'NA'else sprintf('%.17g',v))
matrix_input<-function(x)if(is.null(x))NULL else list(rows=nrow(x),columns=ncol(x),values=bits(x))
value<-function(v)list(type=if(is.logical(v))'logical'else'number',integer=is.integer(v),values=if(is.logical(v))as.list(v)else lapply(as.numeric(v),hex),dimensions=NULL,names=NULL,dimnames=list(NULL,NULL),attributes=list())
ast<-function(e){
  if(is.symbol(e))return(list(type='symbol',name=as.character(e)))
  if(!is.call(e))return(list(type='literal',value=value(e)))
  head<-e[[1L]];name<-if(is.symbol(head))as.character(head)else paste0(as.character(head[[2L]]),'::',as.character(head[[3L]]))
  inputs<-as.list(e)[-1L];labels<-names(inputs)
  list(type='call',function_name=name,arguments=lapply(seq_along(inputs),function(i)list(name=if(is.null(labels)||labels[[i]]=='')NULL else labels[[i]],expression=ast(inputs[[i]]))))
}
descriptor<-function(f,constants=list())list(arguments=as.list(names(formals(f))),expression=ast(body(f)),constants=lapply(constants,value))
normal<-function(y,x,beta,param,offset)stats::dnorm(y,mean=x%*%beta+offset,sd=param,log=TRUE)
noise<-.7
algebra<-function(y,x,beta,param,offset)-.5*((y-x%*%beta-offset)/noise)^2-log(noise)-.5*log(2*pi)
student<-function(y,x,beta,param,offset)stats::dt((y-x%*%beta-offset)/param,df=4,log=TRUE)-log(param)
negative_binomial<-function(y,x,beta,param,offset)stats::dnbinom(y,mu=exp(x%*%beta+offset),size=1/param,log=TRUE)
poisson<-function(y,x,beta,param,offset)y*(x%*%beta+offset)-exp(x%*%beta+offset)-lgamma(y+1)
numeric_product<-function(y,x,beta,param,offset)stats::dnorm(y,mean=as.numeric(x%*%beta)+offset,sd=param,log=TRUE)
bad_shape<-function(y,x,beta,param,offset)stats::dnorm(y,mean=x%*%param,sd=1,log=TRUE)
Y<-rbind(c(-1,-.3,.2,.1,.6,1.2,.8,1),c(1,1.5,.9,1.3,-.1,.2,0,.3),
  c(0,.1,-.1,0,0,.1,-.1,0),c(2,2.2,1.8,2.1,3,3.5,3.1,2.8))
NB<-rbind(rep(0,8),c(0,1,2,1,4,5,3,6),c(10,8,12,9,2,3,1,2),c(2,1,3,2,2,1,3,2))
x<-cbind(Intercept=1,condition=rep(c(0,1),each=4))
offset<-matrix(rep(seq(-.2,.2,length.out=8),4),nrow=4,byrow=TRUE)
weights<-matrix(rep(c(.5,1,2,1,.3,1,1.5,1),4),nrow=4,byrow=TRUE)
pc<-list(no.shrink=1,prior.mean=0,prior.scale=1,prior.df=1,prior.no.shrink.mean=0,prior.no.shrink.scale=15)
profiles<-list(normal=list(),algebra=list(f=algebra,constants=list(noise=noise,pi=pi)),
  weighted_offset=list(weights=weights,offset=offset),student_noise=list(f=student),
  custom_nb=list(f=negative_binomial,Y=NB,param=rep(.2,4)),stock_nb=list(stock=TRUE,Y=NB,param=rep(.2,4)),
  nbinomR=list(method='nbinomR',Y=NB,param=rep(.2,4)),
  nbinomR_no_log=list(method='nbinomR',Y=NB,param=rep(.2,4),log_link=FALSE),
  bounded_nb=list(method='nbinomR',Y=NB,param=rep(.2,4),optimizer='L-BFGS-B',bounds=c(-10,10)),
  bounded_normal=list(optimizer='L-BFGS-B',bounds=c(-.5,.5)),
  student_prior=list(prior_df=4,scale=.7,nuisance_mean=.3),
  noncentral_prior=list(prior_df=3,prior_mean=.7,scale=.7),normal_prior=list(prior_df=Inf),
  adapted_intercept=list(coef=1,mle=cbind(c(-1,0,1,2),rep(.4,4))),
  empty_nuisance=list(nuisance=integer()),duplicate_nuisance=list(nuisance=c(1,1)),
  numeric_product=list(f=numeric_product),poisson=list(f=poisson,Y=NB,param=rep(1,4)),
  unused_missing_param=list(f=algebra,constants=list(noise=noise,pi=pi),param=rep(NA_real_,4)),
  invalid_shape=list(f=bad_shape),null_param=list(null_param=TRUE))
cases<-lapply(names(profiles),function(name){
  p<-profiles[[name]];f<-if(isTRUE(p[["stock", exact=TRUE]]))apeglm::logLikNB else if(is.null(p[["f", exact=TRUE]]))normal else p[["f", exact=TRUE]]
  y<-if(is.null(p[["Y", exact=TRUE]]))Y else p[["Y", exact=TRUE]]
  param<-if(isTRUE(p[["null_param", exact=TRUE]]))NULL else if(is.null(p[["param", exact=TRUE]]))rep(.7,4)else p[["param", exact=TRUE]]
  prior<-pc
  for(pair in list(c('prior_df','prior.df'),c('prior_mean','prior.mean'),c('scale','prior.scale'),c('nuisance_mean','prior.no.shrink.mean')))
    if(!is.null(p[[pair[[1L]]]]))prior[[pair[[2L]]]]<-p[[pair[[1L]]]]
  if(!is.null(p[["nuisance", exact=TRUE]]))prior$no.shrink<-p[["nuisance", exact=TRUE]]
  coefficient<-if(is.null(p[["coef", exact=TRUE]]))2 else p[["coef", exact=TRUE]]
  method<-if(is.null(p[["method", exact=TRUE]]))'general'else p[["method", exact=TRUE]]
  optimizer<-if(is.null(p[["optimizer", exact=TRUE]]))'BFGS'else p[["optimizer", exact=TRUE]]
  bounds<-if(is.null(p[["bounds", exact=TRUE]]))c(-Inf,Inf)else p[["bounds", exact=TRUE]]
  log_link<-if(is.null(p[["log_link", exact=TRUE]]))TRUE else p[["log_link", exact=TRUE]]
  record<-list(name=name,observations=matrix_input(y),design=matrix_input(x),parameters=matrix_input(if(is.null(param))NULL else as.matrix(param)),
    weights=matrix_input(p[["weights", exact=TRUE]]),offsets=matrix_input(p[["offset", exact=TRUE]]),mle=matrix_input(p[["mle", exact=TRUE]]),coefficient=coefficient-1L,
    nuisance=as.list(prior$no.shrink-1L),scale=bits(prior$prior.scale)[[1L]],nuisance_scale=bits(prior$prior.no.shrink.scale)[[1L]],
    prior_df=bits(prior$prior.df)[[1L]],prior_mean=bits(prior$prior.mean)[[1L]],nuisance_mean=bits(prior$prior.no.shrink.mean)[[1L]],
    method=method,optimizer=optimizer,bounds=bits(bounds),log_link=log_link,threshold=bits(.35)[[1L]],
    likelihood=if(isTRUE(p[["stock", exact=TRUE]])||method=='nbinomR')NULL else descriptor(f,if(is.null(p[["constants", exact=TRUE]]))list()else p[["constants", exact=TRUE]]))
  warnings<-character()
  result<-tryCatch({
    fit<-withCallingHandlers(apeglm::apeglm(Y=y,x=x,log.lik=f,param=param,coef=coefficient,mle=p[["mle", exact=TRUE]],
      weights=p[["weights", exact=TRUE]],offset=if(is.null(p[["offset", exact=TRUE]]))matrix(0,4,8)else p[["offset", exact=TRUE]],threshold=.35,prior.control=prior,
      log.link=log_link,method=method,optim.method=optimizer,bounds=bounds),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(mode=expected(fit$map),sd=expected(fit$sd),fsr=expected(fit$fsr),sign_svalues=expected(fit$svalue),
      fsos=expected(fit$thresh),threshold_svalues=expected(apeglm::svalue(fit$thresh)),
      fitted_scale=expected(fit$prior.control$prior.scale)[[1L]],prior_variance=if(is.null(fit$prior.control$prior.var))NULL else expected(fit$prior.control$prior.var)[[1L]])
  },error=function(e)list(expected_error=conditionMessage(e)))
  c(record,result,list(warnings=as.list(warnings)))
})
q<-c(-8,-2,-.1,0,.1,2,8)
leaves<-list()
add_leaf<-function(name,f,y,constants=list(),parameters=1){
  leaves[[length(leaves)+1L]]<<-list(name=name,likelihood=descriptor(f,constants),observations=bits(y),parameters=bits(parameters),
    expected=expected(f(y,matrix(0,length(y),1),0,parameters,rep(0,length(y)))))
}
add_leaf('normal',function(y,x,beta,param,offset)stats::dnorm(y,mean=.3,sd=.7,log=TRUE),c(q,NA_real_,NaN,-Inf,Inf))
add_leaf('nb_large_size',function(y,x,beta,param,offset)stats::dnbinom(y,size=1e8,mu=4,log=TRUE),c(0,1,4,20,1000))
add_leaf('nb_large_count',function(y,x,beta,param,offset)stats::dnbinom(y,size=30,mu=1e6,log=TRUE),c(0,1,1e5,1e6,1e7))
add_leaf('nb_probability',function(y,x,beta,param,offset)stats::dnbinom(y,size=3,prob=.4,log=TRUE),c(0,1,2,10,100))
add_leaf('binomial',function(y,x,beta,param,offset)stats::dbinom(y,size=1000,prob=.2,log=TRUE),c(0,1,199,200,201,1000))
for(df in c(.5,1,3,30,100,1e10,Inf))for(ncp in c(0,.7)){
  f<-local({d<-df;n<-ncp;function(y,x,beta,param,offset)stats::dt(y,df=d,ncp=n,log=TRUE)})
  add_leaf(paste0('student_',df,'_',ncp),f,q,list(d=df,n=ncp))
}
add_leaf('unused_branch',function(y,x,beta,param,offset)ifelse(y>0,stats::dnorm(y,log=TRUE),x%*%param),rep(1,3),parameters=c(1,2))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm')),input_encoding='f64_bits'),cases=cases,leaves=leaves),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
