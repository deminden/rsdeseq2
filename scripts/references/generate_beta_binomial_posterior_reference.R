#!/usr/bin/env Rscript
# Independent original apeglm outputs. Never load candidate rsdeseq2.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
hex<-function(x)paste(format(writeBin(as.double(x),raw(),size=8L,endian='big')),collapse='')
bits<-function(x)lapply(as.numeric(x),function(x)paste0('bits:',hex(x)))
expected<-function(x)lapply(as.numeric(x),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
mat<-function(x)if(is.null(x))NULL else list(rows=nrow(x),columns=ncol(x),values=bits(x))
x<-cbind(Intercept=1,condition=rep(c(0,1),each=4))
n<-matrix(c(10,12,15,18,20,22,25,30),4,8,byrow=TRUE)
y<-rbind(rep(0,8),n[2,],c(1,2,3,2,14,15,20,25),c(5,6,7,9,10,11,12,15))
weights<-matrix(rep(c(.5,1,2,1,.3,1,1.5,1),4),4,8,byrow=TRUE)
profiles<-list(base=list(),weighted=list(weights=weights),
  small_theta=list(theta=.05),large_theta=list(theta=1000),
  ignored_offset=list(offset=matrix(NA_real_,4,8)),
  no_trials=list(n=matrix(0,4,8),y=matrix(0,4,8)),
  three_coefficients=list(x=cbind(x,continuous=c(-1,0,.2,.8,.7,-.2,1,-.5))),
  adapted=list(mle=cbind(c(-1,.2,.5,2),rep(.4,4))),
  intercept=list(coefficient=1,nuisance=2,mle=cbind(c(-1,.2,.5,2),rep(.4,4))),
  all_gaussian=list(nuisance=c(1,2)),no_nuisance=list(nuisance=integer()),
  bounded=list(optimizer='L-BFGS-B',bounds=c(-.5,.5)),
  zero_weights=list(weights=matrix(0,4,8)))
methods<-c('betabinR','betabinCR','betabinC')
cases<-list()
for(profile in names(profiles))for(method in methods){
  p<-profiles[[profile]]
  Y<-if(is.null(p[['y']]))y else p[['y']]
  N<-if(is.null(p[['n']]))n else p[['n']]
  X<-if(is.null(p[['x']]))x else p[['x']]
  theta<-if(is.null(p[['theta']]))8 else p[['theta']]
  nuisance<-if(is.null(p[['nuisance']]))1 else p[['nuisance']]
  coefficient<-if(is.null(p[['coefficient']]))2 else p[['coefficient']]
  optimizer<-if(is.null(p[['optimizer']]))'BFGS'else p[['optimizer']]
  bounds<-if(is.null(p[['bounds']]))c(-Inf,Inf)else p[['bounds']]
  param<-cbind(rep(theta,nrow(Y)),N)
  prior<-list(no.shrink=nuisance,prior.mean=0,prior.scale=.7,prior.df=1,
    prior.no.shrink.mean=0,prior.no.shrink.scale=15)
  record<-list(name=paste(profile,method,sep='_'),observations=mat(Y),design=mat(X),parameters=mat(param),
    weights=mat(p[['weights']]),offsets=mat(p[['offset']]),mle=mat(p[['mle']]),coefficient=coefficient-1L,
    nuisance=as.list(nuisance-1L),scale=bits(.7)[[1L]],nuisance_scale=bits(15)[[1L]],
    prior_df=bits(1)[[1L]],prior_mean=bits(0)[[1L]],nuisance_mean=bits(0)[[1L]],
    method=method,optimizer=optimizer,bounds=bits(bounds),log_link=TRUE,threshold=bits(.35)[[1L]],likelihood=NULL)
  warnings<-character()
  result<-tryCatch({
    fit<-withCallingHandlers(apeglm::apeglm(Y=Y,x=X,param=param,coef=coefficient,mle=p[['mle']],
      weights=p[['weights']],offset=p[['offset']],prior.control=prior,method=method,
      optim.method=optimizer,bounds=bounds,threshold=.35),warning=function(w){
        warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(mode=expected(fit$map),sd=expected(fit$sd),fsr=expected(fit$fsr),sign_svalues=expected(fit$svalue),
      fsos=expected(fit$thresh),threshold_svalues=expected(apeglm::svalue(fit$thresh)),
      fitted_scale=expected(fit$prior.control$prior.scale)[[1L]],
      prior_variance=if(is.null(fit$prior.control$prior.var))NULL else expected(fit$prior.control$prior.var)[[1L]])
  },error=function(e)list(expected_error=conditionMessage(e)))
  cases[[length(cases)+1L]]<-c(record,result,list(warnings=as.list(warnings)))
}
# Frozen analytic likelihood/score observations catch policy and derivative
# mistakes independently of a posterior optimizer or its stopping rule.
leaves<-list()
for(theta in c(.01,8,1000))for(beta in list(c(0,0),c(-5,2),c(3,-1))){
  prior<-list(no.shrink=1,prior.mean=0,prior.scale=.7,prior.df=1,prior.no.shrink.mean=0,prior.no.shrink.scale=15)
  a<-list(beta=beta,x=x,y=y[3,],size=n[3,],theta=theta,weights=weights[3,],sigma=15,S=.7,no.shrink=1,shrink=2,cnst=0)
  leaves[[length(leaves)+1L]]<-list(name=paste(theta,paste(beta,collapse=','),sep='_'),
    design=mat(x),successes=bits(y[3,]),trials=bits(n[3,]),weights=bits(weights[3,]),theta=bits(theta)[[1L]],
    beta=bits(beta),scale=bits(.7)[[1L]],nuisance_scale=bits(15)[[1L]],nuisance=list(0L),
    value=expected(do.call(apeglm:::betabinFn,a))[[1L]],gradient=expected(do.call(apeglm:::betabinGr,a)))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm')),
  input_encoding='f64_bits'),cases=cases,analytic_leaves=leaves),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
