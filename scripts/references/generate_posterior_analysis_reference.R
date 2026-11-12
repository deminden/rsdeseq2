#!/usr/bin/env Rscript
# Independent original-package references; never loads candidate rsdeseq2.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
hex<-function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian='big')),collapse='')
bits<-function(v)lapply(as.numeric(v),function(v)paste0('bits:',hex(v)))
expected<-function(v)lapply(as.numeric(v),function(v)if(is.na(v)&&!is.nan(v))'NA'else sprintf('%.17g',v))
mat<-function(x)if(is.null(x))NULL else list(rows=nrow(x),columns=ncol(x),values=bits(x))
value<-function(v)list(type=if(is.logical(v))'logical'else'number',integer=is.integer(v),values=if(is.logical(v))as.list(v)else lapply(as.numeric(v),hex),dimensions=NULL,names=NULL,dimnames=list(NULL,NULL),attributes=list())
ast<-function(e){
  if(is.symbol(e))return(list(type='symbol',name=as.character(e)))
  if(!is.call(e))return(list(type='literal',value=value(e)))
  head<-e[[1L]];name<-if(is.symbol(head))as.character(head)else paste0(as.character(head[[2L]]),'::',as.character(head[[3L]]))
  inputs<-as.list(e)[-1L];labels<-names(inputs)
  list(type='call',function_name=name,arguments=lapply(seq_along(inputs),function(i)list(name=if(is.null(labels)||labels[[i]]=='')NULL else labels[[i]],expression=ast(inputs[[i]]))))
}
descriptor<-function(f)list(arguments=as.list(names(formals(f))),expression=ast(body(f)),constants=list())
normal<-function(y,x,beta,param,offset)stats::dnorm(y,mean=x%*%beta+offset,sd=param,log=TRUE)
beta_kernel<-function(y,x,beta,param,offset)lgamma(y+param[1]/(1+exp(-x%*%beta)))-lgamma(param[1]/(1+exp(-x%*%beta)))+lgamma(param[-1]-y+param[1]-param[1]/(1+exp(-x%*%beta)))-lgamma(param[1]-param[1]/(1+exp(-x%*%beta)))
Y<-rbind(c(-1,-.3,.2,.1,.6,1.2,.8,1),c(1,1.5,.9,1.3,-.1,.2,0,.3),c(0,.1,-.1,0,0,.1,-.1,0))
NB<-rbind(c(0,1,2,1,4,5,3,6),c(10,8,12,9,2,3,1,2),rep(0,8))
x<-cbind(Intercept=1,condition=rep(c(0,1),each=4))
x3<-cbind(x,continuous=c(-1,0,.2,.8,.7,-.2,1,-.5))
pc<-list(no.shrink=1,prior.mean=0,prior.scale=.7,prior.df=1,prior.no.shrink.mean=0,prior.no.shrink.scale=15)
profiles<-list(
 normal=list(),correlated=list(x=x3),
 weighted=list(weights=matrix(rep(c(.5,1,2,1,.3,1,1.5,1),3),3,8,byrow=TRUE),offset=matrix(rep(seq(-.2,.2,length.out=8),3),3,8,byrow=TRUE)),
 narrow=list(nsd=1),fine=list(ngrid=100L,ngrid.nuis=3L),map_nuisance=list(ngrid.nuis=1L),
 uncertain=list(param.sd=rep(.1,3)),uncertain_nuisance_one=list(param.sd=rep(.1,3),ngrid.nuis=1L),
 all_gaussian=list(no.shrink=1:2),empty_nuisance=list(no.shrink=integer()),
 student_prior=list(prior.df=4),negative_binomial=list(method='nbinomCR',Y=NB,param=rep(.2,3)),
 grid_override=list(method='nbinomCR',Y=NB,param=rep(.2,3),f=function(y,x,beta,param,offset).9*stats::dnbinom(y,mu=exp(x%*%beta+offset),size=1/param,log=TRUE)),
 beta_binomial=list(method='betabinCR',Y=rbind(c(1,2,3,2,14,15,20,25),c(5,6,7,9,10,11,12,15),rep(0,8)),param=cbind(8,matrix(rep(c(10,12,15,18,20,22,25,30),3),3,8,byrow=TRUE))))
cases<-list()
for(name in names(profiles))for(kind in c('laplace','HPD','credible')){
 p<-profiles[[name]]
 if(kind=='laplace'&&!is.null(p[['param.sd']]))next
 getp<-function(key,default)if(is.null(p[[key]]))default else p[[key]]
 X<-getp('x',x);y<-getp('Y',Y);method<-getp('method','general');param<-getp('param',rep(.7,3))
 f<-getp('f',if(method=='nbinomCR')apeglm::logLikNB else if(method=='betabinCR')beta_kernel else normal)
 prior<-pc
 for(key in c('no.shrink','prior.df'))if(!is.null(p[[key]]))prior[[key]]<-p[[key]]
 C<-if(ncol(X)==3)cbind(c(0,1,.5),c(-1,1,0),c(0,0,0))else cbind(c(0,1),c(-1,1),c(0,0))
 offset<-getp('offset',matrix(0,nrow(y),ncol(y)))
 o<-list(kind=kind,level=bits(.9)[[1]],grid_intervals=getp('ngrid',50L),grid_standard_deviations=bits(getp('nsd',5))[[1]],nuisance_points=getp('ngrid.nuis',5L),nuisance_standard_deviations=bits(2)[[1]],contrasts=mat(C),parameter_standard_errors=if(is.null(p[['param.sd']]))NULL else bits(p[['param.sd']]))
 warnings<-character()
 result<-tryCatch({fit<-withCallingHandlers(apeglm::apeglm(Y=y,x=X,log.lik=f,param=param,coef=2,method=method,offset=offset,
   weights=p[['weights']],prior.control=prior,contrasts=C,interval.type=kind,interval.level=.9,ngrid=o$grid_intervals,nsd=getp('nsd',5),ngrid.nuis=o$nuisance_points,nsd.nuis=2,param.sd=p[['param.sd']]),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
   repairs<-list()
   if(name=='beta_binomial'&&kind=='credible') {
     # Original zero-index interval assembly drops the lower endpoint and
     # recycles the upper. Independently invert its fitted exponential tail.
     g<-3L;grid<-seq(fit$map[g,2]-5*fit$sd[g,2],fit$map[g,2]+5*fit$sd[g,2],length.out=51)
     unpo<-apeglm:::integrateOther(cor.mat=diag(2),coef=2,betas=grid,map=fit$map[g,],sd=fit$sd[g,],log.lik=f,log.prior=apeglm:::log.prior,y=y[g,],x=X,param=param[g,],param.sd=NULL,prior.control=prior,weights=rep(1,8),offset=offset[g,],ngrid.nuis=5,nsd.nuis=2)
     slope<-diff(log(unpo[1:2]))/diff(grid[1:2]);area<-fit$diag[g,4]
     # Solve area*exp(slope*(q-grid[1]))=.05 via an independent root finder.
     lower<-uniroot(function(q)area*exp(slope*(q-grid[1]))-.05,c(grid[1]-10,grid[1]),tol=1e-12)$root
     stopifnot(fit$interval[g,1]==fit$interval[g,2],lower<grid[1])
     repairs<-list(list(gene=g-1L,column=0L,expected=expected(lower)[[1L]],reason='original lower endpoint dropped and upper recycled'))
   }
   list(interval_repairs=repairs,mode=expected(fit$map),sd=expected(fit$sd),fsr=expected(fit$fsr),sign_svalues=expected(fit$svalue),intervals=expected(fit$interval),contrast_mode=expected(fit$contrast.map),contrast_sd=expected(fit$contrast.sd),outside=if(kind=='laplace')NULL else expected(fit$diag[,4:5]))
 },error=function(e)list(expected_error=conditionMessage(e)))
 cases[[length(cases)+1L]]<-c(list(name=paste(name,kind,sep='_'),observations=mat(y),design=mat(X),parameters=mat(as.matrix(param)),offsets=mat(offset),weights=mat(p[['weights']]),coefficient=1L,nuisance=as.list(prior$no.shrink-1L),scale=bits(.7)[[1L]],nuisance_scale=bits(15)[[1L]],prior_df=bits(prior$prior.df)[[1L]],prior_mean=bits(0)[[1L]],nuisance_mean=bits(0)[[1L]],method=method,optimizer='BFGS',bounds=bits(c(-Inf,Inf)),log_link=TRUE,threshold=NULL,likelihood=if(method=='general'||method=='betabinCR'||name=='grid_override')descriptor(f)else NULL,analysis=o),result,list(warnings=as.list(warnings)))
}
# Fixed-mode queries isolate grid arithmetic from optimizer and Hessian drift.
# The likelihood is Gaussian in beta[2]; its nuisance contribution cancels.
leaves<-list()
for(mean in c(-1.3,0,.3,2.5))for(kind in c('HPD','credible'))for(nsd in c(1,5)){
 f<-local({m<-mean;function(y,x,beta,param,offset)dnorm(beta[2],m,.7,log=TRUE)})
 prior<-pc;prior$no.shrink<-1:2;prior$prior.no.shrink.scale<-1e8
 out<-apeglm:::gridResults(y=0,x=matrix(0,1,2),log.lik=f,param=1,coef=2,interval.type=kind,interval.level=.9,threshold=.2,contrasts=NULL,weights=1,offset=0,flip.sign=TRUE,prior.control=prior,ngrid=50,nsd=nsd,ngrid.nuis=5,nsd.nuis=2,log.link=FALSE,param.sd=NULL,o=list(convergence=0,counts=c(1,1),value=0),map=c(0,mean),cov.mat=diag(c(1,.7)^2),sd=c(1,.7),out=list())
 leaves[[length(leaves)+1L]]<-list(name=paste(mean,kind,nsd,sep='_'),mean=mean,kind=kind,nsd=nsd,interval=expected(out$ci),fsr=expected(out$fsr),threshold=expected(if(mean<0)1-out$thresh else out$thresh),outside=expected(out$diag[4:5]))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm'))),cases=cases,grid_leaves=leaves),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
