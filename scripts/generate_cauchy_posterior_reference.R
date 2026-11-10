#!/usr/bin/env Rscript
# Independent original-package priors and negative-binomial posterior batches.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',
  as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
bits <- function(x)lapply(as.numeric(x),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
expected <- function(x)lapply(as.numeric(x),function(v)if(is.na(v)&&!is.nan(v))'NA' else sprintf('%.17g',v))
matrix_input <- function(x)if(is.null(x))NULL else list(rows=nrow(x),columns=ncol(x),values=bits(x))
decode <- function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')},numeric(1L))
i<-seq_len(128)
prior_profiles<-list(lower=list(x=rep(0,32),s=rep(.7,32)),
  interior=list(x=1.4*sin(i*.73)+.6*cos(i*.11),s=.2+(i%%11)/13),
  near_lower=list(x=c(.001,.002,.003,.004),s=rep(.0002,4)),
  cap_scale=list(x=4*sin(i*.4),s=rep(.2,128)),
  missing=list(x=c(NA,NaN,-2,-.3,0,.3,2),s=c(NA,Inf,.2,.7,1,.7,.2)),
  endpoint=list(x=rep(sqrt(400+.2^2),8),s=rep(.2,8)))
priors<-lapply(names(prior_profiles),function(name){
  p<-prior_profiles[[name]]
  list(name=name,estimates=bits(p$x),standard_errors=bits(p$s),
    variance=expected(apeglm:::priorVar(cbind(p$x,p$s)))[[1L]])
})
Y<-rbind(rep(0,8),c(0,1,0,1,1,0,0,2),c(1,0,2,1,0,1,0,0),
  c(10,20,30,20,40,30,60,50),c(200,180,230,210,80,90,100,75),
  c(10000,12000,9000,11000,15000,14000,16000,17000),rep(4,8),
  c(0,0,0,0,40,50,60,55),c(30,40,30,50,0,0,0,0),
  c(4,7,3,6,2,3,1,2),c(12,23,15,18,22,35,25,30),c(1,0,0,0,0,0,0,0))
x<-cbind(Intercept=1,condition=rep(c(0,1),each=4))
x3<-cbind(x,covariate=c(-1,-.5,.5,1,-1,-.5,.5,1))
disp<-c(NA,.2,.2,.1,.12,.03,.3,.2,.2,.4,.15,.3)
weights<-matrix(rep(c(1,.25,.5,2,1,.3,1,1),nrow(Y)),nrow=nrow(Y),byrow=TRUE)
weights[4,]<-0
offset<-outer(seq(-.2,.2,length.out=nrow(Y)),seq(-.3,.3,length.out=8),'+')
prior<-list(no.shrink=1,prior.mean=0,prior.scale=1,prior.df=1,
  prior.no.shrink.mean=0,prior.no.shrink.scale=15)
profiles<-list(default=list(),weights=list(weights=weights),offset=list(offset=offset),
  weighted_offset=list(weights=weights,offset=offset),
  three_coefficients=list(x=x3,coef=3),multiple_shrunk=list(x=x3,coef=3),
  narrow_prior=list(scale=.05),wide_prior=list(scale=2),
  narrow_nuisance=list(nuisance_scale=.5),
  adapted=list(mle=cbind(c(NA,-.2,.3,1,-1,.4,0,3,-3,-.7,.5,-.2),rep(.2,12))),
  no_flip=list(flip=FALSE,threshold=-.2),map_only=list(method='nbinomC'))
# Frozen input row from the existing seeded 1200-gene public fixture. It detects
# a mathematically more converged optimizer that violates the source MAP bound.
profiles$stopping_rule <- list(counts=matrix(c(0,6,3,0,4,3,10,11,12,7,5,11),nrow=1),
  x=cbind(Intercept=1,condition=rep(c(0,1),each=6)),
  dispersions=decode('bits:3fa0a721609b5e9d'),scale=decode('bits:3fe69c44127d7bb9'),
  offset=matrix(decode(c('bits:bfd1de98be43a580','bits:bfcd162b7525be21','bits:bfc4a68f88fb569e',
    'bits:bfbcfb332c715f00','bits:bfacd58b828026ff','bits:bf8aba27ddb879e5','bits:3f9ca77675808205',
    'bits:3fb49ac156ff0646','bits:3fc159b3a29247e1','bits:3fc8cbd0a6ac40c2','bits:3fcc6a1defcc33df',
    'bits:3fd144c866c27a31')),nrow=1))
cases<-lapply(names(profiles),function(name){
  p<-profiles[[name]];design<-if(is.null(p$x))x else p$x
  counts<-if(is.null(p$counts))Y else p$counts
  dispersions<-if(is.null(p$dispersions))disp else p$dispersions
  coefficient<-if(is.null(p$coef))2 else p$coef
  pc<-prior
  if(name=='three_coefficients')pc$no.shrink<-c(1,2)
  if(!is.null(p$scale))pc$prior.scale<-p$scale
  if(!is.null(p$nuisance_scale))pc$prior.no.shrink.scale<-p$nuisance_scale
  threshold<-if(is.null(p$threshold)).35 else p$threshold
  flip<-if(is.null(p$flip))TRUE else p$flip
  method<-if(is.null(p$method))'nbinomCR' else p$method
  warnings<-character()
  fit<-withCallingHandlers(apeglm::apeglm(Y=counts,x=design,param=dispersions,coef=coefficient,
    mle=p$mle,threshold=threshold,weights=p$weights,offset=p$offset,
    flip.sign=flip,prior.control=pc,method=method),warning=function(w){
      warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  list(name=name,counts=matrix_input(counts),design=matrix_input(design),dispersions=bits(dispersions),
    weights=matrix_input(p$weights),offsets=matrix_input(p$offset),
    mle=matrix_input(p$mle),coefficient=coefficient-1L,scale=bits(pc$prior.scale)[[1L]],
    nuisance=as.list(pc$no.shrink-1L),nuisance_scale=bits(pc$prior.no.shrink.scale)[[1L]],
    threshold=bits(threshold)[[1L]],flip_sign=flip,map_only=method=='nbinomC',
    prior_variance=if(is.null(fit$prior.control$prior.var))NULL else expected(fit$prior.control$prior.var)[[1L]],
    fitted_scale=expected(fit$prior.control$prior.scale)[[1L]],
    mode=expected(fit$map),sd=expected(fit$sd),fsr=expected(fit$fsr),
    sign_svalues=expected(fit$svalue),fsos=expected(fit$thresh),
    threshold_svalues=expected(apeglm::svalue(fit$thresh)),
    convergence=expected(fit$diag[,'conv']),warnings=as.list(warnings))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),
  apeglm=as.character(packageVersion('apeglm')),input_encoding='f64_bits'),
  priors=priors,cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
