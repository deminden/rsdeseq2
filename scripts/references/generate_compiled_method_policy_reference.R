#!/usr/bin/env Rscript
# Original apeglm only; source-refinement input was selected without native output.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
decode<-function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(sub('bits:','',h),seq(1,15,2),seq(2,16,2)),16L)),'double',n=1,size=8,endian='big'),numeric(1))
mat<-function(m)if(is.null(m))NULL else matrix(decode(m$values),nrow=m$rows,ncol=m$columns)
descriptor<-function(m)list(rows=nrow(m),columns=ncol(m),values=bits(m))
state<-function(v)lapply(v,function(x)paste(format(writeBin(as.integer(x),raw(),size=4L,endian='big')),collapse=''))
base<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/cauchy_posterior_r462.json'),simplifyVector=FALSE)$cases[[1L]]
Y<-mat(base$counts);x<-mat(base$design);dispersion<-decode(base$dispersions)
nonzero<-rowSums(Y)>0
hard<-rbind(c(10000,8000,12000,0,1,0),c(1,1,1,10000,8000,12000),c(0,1,0,1,0,0))
hardx<-cbind(1,c(0,0,0,1,1,1))
.compiled_policy_capture<-new.env()
suppressMessages(trace('optimNbinom',where=asNamespace('apeglm'),print=FALSE,tracer=quote({.compiled_policy_capture$calls<-.compiled_policy_capture$calls+1L})))
profiles<-list(
 list(name='map_only_ignores_bounded_controls',method='nbinomC',optimizer='L-BFGS-B',bounds=c(-.01,.01),log_link=FALSE),
 list(name='map_only_ignores_unknown_controls',method='nbinomC',optimizer='not_an_optimizer',bounds=numeric(),log_link=FALSE),
 list(name='random_map_ignores_unknown_controls',method='nbinomC*',optimizer='not_an_optimizer',bounds=numeric(),log_link=FALSE),
 list(name='valid_prefit_ignores_tight_bounds',method='nbinomCR',optimizer='L-BFGS-B',bounds=c(-.01,.01),log_link=FALSE,nonzero=TRUE),
 list(name='valid_prefit_ignores_unknown_optimizer',method='nbinomCR',optimizer='not_an_optimizer',bounds=numeric(),log_link=TRUE,nonzero=TRUE),
 list(name='valid_prefit_defers_invalid_bounds',method='nbinomCR',optimizer='L-BFGS-B',bounds=numeric(),log_link=TRUE,nonzero=TRUE),
 list(name='false_loglink_zero_row_rejects',method='nbinomCR',optimizer='BFGS',bounds=c(-Inf,Inf),log_link=FALSE,native_error='no prefit'),
 list(name='refinement_bfgs',method='nbinomCR',optimizer='BFGS',bounds=numeric(),log_link=TRUE,hard=TRUE),
 list(name='refinement_lbfgsb',method='nbinomCR',optimizer='L-BFGS-B',bounds=c(-Inf,Inf),log_link=TRUE,hard=TRUE),
 list(name='refinement_lbfgsb_tight',method='nbinomCR',optimizer='L-BFGS-B',bounds=c(-.2,.2),log_link=TRUE,hard=TRUE),
 list(name='refinement_unknown_rejects',method='nbinomCR',optimizer='not_an_optimizer',bounds=numeric(),log_link=TRUE,hard=TRUE,native_error='unsupported refinement optimizer'),
 list(name='refinement_empty_bounds_unbounded',method='nbinomCR',optimizer='L-BFGS-B',bounds=numeric(),log_link=TRUE,hard=TRUE))
cases<-lapply(profiles,function(p){
  y<-Y;design<-x;d<-dispersion;scale<-decode(base$scale)
  if(isTRUE(p$nonzero)){y<-y[nonzero,,drop=FALSE];d<-d[nonzero]}
  if(isTRUE(p$hard)){y<-hard;design<-hardx;d<-rep(.001,nrow(y));scale<-.001}
  pc<-list(no.shrink=1,prior.mean=0,prior.scale=scale,prior.df=1,prior.no.shrink.mean=0,prior.no.shrink.scale=15)
  before<-NULL;after<-NULL
  if(p$method=='nbinomC*'){RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42L);before<-state(.Random.seed)}
  .compiled_policy_capture$calls<-0L;warnings<-character()
  fit<-tryCatch(withCallingHandlers(apeglm::apeglm(Y=y,x=design,param=d,coef=2,threshold=.3,prior.control=pc,
    method=p$method,optim.method=p$optimizer,bounds=p$bounds,log.link=p$log_link),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
  if(!is.null(before))after<-state(.Random.seed)
  if(isTRUE(p$hard))stopifnot(.compiled_policy_capture$calls>0L)
  if(isTRUE(p$nonzero))stopifnot(.compiled_policy_capture$calls==0L)
  list(name=p$name,observations=descriptor(y),design=descriptor(design),parameters=descriptor(matrix(d,ncol=1)),
    offsets=NULL,weights=NULL,mle=NULL,coefficient=1L,nuisance=list(0L),scale=bits(scale)[[1L]],nuisance_scale=bits(15)[[1L]],
    prior_df=bits(1)[[1L]],prior_mean=bits(0)[[1L]],nuisance_mean=bits(0)[[1L]],method=p$method,optimizer=p$optimizer,
    bounds=bits(p$bounds),log_link=p$log_link,threshold=bits(.3)[[1L]],likelihood=NULL,rng_state=before,expected_random_state=after,
    expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,native_expected_error=p$native_error,
    mode=if(inherits(fit,'error'))NULL else expected(fit$map),sd=if(inherits(fit,'error'))NULL else expected(fit$sd),
    fsr=if(inherits(fit,'error'))NULL else expected(fit$fsr),sign_svalues=if(inherits(fit,'error'))NULL else expected(fit$svalue),
    fsos=if(inherits(fit,'error'))NULL else expected(fit$thresh),threshold_svalues=if(inherits(fit,'error'))NULL else expected(apeglm::svalue(fit$thresh)),
    fitted_scale=if(inherits(fit,'error'))NULL else expected(fit$prior.control$prior.scale)[[1L]],prior_variance=NULL,
    original_refinement_calls=.compiled_policy_capture$calls,warnings=as.list(warnings))
})
suppressMessages(untrace('optimNbinom',where=asNamespace('apeglm')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm'))),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
