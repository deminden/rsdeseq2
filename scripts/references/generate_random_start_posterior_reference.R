#!/usr/bin/env Rscript
# Original apeglm only. Reuse frozen INPUT fields, never candidate outputs.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('apeglm'))=='1.34.0',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
hex<-function(v,size)paste(format(writeBin(v,raw(),size=size,endian='big')),collapse='')
state<-function(v)lapply(v,function(x)hex(as.integer(x),4L))
decode<-function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1,15,2),seq(2,16,2)),16L)),'double',n=1,size=8,endian='big')},numeric(1))
mat<-function(m)if(is.null(m))NULL else matrix(decode(m$values),nrow=m$rows,ncol=m$columns)
expected<-function(v)lapply(as.numeric(v),function(x)if(is.na(x)&&!is.nan(x))'NA'else sprintf('%.17g',x))
nb<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/cauchy_posterior_r462.json'),simplifyVector=FALSE)$cases
bb<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/beta_binomial_posterior_r462.json'),simplifyVector=FALSE)$cases
bb<-Filter(function(c)identical(c$method,'betabinC'),bb)
# A wholly zero NB batch must still consume its one batch start vector.
zero<-nb[[1L]];zero$name<-'all_zero_batch';zero$counts$values<-lapply(zero$counts$values,function(x)'bits:0000000000000000');zero$mle<-NULL
nb<-c(nb,list(zero))
nb_fields<-c('name','counts','design','dispersions','offsets','weights','mle','coefficient','scale','nuisance','nuisance_scale','threshold','flip_sign','map_only')
bb_fields<-c('name','observations','design','parameters','offsets','weights','mle','coefficient','nuisance','scale','nuisance_scale',
 'prior_df','prior_mean','nuisance_mean','method','optimizer','bounds','log_link','threshold','likelihood')
record<-function(c,beta,kind,seed,explicit_state=NULL){
 method<-if(beta)'betabinC*'else'nbinomC*'
 input<-c[if(beta)bb_fields else nb_fields]
 input$name<-paste(kind,seed,c$name,sep='_')
 if(beta)input$method<-method else input$map_only<-TRUE
 RNGkind(kind,'Inversion','Rejection');set.seed(seed)
 if(!is.null(explicit_state))assign('.Random.seed',explicit_state,.GlobalEnv)
 before<-.Random.seed
 initial<-rnorm(c$design$columns,0,.5);expected_after<-.Random.seed
 assign('.Random.seed',before,.GlobalEnv)
 input$rng_key<-if(is.null(explicit_state))paste(kind,seed,sep='_')else'MT_integer_NA'
 if(!is.null(explicit_state))input$name<-paste0('integer_NA_',input$name)
 input$rng_state<-state(before)
 pc<-list(no.shrink=as.numeric(unlist(c$nuisance))+1,prior.mean=0,prior.scale=decode(c$scale),prior.df=1,
  prior.no.shrink.mean=0,prior.no.shrink.scale=decode(c$nuisance_scale))
 warnings<-character()
 fit<-withCallingHandlers(apeglm::apeglm(Y=mat(if(beta)c$observations else c$counts),x=mat(c$design),
  param=if(beta)mat(c$parameters)else decode(c$dispersions),coef=c$coefficient+1L,mle=mat(c$mle),
  weights=mat(c$weights),offset=mat(c$offsets),threshold=decode(c$threshold),prior.control=pc,
  flip.sign=if(beta)TRUE else c$flip_sign,method=method),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(identical(.Random.seed,expected_after))
 c(input,list(mode=expected(fit$map),sd=expected(fit$sd),fsr=expected(fit$fsr),sign_svalues=expected(fit$svalue),
  fsos=expected(fit$thresh),threshold_svalues=expected(apeglm::svalue(fit$thresh)),
  fitted_scale=expected(fit$prior.control$prior.scale)[[1L]],
  prior_variance=if(is.null(fit$prior.control$prior.var))NULL else expected(fit$prior.control$prior.var)[[1L]],
  expected_start=lapply(initial,function(v)hex(as.double(v),8L)),expected_random_state=state(.Random.seed),warnings=as.list(warnings)))
}
ncases<-list();bcases<-list()
for(kind in c('Mersenne-Twister',"L'Ecuyer-CMRG"))for(seed in c(1L,42L,781L)){
 ncases<-c(ncases,lapply(nb,record,beta=FALSE,kind=kind,seed=seed))
 bcases<-c(bcases,lapply(bb,record,beta=TRUE,kind=kind,seed=seed))
}
# Integer NA is a legitimate MT word bit pattern, including through the R wire.
edge<-c(10403L,1L,1L,0L,NA_integer_,rep(1L,621L))
ncases<-c(ncases,list(record(nb[[1L]],FALSE,'Mersenne-Twister',0L,edge)))
bcases<-c(bcases,list(record(bb[[1L]],TRUE,'Mersenne-Twister',0L,edge)))
# Share identical state vectors across profiles; preserve every source bit.
rng_states<-list()
pack<-function(cases)lapply(cases,function(c){
 key<-c$rng_key;p<-as.character(c$design$columns)
 if(is.null(rng_states[[key]]))rng_states[[key]]<<-list(before=c$rng_state,after=list(),starts=list())
 stopifnot(identical(rng_states[[key]]$before,c$rng_state))
 if(!is.null(rng_states[[key]]$after[[p]]))stopifnot(identical(rng_states[[key]]$after[[p]],c$expected_random_state))
 rng_states[[key]]$after[[p]]<<-c$expected_random_state
 rng_states[[key]]$starts[[p]]<<-c$expected_start
 c$rng_state<-NULL;c$expected_random_state<-NULL;c$expected_start<-NULL;c
})
ncases<-pack(ncases);bcases<-pack(bcases)
jsonlite::write_json(list(source=list(R=as.character(getRversion()),apeglm=as.character(packageVersion('apeglm'))),
 heavy_tail_cases=ncases,regression_cases=bcases,rng_states=rng_states),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
