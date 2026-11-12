#!/usr/bin/env Rscript
# Direct native/R-free proof plus one paired public workflow. No full suite.
args<-commandArgs(TRUE)
if(!length(args)%in%c(3L,4L))stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR [PAIRED_FITS_RDS]')
driver<-normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
requireNamespace('DESeq2');requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
decode<-function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')},numeric(1L))
encode<-function(v)lapply(as.numeric(v),hex)
encode_optional<-function(v)lapply(as.numeric(v),function(x)if(is.na(x))NULL else hex(x))
matrix_input<-function(m)if(is.null(m))NULL else matrix(decode(m$values),nrow=m$rows,ncol=m$columns)
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
matrix_bits<-function(m)if(is.null(m))NULL else list(rows=nrow(m),columns=ncol(m),values=bits(m))
compare<-function(a,b,bound){
  stopifnot(length(a)==length(b),identical(unname(is.na(a)),unname(is.na(b))),!any(is.infinite(a)),!any(is.infinite(b)))
  finite<-is.finite(a)&is.finite(b);error<-max(abs(a[finite]-b[finite]));stopifnot(error<=bound);error
}
blocked<-function(expr){
  functions<-list(base=c('sqrt','exp','expm1','log','log1p','log2','mean','rank','cumsum','^'),
    stats=c('pnorm','dnorm','dnbinom','dt','integrate','quantile','optim','optimHess','uniroot'),
    apeglm=c('apeglm','priorVar','svalue'))
  on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
  for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
    tracer=substitute(stop(MESSAGE,call.=FALSE),list(MESSAGE=paste('native apeglm executed an R calculation:',package,name))))))
  force(expr)
}
guard<-tryCatch(blocked(apeglm::svalue(.2)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/cauchy_posterior_r462.json'),simplifyVector=FALSE)
cases<-reference$cases
factor_case<-Filter(function(c)identical(c$name,'offset'),cases)[[1L]]
factor_case$name<-'normalization_factors'
factor_case$normalization_factors<-matrix_bits(exp(matrix_input(factor_case$offsets)))
factor_case$offsets<-NULL
cases<-c(cases,list(factor_case))
input<-file.path(out,'input.json')
input_fields<-c('name','counts','design','dispersions','weights','offsets','mle','coefficient','scale','nuisance','nuisance_scale',
  'threshold','flip_sign','map_only','size_factors','normalization_factors','mle_scale','threshold_scale')
jsonlite::write_json(list(heavy_tail_cases=lapply(cases,function(c)c[intersect(names(c),input_fields)]),
  prior_variance_cases=lapply(reference$priors,function(c)c[c('name','estimates','standard_errors')])),input,auto_unbox=TRUE,null='null')
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
factor_errors<-NULL
wrapped<-lapply(seq_along(cases),function(i){
  c<-cases[[i]]
  pc<-list(no.shrink=unlist(c$nuisance)+1L,prior.mean=0,prior.scale=decode(c$scale),prior.df=1,
    prior.no.shrink.mean=0,prior.no.shrink.scale=decode(c$nuisance_scale))
  arguments<-list(Y=matrix_input(c$counts),x=matrix_input(c$design),param=decode(c$dispersions),
    coef=c$coefficient+1L,weights=matrix_input(c$weights),offset=matrix_input(c$offsets),
    mle=matrix_input(c$mle),prior.control=pc,threshold=decode(c$threshold),flip.sign=c$flip_sign,
    method=if(c$map_only)'nbinomC'else'nbinomCR')
  if(!is.null(c$normalization_factors))arguments$offset<-list(factors=matrix_input(c$normalization_factors))
  r<-blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='apeglm_shrinkage',arguments=arguments)))
  if(identical(c$name,'normalization_factors')){
    original_args<-arguments;original_args$offset<-log(matrix_input(c$normalization_factors))
    original<-do.call(apeglm::apeglm,original_args);coefficient<-c$coefficient+1L
    factor_errors<<-list(mode=compare(r$posterior_mode,original$map[,coefficient]/log(2),.01),
      sd=compare(r$posterior_sd,original$sd[,coefficient]/log(2),.01),
      svalue=compare(r$svalue,original$svalue[,1L],.02),
      threshold_svalue=compare(r$threshold_svalue,apeglm::svalue(original$thresh),.05))
  }
  s<-standalone$cases[[i]]$value
  expected<-list(mode=s$summary_mode,sd=s$summary_sd,svalue=s$sign_svalues,
    threshold_svalue=s$threshold_svalues,prior_scale=s$prior_scale,prior_variance=s$prior_variance,
    nuisance=lapply(s$nuisance,function(v)v+1L),convergence=lapply(s$diagnostics,`[[`,'convergence'))
  actual<-list(mode=encode(r$posterior_mode),sd=encode_optional(r$posterior_sd),svalue=encode_optional(r$svalue),
    threshold_svalue=encode_optional(r$threshold_svalue),prior_scale=encode(r$prior_control$prior.scale),
    prior_variance=if(is.null(r$prior_control$prior.var))NULL else encode(r$prior_control$prior.var),
    nuisance=as.list(as.integer(r$prior_control$no.shrink)),
    convergence=lapply(r$convergence,function(v)if(is.na(v))NULL else as.integer(v)))
  stopifnot(identical(.formula_canonical(actual),.formula_canonical(expected)))
  list(name=c$name,value=actual)
})
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
writeLines(output,file.path(out,'standalone.json'))
# Reuse count fitting when supplied; independent reference shrinkage remains R.
cached<-if(length(args)==4L)readRDS(args[[4L]])else NULL
fixture<-.rsdeseq2_large_transform_fixture()
new_dds<-function()DESeq2::DESeqDataSetFromMatrix(fixture$counts,fixture$colData,~condition)
rust<-if(is.null(cached))rsdeseq2::DESeq(new_dds(),fitType='mean')else cached$rust
stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(rust,'counts')),fixture$counts,check.attributes=TRUE,tolerance=0)))
settings<-list(default=list(),sign=list(svalue=TRUE),threshold=list(lfcThreshold=.5))
if(!is.null(cached)&&identical(cached$reference_type,'apeglm')) {
  original<-cached$reference
}else{
  reference_fit<-DESeq2::DESeq(new_dds(),fitType='mean',quiet=TRUE)
  coefficient<-tail(DESeq2::resultsNames(reference_fit),1L)
  original<-lapply(settings,function(options)do.call(DESeq2::lfcShrink,
    c(list(dds=reference_fit,coef=coefficient,type='apeglm',quiet=TRUE),options)))
}
saveRDS(list(rust=rust,reference=original,reference_type='apeglm'),file.path(out,'paired_fits.rds'))
actual<-lapply(settings,function(options)blocked(do.call(rsdeseq2::lfcShrink,
  c(list(dds=rust,coef='conditionB',type='apeglm',quiet=TRUE),options))))
# Check the public workflow's native size-factor/log2 input path directly
# against a separate R-free process, beyond the supplied-natural-input panel.
fit<-rsdeseq2:::.rsdeseq2_fitted_view(rust)
raw<-rsdeseq2::results(rust,name='conditionB')
public_case<-list(name='public_normalization_log2',counts=matrix_bits(fit$counts),design=matrix_bits(fit$design),
  dispersions=bits(fit$dispersions),offsets=NULL,size_factors=if(is.null(fit$normalizationFactors))bits(fit$sizeFactors)else NULL,
  normalization_factors=matrix_bits(fit$normalizationFactors),weights=matrix_bits(fit$weights),
  mle=matrix_bits(cbind(raw$log2FoldChange,raw$lfcSE)),mle_scale='log2',threshold_scale='log2',
  coefficient=1L,scale=bits(1)[[1L]],nuisance=list(0L),nuisance_scale=bits(15)[[1L]],
  threshold=bits(.5)[[1L]],flip_sign=TRUE,map_only=FALSE)
public_input<-file.path(out,'public-input.json')
jsonlite::write_json(list(heavy_tail_cases=list(public_case)),public_input,auto_unbox=TRUE,null='null')
public_output<-system2(driver,shQuote(normalizePath(public_input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(public_output,'status')))stop(paste(public_output,collapse='\n'))
writeLines(public_output,file.path(out,'public-standalone.json'))
public_native<-jsonlite::fromJSON(paste(public_output,collapse='\n'),simplifyVector=FALSE)$cases[[1L]]$value
pc<-DESeq2::priorInfo(actual$threshold)$prior.control
stopifnot(identical(encode(actual$threshold$log2FoldChange),public_native$summary_mode),
  identical(encode_optional(actual$threshold$lfcSE),public_native$summary_sd),
  identical(encode_optional(actual$sign$svalue),public_native$sign_svalues),
  identical(encode_optional(actual$threshold$svalue),public_native$threshold_svalues),
  identical(encode(pc$prior.scale),public_native$prior_scale),
  identical(encode(pc$prior.var),public_native$prior_variance))
errors<-lapply(names(settings),function(name){
  a<-actual[[name]];b<-original[[name]]
  stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'),identical(DESeq2::priorInfo(a)$package,'rsdeseq2'))
  list(profile=name,mode=compare(a$log2FoldChange,b$log2FoldChange,.01),sd=compare(a$lfcSE,b$lfcSE,.01),
    svalue=if(name=='default')NULL else compare(a$svalue,b$svalue,if(name=='sign').02 else .05))
})
jsonlite::write_json(list(R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2')),
  apeglm=as.character(packageVersion('apeglm')),genes=nrow(fixture$counts),
  complete_native_standalone_batches=length(cases),public_normalization_log2_bits_identical=TRUE,
  public_errors=errors,normalization_factor_errors=factor_errors,R_calculations_blocked=TRUE,
  existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
cat(length(cases),'native/standalone batches: computed columns and fitted-prior values match bit for bit; public apeglm default/sign/threshold pass with R calculations blocked\n')
