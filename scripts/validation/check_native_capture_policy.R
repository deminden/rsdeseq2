#!/usr/bin/env Rscript
# Public canonical import and compiled controls; reuse count fitting, block R math.
args<-commandArgs(TRUE)
if(length(args)!=3L)stop('usage: R_LIBRARY FRESH_OUTPUT_DIR CACHED_PAIRED_FITS_RDS')
library(rsdeseq2,lib.loc=normalizePath(args[[1L]]));requireNamespace('ashr');requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/validation/native_posterior_check_helpers.R'))
out<-args[[2L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
cached<-readRDS(args[[3L]]);dds<-cached$rust[seq_len(128L),];raw<-rsdeseq2::results(dds,name='conditionB')
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/normalmix_capture_r462.json'),simplifyVector=FALSE)
selected<-Filter(function(c)c$name%in%c('matrix_equal','matrix_transposed_dims','matrix_weights_vector_sd','vector_weights_flat_sd'),reference$kernel_cases)
matrix_errors<-lapply(selected,function(c){
  m<-c$builtin_likelihood;w<-if(!is.null(m$weights$type)).formula_decode(m$weights)else decode(m$weights)
  lik<-ashr::lik_normalmix(w,.formula_decode(m$scales))
  options<-list(lik=lik,mixcompdist='uniform',mixsd=c(.2,.7,2))
  original<-do.call(ashr::ash,c(list(betahat=raw$log2FoldChange,sebetahat=raw$lfcSE,method='shrink'),options))
  actual<-blocked(do.call(rsdeseq2::lfcShrink,c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,quiet=TRUE),options)))
  stopifnot(identical(DESeq2::priorInfo(actual)$engine,'Rust'))
  saveRDS(list(original=original,actual=actual,options=options),file.path(out,paste0('matrix-',c$name,'.rds')))
  list(profile=c$name,mean=compare(actual$log2FoldChange,original$result$PosteriorMean,1e-3),
    sd=compare(actual$lfcSE,original$result$PosteriorSD,1e-3),svalue=compare(actual$svalue,original$result$svalue,.02))
})
nonzero<-rowSums(SummarizedExperiment::assay(dds,'counts'))>0
positive<-dds[nonzero,]
profiles<-list(map_bounded=list(method='nbinomC',optimizer='L-BFGS-B',bounds=c(-.01,.01),log_link=FALSE),
 map_unknown=list(method='nbinomC',optimizer='not_an_optimizer',bounds=numeric(),log_link=FALSE),
 random_unknown=list(method='nbinomC*',optimizer='not_an_optimizer',bounds=numeric(),log_link=FALSE),
 cr_bounded=list(method='nbinomCR',optimizer='L-BFGS-B',bounds=c(-.01,.01),log_link=FALSE,positive=TRUE))
compiled_errors<-lapply(names(profiles),function(name){
  p<-profiles[[name]];data<-if(isTRUE(p$positive))positive else dds
  fit<-rsdeseq2:::.rsdeseq2_fitted_view(data)
  y<-as.matrix(SummarizedExperiment::assay(data,'counts'))
  offset<-if(is.null(fit$normalizationFactors))matrix(log(fit$sizeFactors),nrow(y),ncol(y),byrow=TRUE)else log(fit$normalizationFactors)
  if(endsWith(p$method,'*')){RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42L);before<-.Random.seed}
  original<-tryCatch(apeglm::apeglm(Y=y,x=fit$design,param=fit$dispersions,offset=offset,weights=fit$weights,
    coef=2,method=p$method,optim.method=p$optimizer,bounds=p$bounds,log.link=p$log_link),error=identity)
  if(endsWith(p$method,'*')){after<-.Random.seed;assign('.Random.seed',before,.GlobalEnv)}
  actual<-tryCatch(blocked(rsdeseq2::lfcShrink(data,coef='conditionB',type='apeglm',apeMethod=p$method,apeAdapt=FALSE,
    svalue=TRUE,optim.method=p$optimizer,bounds=p$bounds,log.link=p$log_link,quiet=TRUE)),error=identity)
  saveRDS(list(original=original,actual=actual,options=p),file.path(out,paste0('compiled-',name,'.rds')))
  if(inherits(original,'error'))stop(original)
  if(inherits(actual,'error'))stop(actual)
  if(endsWith(p$method,'*'))stopifnot(identical(.Random.seed,after))
  stopifnot(identical(DESeq2::priorInfo(actual)$engine,'Rust'))
  list(profile=name,mode=compare(actual$log2FoldChange,original$map[,2]/log(2),.01),
    sd=compare(actual$lfcSE,original$sd[,2]/log(2),.01),svalue=compare(actual$svalue,as.numeric(original$svalue),.02))
})
jsonlite::write_json(list(matrix_profiles=matrix_errors,compiled_profiles=compiled_errors,R_calculations_blocked=TRUE,
  count_fits_reused=TRUE,existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat('Four public matrix-capture and four compiled-policy profiles checked; R calculations blocked\n')
