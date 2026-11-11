#!/usr/bin/env Rscript
# Native/wrapped kernel identity and public workflows, never R-backed fitting.
args<-commandArgs(TRUE)
if(length(args)!=4L)stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS')
driver<-normalizePath(args[[1L]]);library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
requireNamespace('ashr');requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/native_posterior_check_helpers.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/continuous_likelihood_r462.json'),simplifyVector=FALSE)
input<-file.path(out,'input.json');jsonlite::write_json(list(builtin_likelihood_cases=reference$kernel_cases),input,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE);if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'));standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
canonical<-function(v){if(!is.list(v))return(v);if(!is.null(names(v)))v<-v[order(names(v))];lapply(v,canonical)}
wrapped<-lapply(reference$kernel_cases,function(c){
  m<-c$builtin_likelihood;model<-if(m$type=='normalmix')list(type='normalmix',weights=decode(m$weights),scales=.formula_decode(m$scales))else list(type='logF',df1=decode(m$df1),df2=decode(m$df2))
  r<-blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=model,arguments=decode(c$arguments),lower=decode(c$lower),upper=decode(c$upper))))
  list(name=c$name,ok=TRUE,value=lapply(r,encode))
})
wrapped<-.formula_canonical(list(cases=wrapped));jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
stopifnot(identical(canonical(wrapped),canonical(standalone)))
cached<-readRDS(args[[4L]]);dds<-cached$rust[seq_len(128L),];raw<-rsdeseq2::results(dds,name='conditionB')
profiles<-list(normalmix=ashr::lik_normalmix(c(.7,.3),c(.5,2)),normalmix_rows=ashr::lik_normalmix(c(.7,.3),cbind(seq(.3,1.2,length.out=128),seq(1.3,3,length.out=128))),
 logF=ashr::lik_logF(5,8),logF_rows=ashr::lik_logF(seq(2,9,length.out=128),seq(10,3,length.out=128)),
 logF_first_infinite=ashr::lik_logF(Inf,8),logF_second_infinite=ashr::lik_logF(5,Inf),
 threshold=ashr::lik_normalmix(c(.7,.3),c(.5,2)),estimated_mode=ashr::lik_normalmix(c(.7,.3),c(.5,2)))
public_errors<-lapply(names(profiles),function(name){
  options<-list(mixcompdist='uniform',lik=profiles[[name]],mixsd=c(.2,.7,2));if(name=='estimated_mode')options$mode<-'estimate'
  original<-tryCatch(do.call(ashr::ash,c(list(betahat=raw$log2FoldChange,sebetahat=raw$lfcSE,method='shrink'),options)),error=identity)
  threshold<-if(name=='threshold').25 else 0
  actual<-tryCatch(blocked(do.call(rsdeseq2::lfcShrink,c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,lfcThreshold=threshold,quiet=TRUE),options))),error=identity)
  saveRDS(list(original=original,actual=actual,options=options,estimates=raw$log2FoldChange,standard_errors=raw$lfcSE),file.path(out,paste0('public-',name,'.rds')))
  if(inherits(original,'error')){
    stopifnot(inherits(actual,'error'))
    message(name,': both interfaces reject the original invalid likelihood matrix')
    return(list(profile=name,rejected=TRUE,original_error=conditionMessage(original),native_error=conditionMessage(actual)))
  }
  if(inherits(actual,'error'))stop(actual)
  sv<-original$result$svalue
  # DESeq2's positive-threshold queries deliberately use set_data's normal default.
  if(threshold>0){a<-ashr::cdf_post(original$fitted_g,threshold,ashr::set_data(raw$log2FoldChange,raw$lfcSE));b<-ashr::cdf_post(original$fitted_g,-threshold,ashr::set_data(raw$log2FoldChange,raw$lfcSE));rates<-ifelse(original$result$PosteriorMean>0,a,1-b);sv<-DESeq2:::apeglm_svalue(rates)}
  stopifnot(identical(DESeq2::priorInfo(actual)$engine,'Rust'))
  message(name,': mean error ',max(abs(actual$log2FoldChange-original$result$PosteriorMean),na.rm=TRUE))
  list(profile=name,mean=compare(actual$log2FoldChange,original$result$PosteriorMean,1e-3),sd=compare(actual$lfcSE,original$result$PosteriorSD,1e-3),svalue=compare(actual$svalue,sv,.02))
})
jsonlite::write_json(list(complete_kernel_batches=length(wrapped$cases),R_calculations_blocked=TRUE,public_errors=public_errors,existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat(length(wrapped$cases),'built-in kernel batches match native/standalone bits; eight public profiles checked\n')
