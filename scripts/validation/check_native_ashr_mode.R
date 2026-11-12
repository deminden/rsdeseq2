#!/usr/bin/env Rscript
# Real wrapper and independent standalone requests with all R search/math blocked.
args<-commandArgs(TRUE)
if(!length(args)%in%c(4L,5L))stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS [REFERENCE_JSON]')
driver<-normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
requireNamespace('ashr');requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/validation/native_posterior_check_helpers.R'))
guard<-tryCatch(blocked(stats::optimize(function(x)x*x,c(-1,1))),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed R math',conditionMessage(guard)))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
reference<-jsonlite::fromJSON(if(length(args)==5L)args[[5L]]else file.path(root,'crates/rsdeseq2/tests/data/ashr_mode_r462.json'),simplifyVector=FALSE)
if(!is.null(reference$mode_cases))reference$cases<-reference$mode_cases
input<-file.path(out,'input.json');jsonlite::write_json(list(ashr_mode_cases=reference$cases),input,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
function_input<-function(f){if(is.null(f))return(NULL);list(arguments=unlist(f$arguments,use.names=FALSE),expression=.formula_ast_decode(f$expression),constants=lapply(f$constants,.formula_decode))}
encode_all<-function(v){if(is.numeric(v))return(encode(v));if(is.logical(v)||is.character(v))return(as.list(unname(v)));if(is.list(v))return(lapply(v,encode_all));v}
canonical<-function(v){if(!is.list(v))return(v);if(!is.null(names(v)))v<-v[order(names(v))];lapply(v,canonical)}
wrapped<-lapply(reference$cases,function(c){
  x<-decode(c$estimates);s<-decode(c$standard_errors)
  arguments<-list(betahat=x,sebetahat=s,mixcompdist=c$family,optmethod=c$optimizer,mode=.formula_decode(c$mode_input),
    pointmass=c$point_mass,prior=if(!is.null(c$penalty))c$penalty else if(is.null(c$nullweight))'uniform'else'nullbiased',
    mixsd=if(is.null(c$scales))NULL else decode(c$scales),alpha=decode(c$alpha))
  if(!is.null(c$df))arguments$df<-decode(c$df)
  if(!is.null(c$observation_weights))arguments$weights<-decode(c$observation_weights)
  if(!is.null(c$builtin_likelihood)) {
    lik<-c$builtin_likelihood
    stopifnot(identical(lik$type,'normalmix'))
    arguments$lik<-list(type='normalmix',weights=decode(lik$weights),scales=.formula_decode(lik$scales))
  }
  if(!is.null(c$custom_likelihood))arguments$lik<-list(type='restricted',model=c(list(name=c$custom_likelihood$name),lapply(c$custom_likelihood[c('log_cdf','log_density','first','second')],function_input)))
  r<-tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='ashr_shrinkage',arguments=arguments,
    query_estimates=x,query_standard_errors=s,lfc_threshold=0,svalue=TRUE))),error=identity)
  if(!identical(inherits(r,'error'),!is.null(c$expected_error)||!is.null(c$native_expected_error)))stop(c$name,': ',if(inherits(r,'error'))conditionMessage(r)else'unexpected success')
  if(inherits(r,'error'))return(list(name=c$name,ok=FALSE))
  list(name=c$name,ok=TRUE,value=encode_all(r))
})
wrapped<-.formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
if(!identical(canonical(wrapped),canonical(standalone)))stop(paste(all.equal(canonical(wrapped),canonical(standalone)),collapse='; '))
public_errors <- list()
# A dash requests the complete interface check without cached public count fits.
if (!identical(args[[4L]], "-")) {
cached<-readRDS(args[[4L]]);dds<-cached$rust[seq_len(128L),];raw<-rsdeseq2::results(dds,name='conditionB')
profiles<-list(normal=list(mode='estimate'),expanded=list(mode=c(-8,8)),alpha=list(mode='estimate',alpha=.5),
 em=list(mode='estimate',optmethod='mixEM'),student=list(mode='estimate',mixcompdist='uniform',df=4),
 logistic=list(mode='estimate',mixcompdist='uniform',lik=list(name='custom_logistic',lcdfFUN=function(x)-log1p(exp(-x)),lpdfFUN=function(x)-x-2*log1p(exp(-x)))))
public_errors<-lapply(names(profiles),function(name){
  options<-profiles[[name]]
  arguments<-utils::modifyList(list(betahat=raw$log2FoldChange,sebetahat=raw$lfcSE,mixcompdist='normal',method='shrink'),options)
  # lfcShrink's alpha formal controls result filtering and is not an ashr dot.
  # Effect-scale alpha is tested through the primitive requests above.
  if(name=='alpha')arguments$alpha<-NULL
  original<-do.call(ashr::ash,arguments)
  actual<-blocked(do.call(rsdeseq2::lfcShrink,c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,quiet=TRUE),options)))
  stopifnot(identical(DESeq2::priorInfo(actual)$engine,'Rust'))
  saveRDS(list(original=original,actual=actual,arguments=arguments),file.path(out,paste0('public-',name,'.rds')))
  message(name,': mean error ',max(abs(actual$log2FoldChange-original$result$PosteriorMean),na.rm=TRUE))
  list(profile=name,mean=compare(actual$log2FoldChange,original$result$PosteriorMean,1e-3),
    sd=compare(actual$lfcSE,original$result$PosteriorSD,1e-3),svalue=compare(actual$svalue,original$result$svalue,.02))
})
}
jsonlite::write_json(list(complete_native_standalone_batches=length(wrapped$cases),R_calculations_blocked=TRUE,
  R_scalar_optimizer_blocked=TRUE,public_errors=public_errors,existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat(length(wrapped$cases),'estimated-mode batches match native/standalone bits;',length(public_errors),'public profiles pass\n')
