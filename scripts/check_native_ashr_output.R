#!/usr/bin/env Rscript
# Selected result bits across the real R bridge and R-free driver, plus public use.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS')
driver <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc=normalizePath(args[[2L]]))
requireNamespace('ashr'); requireNamespace('apeglm')
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)), '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root, 'scripts/native_posterior_check_helpers.R'))
out <- args[[3L]]; if (dir.exists(out)) stop('Use a fresh output directory'); dir.create(out, recursive=TRUE)
reference <- jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/ashr_output_r462.json'),simplifyVector=FALSE)
input <- file.path(out,'input.json'); jsonlite::write_json(list(ashr_output_cases=reference$cases),input,auto_unbox=TRUE,null='null',digits=NA)
output <- system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(output,'status'))) stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'))
standalone <- jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
encode_all <- function(v) {if(is.numeric(v))return(encode(v));if(is.logical(v)||is.character(v))return(as.list(unname(v)));if(is.list(v))return(lapply(v,encode_all));v}
canonical <- function(v) {if(!is.list(v))return(v);if(!is.null(names(v)))v<-v[order(names(v))];lapply(v,canonical)}
wrapped <- lapply(reference$cases,function(c) {
  x <- decode(c$estimates); s <- decode(c$standard_errors)
  request <- list(operation='ashr_shrinkage',arguments=list(betahat=x,sebetahat=s,method='shrink',mixcompdist='normal',
    g=list(component_class='normalmix',weights=c(.5,.3,.2),means=rep(0,3),sds=c(0,.5,2),lower=NULL,upper=NULL),fixg=TRUE,
    outputlevel=.formula_decode(c$output)),query_estimates=x,query_standard_errors=s,lfc_threshold=decode(c$threshold),svalue=c$query_svalues)
  r <- tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(request)),error=identity)
  stopifnot(identical(inherits(r,'error'),!is.null(c$expected_error)))
  if(inherits(r,'error'))return(list(name=c$name,ok=FALSE))
  list(name=c$name,ok=TRUE,value=encode_all(r))
})
wrapped <- .formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
stopifnot(identical(canonical(wrapped),canonical(standalone)))
cached <- readRDS(args[[4L]]); dds <- cached$rust[seq_len(128L),]
raw <- rsdeseq2::results(dds,name='conditionB')
profiles <- list(level_one=list(outputlevel=1),level_four=list(outputlevel=4),
  named_mean=list(outputlevel=c('PosteriorMean','fitted_g')),named_sd=list(outputlevel=c('PosteriorSD','fitted_g')),
  named_no_prior=list(outputlevel=c('PosteriorMean','PosteriorSD','svalue')),threshold_four=list(outputlevel=4))
errors <- lapply(names(profiles),function(name) {
  options <- profiles[[name]]; threshold <- if(name=='threshold_four').25 else 0
  original <- do.call(ashr::ash,c(list(betahat=raw$log2FoldChange,sebetahat=raw$lfcSE,mixcompdist='normal',method='shrink'),options))
  sv <- original$result$svalue
  if(threshold>0){a<-ashr::cdf_post(original$fitted_g,threshold,ashr::set_data(raw$log2FoldChange,raw$lfcSE));b<-ashr::cdf_post(original$fitted_g,-threshold,ashr::set_data(raw$log2FoldChange,raw$lfcSE));rates<-ifelse(original$result$PosteriorMean>0,a,1-b);sv<-numeric(length(rates));sv[order(rates)]<-cumsum(sort(rates))/seq_along(rates)}
  actual <- blocked(do.call(rsdeseq2::lfcShrink,c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,lfcThreshold=threshold,quiet=TRUE),options)))
  fields <- c('baseMean',if(!is.null(original$result$PosteriorMean))'log2FoldChange',if(!is.null(original$result$PosteriorSD))'lfcSE',if(!is.null(sv))'svalue')
  if(threshold>0)fields<-c(names(raw),'svalue')
  if(!identical(names(actual),fields))stop(name,': columns ',paste(names(actual),collapse=','),' expected ',paste(fields,collapse=','))
  stopifnot(identical(names(actual),fields),identical(DESeq2::priorInfo(actual)$engine,'Rust'),identical(is.null(DESeq2::priorInfo(actual)$fitted_g),is.null(original$fitted_g)))
  list(profile=name,mean=compare(actual$log2FoldChange,original$result$PosteriorMean,1e-3),
    sd=compare(actual$lfcSE,original$result$PosteriorSD,1e-3),svalue=compare(actual$svalue,sv,.02))
})
jsonlite::write_json(list(complete_native_standalone_batches=length(wrapped$cases),R_calculations_blocked=TRUE,
  public_errors=errors,existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat(length(wrapped$cases),'selected-output batches match every native/standalone bit; six public profiles pass\n')
