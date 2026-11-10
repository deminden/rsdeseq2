#!/usr/bin/env Rscript
# One paired dataset fit; exercise the changed public method and expression path.
args <- commandArgs(TRUE)
if(!length(args)%in%c(2L,3L))stop('usage: R_LIBRARY FRESH_OUTPUT_DIR [PAIRED_FITS_RDS]')
library(rsdeseq2,lib.loc=normalizePath(args[[1L]]))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
out <- args[[2L]]
if(dir.exists(out))stop('Use a fresh output directory')
dir.create(out,recursive=TRUE)
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
requireNamespace('DESeq2');requireNamespace('ashr')
if(length(args)==3L) {
  paired <- readRDS(args[[3L]])
  rust<-paired$rust;reference<-paired$reference;settings<-paired$settings;genes<-paired$genes
} else {
  fixture <- .rsdeseq2_large_transform_fixture()
  new_dds <- function()DESeq2::DESeqDataSetFromMatrix(fixture$counts,fixture$colData,~condition)
  rust <- rsdeseq2::DESeq(new_dds(),fitType='mean')
  original <- DESeq2::DESeq(new_dds(),fitType='mean',quiet=TRUE)
  coefficient <- tail(DESeq2::resultsNames(original),1L)
  settings <- list(default=list(),sign=list(svalue=TRUE),threshold=list(lfcThreshold=.5))
  reference <- lapply(settings,function(options)do.call(DESeq2::lfcShrink,c(list(dds=original,coef=coefficient,type='ashr'),options)))
  genes <- nrow(fixture$counts)
  saveRDS(list(rust=rust,reference=reference,settings=settings,genes=genes),file.path(out,'paired_fits.rds'))
}
blocked <- function(expr) {
  # R's dataframe/S4 presentation uses order/sum internally for labels and sizes.
  # The direct native conformance scripts block those as well at the boundary.
  functions <- list(base=c('sqrt','exp','expm1','log','log1p','mean','rank','cumsum','^'),
    stats=c('pnorm','dnorm','pt','dt','integrate','quantile','lowess'),ashr=c('ash','cdf_post','set_data'))
  on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
  for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
    tracer=substitute(stop(MESSAGE,call.=FALSE),list(MESSAGE=paste('native ashr executed an R calculation:',package,name))))))
  force(expr)
}
guard <- tryCatch(blocked(ashr::ash(1,.5)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
actual <- lapply(settings,function(options)blocked(do.call(rsdeseq2::lfcShrink,
  c(list(dds=rust,coef='conditionB',type='ashr',independentFiltering=FALSE),options))))
compare <- function(a,b,bound) {
  stopifnot(length(a)==length(b),identical(is.na(a),is.na(b)),!any(is.infinite(a)),!any(is.infinite(b)))
  finite <- is.finite(a)&is.finite(b)
  error <- max(abs(a[finite]-b[finite]))
  stopifnot(error <= bound)
  error
}
errors <- lapply(names(settings),function(name) {
  a<-actual[[name]];b<-reference[[name]]
  stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'),identical(DESeq2::priorInfo(a)$package,'rsdeseq2'))
  list(profile=name,mean=compare(a$log2FoldChange,b$log2FoldChange,1e-3),
    sd=compare(a$lfcSE,b$lfcSE,1e-3),
    svalue=if(name=='default')NULL else compare(a$svalue,b$svalue,.02))
})
# The optional R functions supply syntax and stored constants only.
scale <- .7;shift <- .2
custom <- list(name='custom_logistic',const=TRUE,
  lcdfFUN=function(x)-log1p(exp(-(x-shift)/scale)),
  lpdfFUN=function(x)-(x-shift)/scale-2*log1p(exp(-(x-shift)/scale))-log(scale))
x<-c(-2,-.3,0,.3,2);se<-c(.2,.7,1,.7,.2)
prior<-ashr::unimix(c(.2,.3,.5),c(0,-1,-4),c(0,1,4))
custom_errors <- lapply(list(custom=custom,student=ashr::lik_t(5)),function(lik){
  arguments<-list(betahat=x,sebetahat=se,mixcompdist='uniform',method='shrink',g=prior,fixg=TRUE,lik=lik)
  expected<-do.call(ashr::ash,arguments)
  result<-blocked(rsdeseq2:::.rsdeseq2_native_ash(arguments,x,se,.4,TRUE))
  data<-ashr::set_data(x,se)
  rates<-ifelse(expected$result$PosteriorMean>0,ashr::cdf_post(expected$fitted_g,.4,data),1-ashr::cdf_post(expected$fitted_g,-.4,data))
  sv<-cumsum(sort(rates))/seq_along(rates);sv<-sv[rank(rates,ties.method='first')]
  list(mean=compare(result$result$PosteriorMean,expected$result$PosteriorMean,2e-6),
    sd=if(all(is.na(expected$result$PosteriorSD))){stopifnot(all(is.na(result$result$PosteriorSD)));NULL}else compare(result$result$PosteriorSD,expected$result$PosteriorSD,2e-6),
    threshold_svalue=compare(result$result$svalue,sv,2e-6))
})
record<-list(R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2')),
  ashr=as.character(packageVersion('ashr')),genes=genes,
  public_errors=errors,likelihood_errors=custom_errors,R_calculations_blocked=TRUE,
  existing_public_bounds_unchanged=TRUE)
jsonlite::write_json(record,file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
cat('Public ashr default/sign/threshold and restricted/Student likelihoods pass; R calculations blocked\n')
