#!/usr/bin/env Rscript
# Original DESeq2 replacement and metadata routines only; no candidate outputs.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',as.character(packageVersion('DESeq2'))=='1.52.0',!file.exists(args[[1L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
expected<-function(v)lapply(as.numeric(v),function(x)sprintf('%.17g',x))
mat<-function(m)if(is.null(m))NULL else list(rows=nrow(m),columns=ncol(m),values=bits(m))
counts<-.rsdeseq2_test_counts();counts[1,8]<-2000;counts[2,]<-0;counts[2,8]<-1000;counts[40,]<-0
x<-cbind(Intercept=1,conditionB=rep(0:1,each=4));sf<-exp(seq(log(.8),log(1.2),length.out=8))
weights<-outer(seq(.2,1.4,length.out=40),c(.3,1,.7,1,1,.4,1,.8))
nf<-outer(seq(.9,1.1,length.out=40),sf)
cooks<-matrix(0,40,8);cooks[1:2,8]<-1000
profiles<-list(unweighted=list(),weighted=list(weights=weights),weighted_factors=list(weights=weights,factors=nf),
 weighted_no_replacement=list(weights=weights,cutoff=1e9))
cases<-lapply(names(profiles),function(name){
 p<-profiles[[name]];cutoff<-if(is.null(p$cutoff))1 else p$cutoff
 dds<-DESeq2::DESeqDataSetFromMatrix(counts,data.frame(condition=factor(rep(c('A','B'),each=4)),row.names=colnames(counts)),~condition)
 DESeq2::sizeFactors(dds)<-sf
 if(!is.null(p$factors))DESeq2::normalizationFactors(dds)<-p$factors
 if(!is.null(p$weights))SummarizedExperiment::assay(dds,'weights',withDimnames=FALSE)<-p$weights
 SummarizedExperiment::assay(dds,'cooks',withDimnames=FALSE)<-cooks;attr(dds,'modelMatrix')<-x
 dds<-DESeq2:::getBaseMeansAndVariances(dds)
 replaced<-DESeq2:::replaceOutliers(dds,trim=.2,cooksCutoff=cutoff,minReplicates=3,
  whichSamples=c(rep(FALSE,7),TRUE))
 nrefit<-sum(S4Vectors::mcols(replaced)$replace,na.rm=TRUE)
 if(nrefit>0L)replaced<-DESeq2:::getBaseMeansAndVariances(replaced)
 metadata<-S4Vectors::mcols(replaced)
 list(name=name,counts=mat(counts),design=mat(x),size_factors=bits(sf),normalization_factors=mat(p$factors),weights=mat(p$weights),
  cooks=mat(cooks),cutoff=bits(cutoff)[[1L]],which_samples=as.list(c(rep(FALSE,7),TRUE)),
  replaced_counts=mat(SummarizedExperiment::assay(replaced,'counts')),base_mean=expected(metadata$baseMean),base_var=expected(metadata$baseVar),
  all_zero=as.list(metadata$allZero),n_refit=nrefit)
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2'))),cases=cases),
 args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
