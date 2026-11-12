#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: FIXTURE_DIR")
suppressPackageStartupMessages(library(DESeq2))
stopifnot(getRversion() == "4.6.2", packageVersion("DESeq2") == "1.52.0")
root <- args[[1L]]
ct <- read.delim(file.path(root,"counts.tsv"), check.names=FALSE)
dm <- read.delim(file.path(root,"design.tsv"), check.names=FALSE)
counts <- as.matrix(ct[,-1,drop=FALSE]);storage.mode(counts)<-"integer";rownames(counts)<-ct[[1L]]
design <- as.matrix(dm[,-1,drop=FALSE]);storage.mode(design)<-"double";rownames(design)<-dm[[1L]]
stopifnot(identical(colnames(counts),rownames(design)))
fit <- DESeq(DESeqDataSetFromMatrix(counts,DataFrame(row.names=colnames(counts)),design),
             fitType="parametric",betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE,quiet=TRUE)
if (attr(dispersionFunction(fit),"fitType") != "parametric") stop("Reference changed trend family")
writeResult <- function(object, name) {
 out <- results(object,name=tail(resultsNames(object),1L),cooksCutoff=FALSE,independentFiltering=FALSE)
 write.table(data.frame(gene=rownames(out),dispersion=dispersions(object),baseMean=out$baseMean,
   log2FoldChange=out$log2FoldChange,lfcSE=out$lfcSE,stat=out$stat,pvalue=out$pvalue,padj=out$padj),
   file.path(root,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
writeResult(fit,"reference_wald.tsv")
lrt <- nbinomLRT(fit,full=design,reduced=design[,1:2,drop=FALSE],quiet=TRUE)
writeResult(lrt,"reference_lrt.tsv")
writeLines(capture.output(sessionInfo()),file.path(root,"R-session.txt"))
