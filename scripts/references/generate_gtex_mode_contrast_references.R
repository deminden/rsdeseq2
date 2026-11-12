#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: Rscript generate_references.R CONTRAST_DIR")
d <- normalizePath(args[[1]], mustWork=TRUE)
suppressPackageStartupMessages({library(DESeq2); library(S4Vectors)})
stopifnot(as.character(packageVersion("DESeq2")) == "1.52.0")
ct <- read.delim(file.path(d,"counts.tsv"),check.names=FALSE)
ds <- read.delim(file.path(d,"design_matrix.tsv"),check.names=FALSE)
cm <- as.matrix(ct[,-1,drop=FALSE]); storage.mode(cm)<-"integer"; rownames(cm)<-ct[[1]]
dm <- as.matrix(ds[,-1,drop=FALSE]); storage.mode(dm)<-"double"; rownames(dm)<-ds[[1]]; colnames(dm)<-names(ds)[-1]
dds <- DESeqDataSetFromMatrix(cm,DataFrame(row.names=ds[[1]]),dm)
dds <- DESeq(dds,fitType="parametric",betaPrior=FALSE,minReplicatesForReplace=Inf,quiet=TRUE)
coef <- resultsNames(dds)[length(resultsNames(dds))]
r <- results(dds,name=coef,cooksCutoff=FALSE,independentFiltering=FALSE)
write.table(data.frame(gene=rownames(r),dispersion=dispersions(dds),baseMean=r$baseMean,
  log2FoldChange=r$log2FoldChange,lfcSE=r$lfcSE,stat=r$stat,pvalue=r$pvalue,padj=r$padj),
  file.path(d,"reference_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
x <- as.data.frame(mcols(dds)); x$gene <- rownames(dds)
keep <- intersect(c("gene","dispGeneEst","dispFit","dispMAP","dispersion","dispOutlier","dispIter","betaConv","betaIter"),names(x))
write.table(x[,keep,drop=FALSE],file.path(d,"reference_fit_diagnostics.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(data.frame(sample=colnames(dds),sizeFactor=sizeFactors(dds)),
  file.path(d,"reference_size_factors.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(d,"sessionInfo.txt"))
cat(sprintf("Wrote DESeq2 reference for %d genes x %d samples (%s)\n",nrow(r),ncol(cm),coef))
