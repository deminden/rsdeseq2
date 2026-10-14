#!/usr/bin/env Rscript
# Reference-side process benchmark; never used by the Rust runtime.
args <- commandArgs(TRUE)
suppressPackageStartupMessages({library(DESeq2); library(S4Vectors)})
stopifnot(as.character(packageVersion("DESeq2")) == "1.52.0")
if (identical(args, "--provenance")) {
  writeLines(capture.output(sessionInfo()))
  quit(save="no")
}
if (length(args) != 3L) stop("usage: Rscript benchmark_gtex_native_wald.R COUNTS DESIGN OUTPUT")
ct <- read.delim(args[[1]],check.names=FALSE)
ds <- read.delim(args[[2]],check.names=FALSE)
cm <- as.matrix(ct[,-1,drop=FALSE]); storage.mode(cm) <- "integer"; rownames(cm) <- ct[[1]]
dm <- as.matrix(ds[,-1,drop=FALSE]); storage.mode(dm) <- "double"; rownames(dm) <- ds[[1]]
stopifnot(identical(colnames(cm),ds[[1]]))
dds <- DESeqDataSetFromMatrix(cm,DataFrame(row.names=ds[[1]]),dm)
dds <- DESeq(dds,fitType="parametric",betaPrior=FALSE,minReplicatesForReplace=Inf,
             parallel=FALSE,quiet=TRUE)
r <- results(dds,name=tail(resultsNames(dds),1L),cooksCutoff=FALSE,independentFiltering=FALSE)
write.table(data.frame(gene=rownames(r),dispersion=dispersions(dds),baseMean=r$baseMean,
  log2FoldChange=r$log2FoldChange,lfcSE=r$lfcSE,stat=r$stat,pvalue=r$pvalue,padj=r$padj),
  args[[3]],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
