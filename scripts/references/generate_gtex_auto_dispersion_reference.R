#!/usr/bin/env Rscript
# Preserve full-dataset normalization while replaying a compact GTEx regression.
args <- commandArgs(TRUE)
fixture <- if (length(args)) args[[1]] else 'crates/rsdeseq2/tests/data/gtex_auto_dispersion'
suppressPackageStartupMessages(library(DESeq2))
stopifnot(as.character(packageVersion('DESeq2')) == '1.52.0')
read_tsv <- function(name) read.delim(file.path(fixture,name),check.names=FALSE)
ct <- read_tsv('counts.tsv'); dt <- read_tsv('design_matrix.tsv')
sf <- read_tsv('reference_size_factors.tsv')
cm <- as.matrix(ct[,-1]); storage.mode(cm)<-'integer'; rownames(cm)<-ct[[1]]
dm <- as.matrix(dt[,-1]); storage.mode(dm)<-'double'
stopifnot(identical(colnames(cm),dt[[1]]),identical(colnames(cm),sf[[1]]))
make_dds <- function() {
  dds <- DESeqDataSetFromMatrix(cm,S4Vectors::DataFrame(row.names=colnames(cm)),dm)
  sizeFactors(dds) <- sf[[2]]
  dds
}
write_tsv <- function(x,name) write.table(x,file.path(fixture,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
for (mode in c('auto','glm','weighted')) {
  dds <- make_dds()
  if (mode == 'weighted') assays(dds,withDimnames=FALSE)[['weights']] <- matrix(1,nrow(cm),ncol(cm))
  dds <- estimateDispersionsGeneEst(dds,quiet=TRUE,linearMu=if(mode=='glm') FALSE else NULL)
  write_tsv(data.frame(gene=rownames(dds),dispGeneEst=mcols(dds)$dispGeneEst,
      dispGeneIter=mcols(dds)$dispGeneIter),paste0('reference_',mode,'.tsv'))
  mu<-assays(dds)[['mu']]; colnames(mu)<-colnames(cm)
  write_tsv(data.frame(gene=rownames(dds),mu,check.names=FALSE),paste0('reference_',mode,'_mu.tsv'))
}
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(fixture,'sessionInfo.txt'))
