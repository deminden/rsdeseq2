#!/usr/bin/env Rscript
# Focused frozen-VST process benchmark. R math below prepares identical fixture
# inputs for both packages; only the public wrapper call is the timed stage.
args <- commandArgs(TRUE)
if (length(args) == 3L && args[[1L]] == "compare") {
  baseline <- readRDS(args[[2L]])
  candidate <- readRDS(args[[3L]])
  stopifnot(identical(baseline$counts, candidate$counts),
    identical(baseline$sizeFactors, candidate$sizeFactors),
    identical(baseline$normalizationFactors, candidate$normalizationFactors),
    identical(attributes(baseline$transformed), attributes(candidate$transformed)),
    identical(is.na(baseline$transformed), is.na(candidate$transformed)),
    identical(is.nan(baseline$transformed), is.nan(candidate$transformed)))
  finite <- is.finite(baseline$transformed) & is.finite(candidate$transformed)
  stopifnot(all(finite), any(finite))
  error <- max(abs(baseline$transformed - candidate$transformed))
  stopifnot(error <= 2e-6)
  cat("Frozen VST inputs identical; maximum absolute transform change:", format(error, digits=17), "\n")
  quit(status=0)
}
if (length(args) != 7L) stop("usage: R_LIBRARY saved|expression size|matrix GENES SAMPLES OUTPUT_RDS STAGES_TSV; or compare BASELINE_RDS CANDIDATE_RDS")
library(rsdeseq2, lib.loc=normalizePath(args[[1L]]))
suppressPackageStartupMessages(library(DESeq2))
profile <- match.arg(args[[2L]], c("saved", "expression"))
normalization <- match.arg(args[[3L]], c("size", "matrix"))
genes <- as.integer(args[[4L]]); samples <- as.integer(args[[5L]])
stopifnot(genes >= 12L, samples >= 4L, samples %% 2L == 0L)
set.seed(20261110)
means <- exp(seq(log(10),log(1000),length.out=genes))
sf <- exp(seq(log(.7),log(1.4),length.out=samples))
counts <- matrix(rnbinom(genes*samples,mu=as.vector(outer(means,sf)),size=5),genes,samples)
dimnames(counts) <- list(paste0("gene",seq_len(genes)),paste0("sample",seq_len(samples)))
dds <- DESeqDataSetFromMatrix(counts,
  data.frame(condition=factor(rep(c("A","B"),each=samples/2L)),row.names=colnames(counts)),~condition)
nf <- NULL
if (normalization == "matrix") {
  nf <- outer(seq(.8,1.2,length.out=genes),sf)
  normalizationFactors(dds) <- nf
  sf <- NULL
} else sizeFactors(dds) <- sf
if (profile == "saved") {
  training <- exp(seq(log(.3),log(5000),length.out=60L))
  curve <- DESeq2:::localDispersionFit(training,.03+1.5/training+.015*(1+sin(log(training))),1e-8)
} else curve <- local({ asymptote <- .15; function(mu) asymptote+.3/sqrt(base::pmax(mu,1e-8)) })
attr(curve,"fitType") <- "local"
dispersionFunction(dds) <- curve
start <- proc.time()[["elapsed"]]
transformed <- rsdeseq2::varianceStabilizingTransformation(dds,blind=FALSE,fitType="local",nproc=1L)
elapsed <- proc.time()[["elapsed"]] - start
write.table(data.frame(stage="frozen_vst",elapsed_s=elapsed),args[[7L]],sep="\t",row.names=FALSE,quote=FALSE)
saveRDS(list(counts=counts,sizeFactors=sf,normalizationFactors=nf,
  transformed=SummarizedExperiment::assay(transformed)),args[[6L]],compress=FALSE)
