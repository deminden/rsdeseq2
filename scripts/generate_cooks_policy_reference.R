#!/usr/bin/env Rscript
# Independent oracle: original DESeq2 fits/results, never the rsdeseq2 wrapper.
suppressPackageStartupMessages(library(DESeq2))
stopifnot(getRversion() == "4.6.2", packageVersion("DESeq2") == "1.52.0")
arguments <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
script <- gsub("~+~", " ", script, fixed=TRUE)
root <- normalizePath(file.path(dirname(script), ".."))
paths <- if(length(arguments)) arguments else c(
  file.path(root,"crates/rsdeseq2/tests/data/cooks_policy_r462.json"),
  file.path(root,"r-pkg/rsdeseq2/tests/testthat/data/cooks_policy_r462.json"))
if(any(file.exists(paths))) stop("Use fresh output paths; frozen references are not overwritten")
hex <- function(v) if(is.na(v)) NULL else paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
hexes <- function(v) lapply(unname(v),hex)
strings <- function(v) as.list(sprintf("%.17g",v))
set.seed(20261030)
genes <- 40L; samples <- 6L
counts <- matrix(rnbinom(genes*samples,mu=rep(seq(20,500,length.out=genes),samples),size=5),genes,samples)
counts[1,] <- c(1,5,6,9,7,8)
counts[2,] <- c(9,5,6,1,7,8)
counts[3:4,] <- rep(counts[1,],each=2)
rownames(counts) <- paste0("gene",seq_len(genes)); colnames(counts) <- paste0("sample",seq_len(samples))
cooks <- matrix(0,genes,samples,dimnames=dimnames(counts))
cooks[1:2,] <- rep(c(10,.1,.2,10,.3,.4),each=2)
cooks[3,] <- c(10,.1,.2,.3,.4,.5)
cooks[4,] <- c(.1,.2,.3,10,.4,.5)
maximum <- c(rep(10,4),rep(0,genes-4))
cd <- data.frame(condition=factor(rep(c("A","B"),each=3)),numeric=rep(c(0,1),each=3),
  exposure=rep(1,samples),nuisance=seq_len(samples),batch=factor(rep(c("X","Y"),3)),
  tri=factor(rep(c("L","M","N"),2)),row.names=colnames(counts))
specifications <- list(
  factor=list(formula=~condition),
  releveled=list(formula=~condition,relevel=TRUE),
  function_factor=list(formula=~relevel(condition,"B")),
  no_intercept=list(formula=~0+condition),
  numeric=list(formula=~numeric),
  numeric_identity=list(formula=~I(numeric)),
  offset=list(formula=~condition+offset(log(exposure))),
  removed=list(formula=~condition+nuisance-nuisance),
  matrix=list(formula=~condition,matrix=TRUE),
  multilevel=list(formula=~tri),
  additive=list(formula=~batch+condition),
  factor_numeric=list(formula=~factor(numeric)))
cases <- lapply(names(specifications),function(name) {
  spec <- specifications[[name]]; data <- cd
  if(isTRUE(spec$relevel)) data$condition <- relevel(data$condition,"B")
  design <- if(isTRUE(spec$matrix)) model.matrix(spec$formula,data) else spec$formula
  dds <- DESeqDataSetFromMatrix(counts,data,spec$formula)
  sizeFactors(dds) <- rep(1,samples); dispersions(dds) <- rep(.2,genes)
  fitted <- nbinomWaldTest(dds,betaPrior=FALSE,quiet=TRUE)
  raw <- DESeq2::results(fitted,cooksCutoff=FALSE,independentFiltering=FALSE)
  stopifnot(all(is.finite(raw$pvalue[1:4])))
  # Matrix provenance uses the identical fitted numeric model, changing only
  # the design representation inspected by results(), not its coefficients.
  if(is.matrix(design)) design(fitted) <- design
  # Supply diagnostics explicitly to isolate row policy, independent of how an
  # optimizer estimated Cook's values. Counts/coefficients/p-values stay fitted.
  assays(fitted)[["cooks"]] <- cooks
  mcols(fitted)$maxCooks <- maximum
  expected <- DESeq2::results(fitted,cooksCutoff=5,independentFiltering=FALSE)
  data <- as.data.frame(colData(fitted))
  numericDesign <- if(is.matrix(design)) design else model.matrix(design,data)
  factors <- lapply(names(data)[vapply(data,is.factor,logical(1))],function(field)
    list(name=field,levels=nlevels(data[[field]]),codes=as.list(as.integer(data[[field]]))))
  # Eligibility is observed from the uniquely low-count third row's result,
  # rather than reimplementing the policy predicate as our own reference.
  list(name=name,kind=if(is.matrix(design))"matrix"else"formula",
    formula=if(is.matrix(design))NULL else paste(deparse(design),collapse=" "),
    variables=if(is.matrix(design))list()else as.list(all.vars(design)),
    coefficients=ncol(numericDesign),design=strings(as.vector(t(numericDesign))),factors=factors,
    raw_pvalue=hexes(raw$pvalue),pvalue=hexes(expected$pvalue),padj=hexes(expected$padj),
    observed_eligible=!is.na(expected$pvalue[3]),
    masked=as.list(is.na(expected$pvalue)&!is.na(raw$pvalue)))
})
stopifnot(identical(vapply(cases,function(x)x$observed_eligible,logical(1)),
  c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)))
stopifnot(!cases[[1]]$masked[[1]],cases[[1]]$masked[[2]])
reference <- list(metadata=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion("DESeq2")),
  seed=20261030L,diagnostics_supplied=TRUE,purpose="Cook row policy on original DESeq2 fits/results; no Rust reference path"),
  genes=genes,samples=samples,counts=as.list(as.integer(t(counts))),cooks=strings(as.vector(t(cooks))),
  maximum=strings(maximum),cutoff="5",cases=cases)
for(path in paths) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(reference,path,auto_unbox=TRUE,null="null",pretty=TRUE)
  cat("Wrote",path,"\n")
}
