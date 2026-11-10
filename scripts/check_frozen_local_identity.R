#!/usr/bin/env Rscript
# Compare actual production wrapper calls with a process that cannot load R.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR")
driver <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc = normalizePath(args[[2L]]))
suppressPackageStartupMessages(library(DESeq2))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), ".."))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
out <- args[[3L]]
if (dir.exists(out)) stop("Use a fresh identity output directory")
dir.create(out, recursive=TRUE)
numbers <- function(x) as.list(sprintf("%.17g", as.numeric(x)))
bits <- function(x) lapply(unname(x), hex)
ast <- function(x) {
  if (x$type == "literal") x["value"] <- list(.formula_descriptor(x$value))
  if (x$type == "call") x$arguments <- lapply(x$arguments, function(a) { a$expression <- ast(a$expression); a })
  x
}
modelInput <- function(model) {
  if (model$type == "saved") return(list(type="saved", knots=numbers(model$knots),
    residual_values=numbers(model$coefficients[,1L]), residual_slopes=numbers(model$coefficients[,2L]),
    polynomial=numbers(model$polynomial[1:4])))
  list(type="expression", argument=model$argument, expression=ast(model$expression),
       constants=lapply(model$constants,.formula_descriptor))
}
blocked <- function(expression) {
  functions <- list(base=c("sinh", "asinh", "exp", "sqrt", "mean", "sweep", "/"), stats="predict")
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name, where=asNamespace(package))), add=TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(trace(name, where=asNamespace(package), print=FALSE,
      tracer=quote(stop("production frozen VST executed an R calculation", call.=FALSE))))
  force(expression)
}
guard <- tryCatch(blocked(sinh(1)), error=identity)
stopifnot(inherits(guard,"error"),grepl("executed an R calculation",conditionMessage(guard)))

counts <- matrix(as.integer(round(exp(seq(0,log(1e6),length.out=48L)))),12L,4L)
counts[1L,1L] <- 0L
dimnames(counts) <- list(paste0("gene",1:12),paste0("sample",1:4))
storage.mode(counts) <- "double" # The private grid adapter receives validated double storage.
means <- exp(seq(log(.3),log(5000),length.out=60L))
trend <- DESeq2:::localDispersionFit(means,.03+1.5/means+.015*(1+sin(log(means))),1e-8)
serialized <- unserialize(serialize(trend,NULL))
custom <- local({ a <- .15; function(mu) a+.3/sqrt(base::pmax(mu,1e-8)) })
floor <- function(mu) rep(1e-8,length(mu))
curves <- list(saved=trend, reloaded=serialized, custom=custom, minimum=floor)
inputs <- list(); actual <- list()
for (name in names(curves)) for (normalization in c("size_factors", "matrix")) {
  curve <- curves[[name]]
  attr(curve,"fitType") <- "local"
  dds <- DESeqDataSetFromMatrix(counts,
    data.frame(condition=factor(c("A","A","B","B")),row.names=colnames(counts)),~condition)
  sf <- c(.7,1.1,1.4,.9); nf <- NULL
  if (normalization == "matrix") {
    nf <- outer(seq(.8,1.2,length.out=12L),sf)
    normalizationFactors(dds) <- nf
    sf <- NULL
  } else sizeFactors(dds) <- sf
  dispersionFunction(dds) <- curve
  model <- rsdeseq2:::.rsdeseq2_frozen_local_model(curve)
  identifier <- paste(name,normalization,sep="_")
  inputs[[length(inputs)+1L]] <- list(name=identifier,genes=nrow(counts),samples=ncol(counts),
    counts=as.list(as.integer(t(counts))),size_factors=if(is.null(sf))NULL else numbers(sf),
    normalization_factors=if(is.null(nf))NULL else numbers(t(nf)),model=modelInput(model))
  output <- blocked({
    grid <- rsdeseq2:::.rsdeseq2_frozen_local_curve(model,counts,sf,nf)
    transformed <- rsdeseq2::varianceStabilizingTransformation(dds,blind=FALSE,fitType="local")
    list(grid=grid,transformed=SummarizedExperiment::assay(transformed))
  })
  actual[[length(actual)+1L]] <- list(name=identifier,curve=bits(output$grid),transformed=bits(as.vector(t(output$transformed))))
}
path <- file.path(out,"input.json")
jsonlite::write_json(list(frozen_local_cases=inputs),path,auto_unbox=TRUE,null="null",digits=NA)
response <- system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(response,"status"))) stop(paste(response,collapse="\n"))
expected <- jsonlite::fromJSON(paste(response,collapse="\n"),simplifyVector=FALSE)
actual <- .formula_canonical(list(cases=actual))
canonical <- function(x) {
  if (!is.list(x)) return(x)
  if (!is.null(names(x))) x <- x[sort(names(x))]
  lapply(x, canonical)
}
writeLines(response,file.path(out,"standalone.json"))
jsonlite::write_json(actual,file.path(out,"wrapped.json"),auto_unbox=TRUE,null="null",pretty=TRUE)
if (!identical(canonical(actual),canonical(expected))) stop(paste(all.equal(canonical(actual),canonical(expected),tolerance=0),collapse="\n"))
cat(length(actual$cases),"complete frozen VST and 999-point curve cases: every output bit identical; R calculations blocked\n")
