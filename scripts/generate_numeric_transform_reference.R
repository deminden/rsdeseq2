#!/usr/bin/env Rscript
# Independent formula-transform references; never loads rsdeseq2.
stopifnot(getRversion()=="4.6.2", packageVersion("DESeq2")=="1.52.0")
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), ".."))
paths <- if(length(args)) args else file.path(root,
  "crates/rsdeseq2/tests/data/numeric_transform_r462.json")
if(any(file.exists(paths))) stop("Use fresh paths; references are never overwritten")
hex <- function(v) paste(format(writeBin(as.double(v), raw(), size=8, endian="big")), collapse="")
hexes <- function(v) lapply(unname(v), hex)
coefs <- function(x) if(is.null(x)) NULL else list(alpha=hexes(x$alpha), norm2=hexes(x$norm2))
cases <- list()
record <- function(name, specification, evaluate) {
  value <- tryCatch(suppressWarnings(evaluate()), error=identity)
  valid <- !inherits(value,"error") && all(is.finite(value))
  cases[[length(cases)+1L]] <<- c(list(name=name), specification, list(
    expected=if(valid) hexes(as.vector(t(as.matrix(value)))) else NULL,
    coefficients=if(valid) coefs(attr(value,"coefs")) else NULL,
    center=if(valid && !is.null(attr(value,"scaled:center"))) hexes(attr(value,"scaled:center")) else NULL,
    scale=if(valid && !is.null(attr(value,"scaled:scale"))) hexes(attr(value,"scaled:scale")) else NULL))
}
polynomial <- function(name, x, degree, raw=FALSE, fit=NULL) {
  record(name, list(operation="poly", input=hexes(x), rows=length(x), columns=1L,
    degree=degree, raw=raw, supplied_coefficients=coefs(fit)),
    function() stats::poly(x, degree=degree, raw=raw, coefs=fit))
}
irregular <- c(-3.5,-2,-1.125,-.5,0,.0625,.25,.75,1.5,2.25,4,7,9,11,13,15,18,24,30,41)
for(d in 1:16) polynomial(paste0("irregular_degree_",d), irregular, d)
for(d in 1:5) polynomial(paste0("shifted_degree_",d), 1e12+c(0,.125,.25,.5,1,1.75,3,5), d)
for(d in 1:3) polynomial(paste0("repeated_degree_",d), c(-1,-1,0,0,2,2,5,5), d)
for(n in c(32L,53L,64L,68L,85L,91L,95L)) polynomial(paste0("long_degree16_",n), seq_len(n),16L)
polynomial("cancellation",c(1e16,1,-1e16,3),2L)
polynomial("raw_fractional",irregular,7L,TRUE)
polynomial("raw_negative_zero",c(-0,-.5,0,.5,1,2),4L,TRUE)
polynomial("raw_degree16",c(-1.25,-.75,-.1,.1,.75,1.25),16L,TRUE)
polynomial("raw_overflow",c(1e200,2e200),3L,TRUE)
polynomial("too_many_degrees",c(0,0,1,2),3L)
polynomial("zero_degree",irregular,0L)
polynomial("nonfinite",c(0,1,Inf,3),2L)
polynomial("numerical_rank",seq_len(30L),25L)
fit <- attr(stats::poly(irregular,degree=5L),"coefs")
polynomial("reuse_on_new_samples",c(-2.75,.33,1.25,6.125,44.5),5L,fit=fit)
polynomial("reuse_lower_degree",c(-2.75,.33,1.25,6.125,44.5),2L,fit=fit)
polynomial("reuse_single_sample",.33,5L,fit=fit)
policy <- function(x) if(is.logical(x)) list(mode=if(x)"automatic" else "disabled",values=NULL) else list(mode="supplied",values=hexes(x))
scaling <- function(name,x,center=TRUE,scale=TRUE) {
  x <- as.matrix(x)
  record(name,list(operation="scale",input=hexes(as.vector(t(x))),rows=nrow(x),columns=ncol(x),
    center_policy=policy(center),scale_policy=policy(scale)),function() base::scale(x,center=center,scale=scale))
}
for(i in 1:4) {
  x <- switch(i,irregular,1e12+c(0,.25,.5,1.25),c(1e16,1,-1e16,3),c(-0,.25,-.5,2))
  scaling(paste0("scale_default_",i),x)
  scaling(paste0("scale_uncentered_",i),x,FALSE,TRUE)
  scaling(paste0("scale_center_only_",i),x,TRUE,FALSE)
  scaling(paste0("scale_identity_",i),x,FALSE,FALSE)
}
m <- cbind(c(-3,-.5,1,2,4,8),c(1e12,1e12+.25,1e12+.5,1e12+1,1e12+2,1e12+4))
scaling("scale_matrix",m)
scaling("scale_explicit",m,c(1,1e12),c(-2,4))
scaling("scale_explicit_center",m,c(1,1e12),TRUE)
scaling("scale_explicit_divisor",m,FALSE,c(-2,4))
scaling("scale_single_uncentered",3,FALSE,TRUE)
scaling("scale_zero",rep(4,5))
scaling("scale_zero_divisor",m,FALSE,c(0,1))
scaling("scale_bad_controls",m,c(1,2,3),FALSE)
scaling("scale_large_identity",c(1e308,1e308),FALSE,FALSE)
for(x in c(-3.9,-1.1,-.9,-0,0,.9,1.1,3.9,-2147483647.9,2147483647.9,-2147483648,2147483648,Inf,NaN)) {
  record(paste0("integer_",length(cases)+1L),list(operation="integer",input=hexes(x),rows=1L,columns=1L),function() as.integer(x))
}
scaling("scale_negative_divisor",c(-3,-.5,1,2,4,8),TRUE,-2)
scaling("scale_negative_zeros_center_only",c(-0,-0),TRUE,FALSE)
scaling("scale_negative_zeros_supplied",c(-0,-0),FALSE,-2)
reference <- list(metadata=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion("DESeq2")),
  purpose="Independent finite formula covariates, polynomial QR/recurrence, scaling and integer coercion; invalid fitting outputs reject"),cases=cases)
for(path in paths) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(reference,path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA)
  cat("Wrote",length(cases),"cases to",path,"\n")
}
