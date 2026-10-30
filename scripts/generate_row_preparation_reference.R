#!/usr/bin/env Rscript
# Independent R arithmetic/row decisions; never loads rsdeseq2.
stopifnot(getRversion()=="4.6.2",packageVersion("DESeq2")=="1.52.0")
args <- commandArgs(TRUE)
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
script <- gsub("~+~"," ",script,fixed=TRUE)
root <- normalizePath(file.path(dirname(script),".."))
paths <- if(length(args))args else c(file.path(root,"crates/rsdeseq2/tests/data/row_preparation_r462.json"),
  file.path(root,"r-pkg/rsdeseq2/tests/testthat/data/row_preparation_r462.json"))
if(any(file.exists(paths)))stop("Use fresh paths; reference files are never overwritten")
hex <- function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
hexes <- function(v)lapply(unname(v),hex)
specs <- list(
 size=list(sf=c(.7,1,1.4,2),offset=c(log(.8),0,log(1.2),1e-8)),
 matrix=list(nf=matrix(seq(.7,2,length.out=12),3,4),offset=c(-.2,.1,0,.3)),
 precedence=list(sf=c(2,3,4,5),nf=matrix(seq(.7,2,length.out=12),3,4),offset=c(-.2,.1,0,.3)),
 zero_offset=list(sf=c(.7,1,1.4,2),offset=rep(0,4)),
 tiny_base_rescued=list(sf=rep(1e-310,4),offset=rep(700,4)),
 subnormal=list(sf=rep(1,4),offset=rep(-745,4)),
 overflow_product=list(sf=rep(1e300,4),offset=rep(100,4)),
 underflow_product=list(sf=rep(1e-300,4),offset=rep(-700,4)),
 offset_underflow=list(sf=rep(1e300,4),offset=rep(-1000,4)),
 offset_overflow=list(sf=rep(1e-300,4),offset=rep(1000,4)),
 infinite_offset=list(sf=rep(1,4),offset=c(0,0,Inf,0)),
 missing_offset=list(sf=rep(1,4),offset=c(0,0,NA_real_,0)))
offsets <- lapply(names(specs),function(name) {
 x<-specs[[name]];base<-if(is.null(x$nf))matrix(rep(x$sf,each=3),3,4)else x$nf
 value<-sweep(base,2,exp(x$offset),"*")
 valid<-all(is.finite(x$offset))&&all(is.finite(value)&value>0)
 list(name=name,genes=3L,samples=4L,size_factors=if(is.null(x$sf))NULL else hexes(x$sf),
   factors=if(is.null(x$nf))NULL else hexes(as.vector(t(x$nf))),offsets=hexes(x$offset),
   expected=if(valid)hexes(as.vector(t(value)))else NULL)
})
ct<-matrix(c(0,0,0,0,1,2,3,4,0,0,0,5,5,6,7,8,4294967295,0,1,0,0,0,0,0,1,1,1,1),7,4,byrow=TRUE)
d<-c(NA,.1,.2,.3,.4,NaN,.5);f<-c(-Inf,.2,.3,.4,.5,NA,.6)
rowSpecs<-list(raw=list(ct=ct,d=d,f=f),
 excluded=list(ct=ct,d=replace(d,4,NA),f=replace(f,4,NaN),excluded=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)),
 reordered=list(ct=ct[c(4,2,1,7),],d=c(NA,.1,NA,.5),f=c(NaN,.2,-Inf,.6),excluded=c(TRUE,FALSE,FALSE,FALSE)),
 bad_dispersion=list(ct=ct,d=replace(d,2,NA),f=f),
 all_zero=list(ct=matrix(0,3,4),d=rep(NA_real_,3),f=rep(NA_real_,3)),
 all_excluded=list(ct=matrix(1,3,4),d=rep(NA_real_,3),f=rep(NA_real_,3),excluded=rep(TRUE,3)))
rows<-lapply(names(rowSpecs),function(name) {
 x<-rowSpecs[[name]];zero<-rowSums(x$ct)==0
 if(!is.null(x$excluded))zero<-zero|x$excluded
 index<-which(!zero)
 valid<-length(index)>0&&all(is.finite(x$d[index])&x$d[index]>0&is.finite(x$f[index])&x$f[index]>0)
 list(name=name,genes=nrow(x$ct),samples=ncol(x$ct),counts=as.list(as.vector(t(x$ct))),
  dispersions=hexes(x$d),fitted=hexes(x$f),excluded=if(is.null(x$excluded))NULL else as.list(x$excluded),
  expected=if(valid)as.list(index-1L)else NULL)
})
reference<-list(metadata=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion("DESeq2")),
 purpose="Independent base R offset products and row selection; original DESeq2 weights exclusion is separately tested"),
 offsets=offsets,rows=rows)
for(path in paths){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(reference,path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA);cat("Wrote",path,"\n")}
