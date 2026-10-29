#!/usr/bin/env Rscript
# Independent scalar probability references. No production Rust/R wrapper runs.
if (getRversion() != "4.6.2") stop("Requires the pinned R 4.6.2 reference runtime")
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: OUTPUT_TSV")
hex <- function(x) paste(format(writeBin(as.double(x), raw(), 8L, endian = "big")), collapse = "")
rows <- list()
add <- function(family, q, df = NA_real_, samples = 0L, coefficients = 0L, expected) {
  rows[[length(rows)+1L]] <<- data.frame(family = family, input_bits = hex(q),
    df_bits = hex(df), samples = samples, coefficients = coefficients, expected_bits = hex(expected))
}
q <- unique(c(seq(-12,12,by=.125), -37.5,-30,-20, -1e-8, 1e-8,
              1/sqrt(.5), -1/sqrt(.5), 20,30,37.5,
              unlist(lapply(sqrt(2)*c(.84375,1.25),function(b) b+c(-1,0,1)*.Machine$double.eps))))
for (x in q) {
  add("normal_upper", x, expected = pnorm(x, lower.tail = FALSE))
  add("normal_two_sided", x, expected = 2*pnorm(abs(x), lower.tail = FALSE))
}
tq <- unique(c(-1e6,-100,-12,-4,-1/sqrt(.5),-.125,0,1e-8,.125,.25,.5,1,
               1/sqrt(.5),2,4,8,12,30,100,1e6,1e154,1e308))
for (df in c(1e-12,1e-6,.1,.5,1,1.5,2,3,4,5,8,10,12,20,30,63.5,63.999999,64,64.000001,65,100,1000,1e6,1e12,1e30,1e300))
  for (x in tq) add("t_upper", x, df, expected = pt(x,df,lower.tail=FALSE))
for (df in c(.1,.5,1,2,3,5,10,30,100,1000,1e6))
  for (x in unique(c(0,.01,.1,1,2,6,8,40,100,df,df+sqrt(2*df),pmax(0,df-c(1,2,4,8,12)*sqrt(2*df)))))
    add("chi_square",x,df,expected=pchisq(x,df,lower.tail=FALSE))
for (m in c(4L,8L,16L,52L,61L,1000L,400001L,400002L,800001L,800002L,1000000L))
  for (p in unique(c(1L,2L,as.integer(m/2L),m-1L)))
    if (p<m) add("cooks",0,samples=m,coefficients=p,expected=qf(.99,p,m-p))
write.table(do.call(rbind,rows),args[[1L]],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
writeLines(c(trimws(capture.output(sessionInfo()),which="right"),paste("R executable:",file.path(R.home(),"bin","R"))),paste0(args[[1L]],".session.txt"))
