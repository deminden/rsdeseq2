#!/usr/bin/env Rscript
# Independent base-R probability equations for the frozen uniform-prior panel.
# Production rsdeseq2 is never loaded. This calibrates platform rounding in
# near-coincident pnorm log probabilities without changing comparison bounds.
args <- commandArgs(TRUE)
if (length(args) != 1L || file.exists(args[[1L]])) stop("Supply one fresh output TSV")
stopifnot(as.character(getRversion()) == "4.6.2")
profiles <- list(
  broad=list(w=1,a=-100,b=100), single_interval=list(w=1,a=-1,b=1),
  point_zero=list(w=1,a=0,b=0), point_shifted=list(w=1,a=1.25,b=1.25),
  symmetric=list(w=c(.1,.3,.4,.2),a=c(0,-.1,-2,-100),b=c(0,.1,2,100)),
  halfuniform=list(w=c(.05,.2,.3,.3,.15),a=c(0,-.2,-3,0,0),b=c(0,0,0,.2,100)),
  shifted=list(w=c(.1,.3,.4,.2),a=c(-1.5,-2,1,-100),b=c(-1.5,-.2,3,100)),
  zero_weights=list(w=c(0,.4,0,.6),a=c(0,-.5,4,-100),b=c(0,.5,5,100)),
  narrow=list(w=c(.2,.3,.5),a=c(-1e-6,.5,-100),b=c(1e-6,.500002,100)),
  ultra_narrow=list(w=c(.4,.6),a=c(-1e-12,-100),b=c(1e-12,100)))
x <- rep(c(-40,-12,-6,-2,-.5,-1e-8,0,1e-8,.5,2,6,12,40),each=5)
s <- rep(c(.01,.15,.7,2,8),times=13)
cutoffs <- c(-Inf,-40,-2,-.5,0,.5,1.25,2,40,Inf)
rows <- list()
for (name in names(profiles)) {
  p <- profiles[[name]]
  for (i in seq_along(x)) {
    la <- pnorm((x[i]-p$a)/s[i],log.p=TRUE)
    lb <- pnorm((x[i]-p$b)/s[i],log.p=TRUE)
    density <- la + log(1-exp(lb-la)) - log(p$b-p$a)
    point <- p$a == p$b
    density[point] <- dnorm((x[i]-p$a[point])/s[i],log=TRUE)-log(s[i])
    logweights <- density + log(p$w)
    weights <- exp(logweights-max(logweights))
    weights <- weights/sum(weights)
    for (j in seq_along(cutoffs)) {
      cutoff <- cutoffs[j]
      probability <- rep(1,length(p$w))
      probability[p$a>cutoff] <- 0
      inside <- p$a<=cutoff & p$b>cutoff
      lm <- pnorm((x[i]-cutoff)/s[i],log.p=TRUE)
      probability[inside] <- (exp(lm-la[inside])-1)/(exp(lb[inside]-la[inside])-1)
      fallback <- (cutoff-p$a)/(p$b-p$a)
      probability[is.nan(probability)] <- fallback[is.nan(probability)]
      probability[point] <- as.double(p$a[point]<=cutoff)
      rows[[length(rows)+1L]] <- data.frame(case=name,row=i-1L,cutoff=j-1L,
        estimate=sprintf("%.17g",x[i]),se=sprintf("%.17g",s[i]),threshold=sprintf("%.17g",cutoff),
        cdf=sprintf("%.17g",sum(weights*probability)))
    }
  }
}
write.table(do.call(rbind,rows),args[[1L]],sep="\t",quote=FALSE,row.names=FALSE)
cat("Independent uniform CDFs:",R.version$platform,as.character(getRversion()),"\n")
