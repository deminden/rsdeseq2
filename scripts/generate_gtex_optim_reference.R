#!/usr/bin/env Rscript
# Portable objective-only R optim oracle, independent of the Rust optimizer.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, getRversion()=='4.6.1')
root <- normalizePath(args[[1]])
dir <- file.path(root,'crates/rsdeseq2/tests/data/gtex_v8_real_wald')
ref <- read.delim(file.path(dir,'compact/reference_results.tsv'))
cts <- read.delim(file.path(dir,'compact/counts.tsv'),check.names=FALSE)
x <- as.matrix(read.delim(file.path(dir,'design_matrix.tsv'),check.names=FALSE)[,-1])
sf <- read.delim(file.path(dir,'reference_size_factors.tsv'))$sizeFactor
eligible <- which(is.finite(ref$dispersion) & ref$baseMean > 1)
indices <- unique(c(eligible[which.max(ref$dispersion[eligible])],
                    eligible[which.max(abs(ref$stat[eligible]))],
                    eligible[which.max(ref$baseMean[eligible])]))
stopifnot(length(indices)==3L)
cases <- unlist(lapply(indices, function(i) lapply(c(TRUE,FALSE), function(compat) {
 y <- as.numeric(cts[i,-1]); alpha <- ref$dispersion[i]
 start <- c(log2(ref$baseMean[i]), rep(0,ncol(x)-1))
 control <- list(maxit=100L,factr=1e7,pgtol=0,lmm=5L,ndeps=rep(1e-3,ncol(x)))
 objective <- function(beta) {
   mu <- sf*2^as.vector(x%*%beta)
   if (compat) -sum(dnbinom(y, mu=mu, size=1/alpha, log=TRUE)) else
     sum((y+1/alpha)*log(1/alpha+mu) - log(1/alpha)/alpha - y*log(mu))
 }
 gradient <- function(beta) {
   mu <- sf*2^as.vector(x%*%beta)
   as.vector(-log(2)*crossprod(x, (y-mu)/(1+alpha*mu)))
 }
 result <- optim(start, objective, gr=if(compat) NULL else gradient, method='L-BFGS-B',lower=-30,upper=30,control=control)
 stopifnot(result$convergence==0L)
 list(r_optim_compat=compat,fixture=ref$gene[i],case_kind='GTEx_sparse_effect_expression_extreme',dispersion=alpha,
      counts=y,initial_par=start,lower=rep(-30,ncol(x)),upper=rep(30,ncol(x)),control=control,
      result=list(par=result$par,value=result$value,counts=as.list(result$counts)))
})),recursive=FALSE)
fixture <- list(provenance=list(R=R.version.string,method='stats::optim objective-only L-BFGS-B',
 input='gtex_v8_real_wald/compact; full-cohort dispersions and size factors'),
 contrasts=list(list(contrast='GTEx_kidney_blocked',samples=nrow(x),coefficients=ncol(x),
 design=lapply(seq_len(nrow(x)),function(i)unname(x[i,])),size_factors=sf,cases=cases)))
jsonlite::write_json(fixture,file.path(root,'crates/rsdeseq2/src/glm/irls/fixtures/gtex_optim_reference.json'),
 pretty=TRUE,auto_unbox=TRUE,digits=17)
