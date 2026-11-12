#!/usr/bin/env Rscript
# Independent reference only: never load rsdeseq2.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: generate_truncated_posterior_reference.R FRESH_OUTPUT_JSON')
if (file.exists(args[[1L]])) stop('Use a fresh path; references are never overwritten')
stopifnot(as.character(getRversion()) == '4.6.2',
          as.character(packageVersion('ashr')) == '2.2.63')
numbers <- function(x) as.list(sprintf('%.17g', as.numeric(x)))
bits <- function(x) lapply(unname(as.numeric(x)), function(v) {
  paste0('bits:', paste(format(writeBin(v, raw(), size=8, endian='big')), collapse=''))
})
profiles <- list(
  half_normal=list(w=c(.1,.2,.35,.35), m=c(0,0,0,0), sd=c(0,.2,2,2),
                   a=c(-Inf,-Inf,-Inf,0), b=c(0,0,0,Inf)),
  asymmetric_halves=list(w=c(.2,.4,.4), m=c(0,0,0), sd=c(.1,1,8),
                         a=c(-Inf,0,0), b=c(0,Inf,Inf)),
  bounded=list(w=c(.3,.7), m=c(0,0), sd=c(.5,3), a=c(-1,-2), b=c(1,4)),
  shifted=list(w=c(.3,.7), m=c(-.7,1.2), sd=c(.5,3), a=c(-2,-1), b=c(1,4)),
  zero_sd=list(w=c(.2,.8), m=c(-.5,0), sd=c(0,2), a=c(-1,-Inf), b=c(0,Inf)),
  equal_bounds=list(w=c(.2,.8), m=c(0,0), sd=c(1,2), a=c(.3,-Inf), b=c(.3,Inf)),
  zero_weights=list(w=c(0,.3,.7), m=c(0,0,0), sd=c(.01,.4,4),
                    a=c(-Inf,-Inf,0), b=c(0,0,Inf)))
x <- rep(c(-40,-6,-.5,0,.5,6,40), each=4L)
s <- rep(c(.01,.15,1,8), times=7L)
# A wide untruncated component keeps every cutoff inside some component. The
# installed component CDF otherwise raises an upstream uninitialized-B error.
cutoffs <- c(-Inf,-2,-.3,0,.3,2,Inf)
cases <- lapply(names(profiles), function(name) {
  p <- profiles[[name]]
  g <- ashr::tnormalmix(p$w,p$m,p$sd,p$a,p$b)
  data <- ashr::set_data(x,s)
  component_cdf <- function(c) {
    # Invoke each supplied component in a common augmented prior, without
    # changing its posterior weight or its component-level arithmetic.
    augmented <- ashr::tnormalmix(c(p$w,0),c(p$m,0),c(p$sd,1),
                                 c(p$a,-Inf),c(p$b,Inf))
    cp <- ashr:::comp_cdf_post.tnormalmix(augmented,c,data)
    colSums(cp[seq_along(p$w),,drop=FALSE] * ashr:::comp_postprob(g,data))
  }
  list(name=name,weights=bits(p$w),means=bits(p$m),sds=bits(p$sd),
       lower=bits(p$a),upper=bits(p$b),estimates=bits(x),standard_errors=bits(s),
       thresholds=bits(cutoffs),
       log_likelihoods=numbers(ashr:::log_comp_dens_conv.tnormalmix(g,data)),
       posterior_weights=numbers(ashr:::comp_postprob(g,data)),
       posterior_mean=numbers(ashr::postmean(g,data)),
       posterior_sd=numbers(ashr::postsd(g,data)),
       cdf=numbers(t(vapply(cutoffs,component_cdf,numeric(length(x))))))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),
  ashr=as.character(packageVersion('ashr')),input_encoding='f64_bits',
  cdf_adapter='zero-weight untruncated sentinel avoids upstream uninitialized-B error'),
  cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
