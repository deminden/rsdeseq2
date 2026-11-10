#!/usr/bin/env Rscript
# Independent posterior/tail references: original packages only, no rsdeseq2.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: generate_normal_mixture_posterior_reference.R OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2',
          as.character(packageVersion('ashr')) == '2.2.63',
          as.character(packageVersion('apeglm')) == '1.34.0')
numbers <- function(x) as.list(sprintf('%.17g', as.numeric(x)))
priors <- list(
  single = ashr::normalmix(1, 0, 1),
  point_zero = ashr::normalmix(1, 0, 0),
  point_shifted = ashr::normalmix(1, 1.25, 0),
  shrink = ashr::normalmix(c(.05,.2,.4,.35), c(0,0,0,0), c(0,.1,.7,4)),
  shifted = ashr::normalmix(c(.1,.3,.6), c(-2,.3,1.2), c(.2,0,3)),
  zero_weights = ashr::normalmix(c(0,.4,0,.6), c(0,-.7,4,1.1), c(0,.5,2,3)),
  tiny_weight = ashr::normalmix(c(1e-180,.25,.75), c(0,0,0), c(0,.3,5)),
  scale_range = ashr::normalmix(c(.2,.3,.5), c(-.2,.1,0), c(1e-5,.01,10)))
beta <- c(-40,-12,-6,-2,-.5,-1e-8,0,1e-8,.5,2,6,12,40)
se <- c(.01,.15,.7,2,8)
x <- rep(beta, each=length(se)); s <- rep(se, times=length(beta))
cutoffs <- c(-Inf,-40,-2,-.5,0,.5,1.25,2,40,Inf)
cases <- lapply(names(priors), function(name) {
  g <- priors[[name]]; data <- ashr::set_data(x,s)
  list(name=name, weights=numbers(g$pi), means=numbers(g$mean), sds=numbers(g$sd),
       estimates=numbers(x), standard_errors=numbers(s), thresholds=numbers(cutoffs),
       posterior_weights=numbers(ashr:::comp_postprob(g,data)),
       posterior_mean=numbers(ashr::postmean(g,data)), posterior_sd=numbers(ashr::postsd(g,data)),
       cdf=numbers(t(vapply(cutoffs, function(c) ashr::cdf_post(g,c,data), numeric(length(x))))),
       upper_tail=numbers(t(vapply(cutoffs, function(c) colSums(ashr:::comp_postprob(g,data) *
         stats::pnorm(c,ashr:::comp_postmean(g,data),ashr:::comp_postsd(g,data),lower.tail=FALSE)),numeric(length(x))))))
})
svalue_inputs <- list(empty=numeric(), single=.2, sorted=c(0,.01,.1,.5,1),
  ties=c(.2,.2,.1,.2,.1,0,1,.2), missing=c(NA,.3,NA,.1,.3,NA,0),
  all_missing=c(NA_real_,NA_real_), signed_zero=c(0,-0,.1,-0,0),
  fractional=rev((seq_len(257L) %% 37L)/37))
svalues <- lapply(names(svalue_inputs),function(name) {
  p <- svalue_inputs[[name]]
  list(name=name, rates=numbers(p), values=numbers(apeglm::svalue(p)))
})
jsonlite::write_json(list(source=list(R=as.character(getRversion()),
  ashr=as.character(packageVersion('ashr')),apeglm=as.character(packageVersion('apeglm'))),
  cases=cases,svalues=svalues),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
