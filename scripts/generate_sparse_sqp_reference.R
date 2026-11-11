#!/usr/bin/env Rscript
# Original DESeq2/ashr/mixsqp only; no candidate package or native engine.
args <- commandArgs(TRUE)
if (length(args) != 1L || file.exists(args[[1L]])) stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2',
  as.character(packageVersion('DESeq2')) == '1.52.0',
  as.character(packageVersion('ashr')) == '2.2.63',
  as.character(packageVersion('mixsqp')) == '0.3.54')
bits <- function(v) lapply(as.numeric(v), function(x) {
  paste0('bits:', paste(sprintf('%02x', as.integer(writeBin(x, raw(), size=8L,
    endian='big'))), collapse=''))
})
expected <- function(v) lapply(as.numeric(v), function(x) {
  if (is.na(x) && !is.nan(x)) 'NA' else sprintf('%.17g', x)
})
# The existing fitted-view regression's deterministic 40-gene input, generated
# here without importing a helper or any code from the candidate R package.
set.seed(20261014)
n <- 40L
mu <- exp(seq(log(10), log(500), length.out=n))
dispersion <- seq(.03, .35, length.out=n)
fold <- rep(c(1.5, 1, 2, .75), length.out=n)
counts <- t(vapply(seq_len(n), function(i) stats::rnbinom(8L,
  mu=c(rep(mu[[i]], 4L), rep(mu[[i]]*fold[[i]], 4L)),
  size=1/dispersion[[i]]), numeric(8L)))
rownames(counts) <- paste0('gene', seq_len(n))
colnames(counts) <- paste0('sample', seq_len(8L))
frame <- data.frame(condition=factor(rep(c('A','B'), each=4L)),
  row.names=colnames(counts))
fit <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(counts, frame, ~condition),
  fitType='mean', minReplicatesForReplace=Inf, quiet=TRUE)
raw <- DESeq2::results(fit, independentFiltering=FALSE)
subset <- seq.int(35L, 1L, by=-2L)
profiles <- list(
  default=list(rows=subset, control=list()),
  no_em=list(rows=subset, control=list(numiter.em=0)),
  explicit_14=list(rows=subset, control=list(maxiter.activeset=14)),
  explicit_20=list(rows=subset, control=list(maxiter.activeset=20)),
  automatic=list(rows=subset, control=list(maxiter.activeset=NULL)),
  one_sqp=list(rows=subset, control=list(maxiter.sqp=1)),
  reversed=list(rows=rev(subset), control=list()),
  full=list(rows=seq_len(n), control=list()),
  full_explicit_20=list(rows=seq_len(n), control=list(maxiter.activeset=20)))
.sparse_sqp_capture <- new.env()
suppressMessages(trace('estimate_mixprop', where=asNamespace('ashr'),
  print=FALSE, tracer=quote({
    capture <- get('.sparse_sqp_capture', .GlobalEnv)
    capture$data <- data; capture$g <- g; capture$prior <- prior
  })))
cases <- direct <- list()
for (name in names(profiles)) {
  p <- profiles[[name]]
  x <- raw$log2FoldChange[p$rows]; se <- raw$lfcSE[p$rows]
  warnings <- character()
  result <- withCallingHandlers(ashr::ash(x, se, mixcompdist='normal',
    method='shrink', optmethod='mixSQP', control=p$control, outputlevel=3),
    warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    })
  c <- .sparse_sqp_capture
  ll <- t(ashr:::log_comp_dens_conv.normalmix(c$g, c$data))
  cases[[length(cases)+1L]] <- list(name=paste0('sparse_',name),
    estimates=bits(x), standard_errors=bits(se), family='normal', point_mass=FALSE,
    nullweight=NULL, scales=NULL, mode=bits(0)[[1L]], range=bits(c(-Inf,Inf)),
    alpha=bits(0)[[1L]], prune_threshold=bits(1e-10)[[1L]], observation_weights=NULL,
    fixed_prior=FALSE, supplied_prior=NULL, df=NULL, thresholds=list(),
    optimizer='mixSQP', control=p$control, expected_error=NULL,
    posterior_mean=expected(result$result$PosteriorMean),
    posterior_sd=expected(result$result$PosteriorSD), svalues=expected(result$result$svalue),
    log_likelihood=expected(result$loglik)[[1L]], warnings=as.list(warnings))
  direct[[length(direct)+1L]] <- list(name=paste0('sparse_',name), optimizer='mixSQP',
    observations=nrow(ll), components=ncol(ll), log_likelihoods=bits(t(ll)),
    observation_weights=NULL, dirichlet=bits(c$prior), initial_weights=bits(c$g$pi),
    control=p$control, expected_error=NULL, expected_weights=expected(result$fitted_g$pi),
    expected_log_likelihood=expected(result$loglik)[[1L]],
    expected_penalized_log_likelihood=expected(result$fit_details$penloglik)[[1L]],
    converged=result$fit_details$optreturn$converged)
}
suppressMessages(untrace('estimate_mixprop', where=asNamespace('ashr')))
jsonlite::write_json(list(source=list(R=as.character(getRversion()),
  platform=R.version$platform, DESeq2=as.character(packageVersion('DESeq2')),
  ashr=as.character(packageVersion('ashr')), mixsqp=as.character(packageVersion('mixsqp')),
  optimizer_family='SQP'),
  cases=cases, em_fit_cases=direct), args[[1L]], auto_unbox=TRUE,
  pretty=TRUE, null='null', digits=NA)
