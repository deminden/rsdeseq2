#!/usr/bin/env Rscript
# Complete real-interface comparisons, with all R numerical entry points blocked.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS')
driver <- normalizePath(args[[1L]]); library(rsdeseq2, lib.loc=normalizePath(args[[2L]]))
requireNamespace('ashr'); requireNamespace('apeglm')
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~', ' ', script, fixed=TRUE)), '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root, 'scripts/native_posterior_check_helpers.R'))
out <- args[[3L]]; if (dir.exists(out)) stop('Use a fresh output directory'); dir.create(out, recursive=TRUE)
model <- function(m) {
  for (name in intersect(c('y','scale','trials','capture_y','capture_trials'), names(m))) m[[name]] <- decode(m[[name]])
  m
}
encode_all <- function(v) {
  if (is.numeric(v)) return(encode(v))
  if (is.logical(v) || is.character(v)) return(as.list(unname(v)))
  if (is.list(v)) return(lapply(v, encode_all))
  v
}
canonical <- function(v) {
  if (!is.list(v)) return(v)
  if (!is.null(names(v))) v <- v[order(names(v))]
  lapply(v, canonical)
}
check <- function(input, wrapped, prefix) {
  directory <- file.path(out, prefix); dir.create(directory)
  path <- file.path(directory, 'input.json'); jsonlite::write_json(input, path, auto_unbox=TRUE, null='null', digits=NA)
  output <- system2(driver, shQuote(normalizePath(path)), stdout=TRUE, stderr=TRUE)
  if (!is.null(attr(output, 'status'))) stop(paste(output, collapse='\n'))
  writeLines(output, file.path(directory, 'standalone.json'))
  standalone <- jsonlite::fromJSON(paste(output, collapse='\n'), simplifyVector=FALSE)
  wrapped <- .formula_canonical(list(cases=wrapped()))
  jsonlite::write_json(wrapped, file.path(directory, 'wrapped.json'), auto_unbox=TRUE, pretty=TRUE, null='null')
  stopifnot(identical(canonical(wrapped), canonical(standalone)))
}
reference <- jsonlite::fromJSON(file.path(root, 'crates/rsdeseq2/tests/data/count_likelihood_r462.json'), simplifyVector=FALSE)
check(list(builtin_likelihood_cases=reference$kernel_cases), function() lapply(reference$kernel_cases, function(c) {
  r <- blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood', likelihood=model(c$builtin_likelihood),
    arguments=decode(c$arguments), lower=decode(c$lower), upper=decode(c$upper))))
  list(name=c$name, ok=TRUE, value=encode_all(r))
}), 'kernels')
modes <- jsonlite::fromJSON(file.path(root, 'crates/rsdeseq2/tests/data/count_mode_r462.json'), simplifyVector=FALSE)$cases
check(list(ashr_mode_cases=modes), function() lapply(modes, function(c) {
  x <- decode(c$estimates); se <- decode(c$standard_errors)
  arguments <- list(betahat=x, sebetahat=se, mixcompdist=c$family, pointmass=c$point_mass, prior='uniform',
    mixsd=if (is.null(c$scales)) NULL else decode(c$scales), mode=.formula_decode(c$mode_input), lik=model(c$builtin_likelihood))
  r <- blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='ashr_shrinkage', arguments=arguments,
    query_estimates=x, query_standard_errors=se, lfc_threshold=0, svalue=TRUE)))
  list(name=c$name, ok=TRUE, value=encode_all(r))
}), 'modes')
cached <- readRDS(args[[4L]]); dds <- cached$rust[seq_len(128L),]; n <- nrow(dds)
y <- rep(c(0,1,3,7,2), length.out=n); scale <- seq(.5,2,length.out=n)
likelihoods <- list(pois_identity=ashr::lik_pois(y,scale), pois_log=ashr::lik_pois(y,scale,'log'),
  binom_identity=ashr::lik_binom(y,9), binom_logit=ashr::lik_binom(y,9,'logit'))
public_errors <- list()
for (name in names(likelihoods)) for (variant in c('fixed_mode','default_mode','automatic_grid','default_fdr')) {
  options <- list(betahat=rep(0,n), sebetahat=rep(1,n), mixcompdist='uniform', lik=likelihoods[[name]])
  if (variant != 'automatic_grid') options$mixsd <- c(.2,.7,2)
  if (variant == 'fixed_mode') options$mode <- 0
  if (variant == 'default_fdr') options$method <- 'fdr'
  original <- do.call(ashr::ash, utils::modifyList(list(method='shrink'), options))
  actual <- blocked(do.call(rsdeseq2::lfcShrink, c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,quiet=TRUE), options)))
  stopifnot(identical(DESeq2::priorInfo(actual)$engine, 'Rust'))
  profile <- paste(name,variant,sep='_')
  saveRDS(list(original=original,actual=actual,options=options), file.path(out,paste0('public-',profile,'.rds')))
  message(profile, ': mean error ', max(abs(actual$log2FoldChange-original$result$PosteriorMean),na.rm=TRUE))
  public_errors[[profile]] <- list(mean=compare(actual$log2FoldChange,original$result$PosteriorMean,1e-3),
    sd=compare(actual$lfcSE,original$result$PosteriorSD,1e-3), svalue=compare(actual$svalue,original$result$svalue,.02))
}
jsonlite::write_json(list(kernel_batches=length(reference$kernel_cases),mode_batches=length(modes),public_errors=public_errors,
  R_calculations_blocked=TRUE,existing_bounds_unchanged=TRUE), file.path(out,'validation.json'), auto_unbox=TRUE, pretty=TRUE, digits=NA)
cat(length(reference$kernel_cases),'kernel and',length(modes),'mode batches have exact wrapped/native bits; sixteen public profiles pass\n')
