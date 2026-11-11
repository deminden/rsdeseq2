#!/usr/bin/env Rscript
# Original ashr only; freeze selection and selected numeric outputs independently.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion()) == '4.6.2', as.character(packageVersion('ashr')) == '2.2.63', !file.exists(args[[1L]]))
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~', ' ', script, fixed=TRUE)), '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
bits <- function(v) lapply(as.numeric(v), function(x) paste0('bits:', hex(x)))
expected <- function(v) if (is.null(v)) NULL else lapply(v, function(x) if (is.na(x) && !is.nan(x)) 'NA' else sprintf('%.17g', x))
g <- ashr::normalmix(c(.5,.3,.2), rep(0,3), c(0,.5,2))
profiles <- list(negative_infinity=-Inf, negative=-2, zero=0, fractional=.25, one=1,
  one_and_half=1.5, two=2, three=3, four=4, five=5, six=6, infinity=Inf,
  missing=NA_real_, nan=NaN, mean_only='PosteriorMean', sd_only='PosteriorSD',
  sign_only='svalue', prior_only='fitted_g', empty=character(), unknown='unknown',
  duplicate=c('PosteriorMean','PosteriorMean','fitted_g'), case_sensitive='posteriormean',
  all=c('PosteriorMean','PosteriorSD','svalue','qvalue','fitted_g'),
  missing_name=c(NA_character_,'PosteriorSD'), null=NULL, logical=TRUE,
  qvalue_only='qvalue', likelihood_only='loglik', flash_only='flash_data',
  sampler_only='post_sampler', numeric_empty=numeric(), numeric_multiple=c(1,2))
datasets <- list(complete=list(x=c(-3,-.2,0,.7,4), s=c(.3,.2,.5,1,.7)),
  missing=list(x=c(-1,NA_real_,NaN,.6,2), s=c(.2,.4,.8,0,Inf)))
cases <- list()
for (dataset in names(datasets)) for (name in names(profiles)) {
  d <- datasets[[dataset]]; output <- profiles[[name]]; warnings <- character()
  fit <- tryCatch(withCallingHandlers(ashr::ash(d$x,d$s,method='shrink',mixcompdist='normal',g=g,fixg=TRUE,outputlevel=output),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}), error=identity)
  cases[[length(cases)+1L]] <- list(name=paste(dataset,name,sep='_'), estimates=bits(d$x), standard_errors=bits(d$s),
    output=.formula_descriptor(output), threshold=bits(0)[[1L]], query_svalues=TRUE,
    expected_error=if(inherits(fit,'error'))conditionMessage(fit)else NULL,
    posterior_mean=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorMean),
    posterior_sd=if(inherits(fit,'error'))NULL else expected(fit$result$PosteriorSD),
    svalues=if(inherits(fit,'error'))NULL else expected(fit$result$svalue),
    fitted_prior=if(inherits(fit,'error'))NULL else !is.null(fit$fitted_g), warnings=as.list(warnings))
}
# Thresholded queries use the fitted prior even when ordinary svalue is absent.
for (level in c(1,4,5)) {
  x <- datasets$complete$x; s <- datasets$complete$s
  fit <- ashr::ash(x,s,method='shrink',mixcompdist='normal',g=g,fixg=TRUE,outputlevel=level)
  valid <- !is.null(fit$fitted_g) && !is.null(fit$result$PosteriorMean)
  sv <- if(valid) { a<-ashr::cdf_post(fit$fitted_g,.25,ashr::set_data(x,s));b<-ashr::cdf_post(fit$fitted_g,-.25,ashr::set_data(x,s));
    rates<-ifelse(fit$result$PosteriorMean>0,a,1-b);v<-numeric(length(rates));v[order(rates)]<-cumsum(sort(rates))/seq_along(rates);v } else NULL
  cases[[length(cases)+1L]] <- list(name=paste0('threshold_level_',level),estimates=bits(x),standard_errors=bits(s),
    output=.formula_descriptor(level),threshold=bits(.25)[[1L]],query_svalues=FALSE,
    expected_error=if(valid)NULL else 'threshold queries require PosteriorMean and fitted_g outputs',
    posterior_mean=expected(fit$result$PosteriorMean),posterior_sd=expected(fit$result$PosteriorSD),svalues=expected(sv),
    fitted_prior=!is.null(fit$fitted_g),warnings=list())
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),
  scope='lfcShrink consumed output fields; unavailable threshold dependencies reject explicitly'),cases=cases),
  args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
