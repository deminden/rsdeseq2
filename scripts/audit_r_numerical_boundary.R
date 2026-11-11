#!/usr/bin/env Rscript
# Lexical guard, supplemented by runtime numerical-call firewalls in identity tests.
args <- commandArgs(TRUE)
if (!length(args) %in% c(1L,2L)) stop('usage: R_SOURCE_DIRECTORY [FRESH_REPORT_JSON]')
forbidden <- c('sum','prod','mean','median','var','sd','quantile','lowess','loess','smooth.spline',
  'p.adjust','optim','optimize','optimise','integrate','uniroot','nlm','nlminb','log','log1p','log2','log10',
  'exp','expm1','sqrt','lgamma','gamma','digamma','trigamma','crossprod','tcrossprod','%*%','^',
  'rowSums','colSums','rowMeans','colMeans','weighted.mean','cov','cov.wt','cor',
  'qr','qr.solve','qr.coef','qr.fitted','qr.resid','solve','chol','chol2inv','svd','La.svd',
  'eigen','det','determinant','model.matrix','model.frame','poly','scale','eval','evalq','eval.parent','parse','match.fun')
distribution <- '^(d|p|q|r)(norm|t|chisq|f|gamma|beta|binom|nbinom|pois|unif|exp|cauchy|lnorm|logis|hyper|geom|weibull|signrank|wilcox)$'
name <- function(x) {
  if (is.symbol(x) || is.character(x)) return(as.character(x))
  if (is.call(x) && as.character(x[[1L]]) %in% c('::',':::')) return(paste(as.character(x[[2L]]),as.character(x[[3L]]),sep='::'))
  '<dynamic>'
}
bad <- function(x) {
  leaf <- sub('^.*::','',x)
  leaf %in% forbidden || grepl(distribution,leaf) || x %in% c('ashr::ash','ashr::ash.workhorse','apeglm::apeglm','DESeq2::DESeq','DESeq2::results')
}
scan <- function(expressions, file) {
  findings <- list(); templates <- 0L; calls <- 0L
  add <- function(call, scope, line) findings[[length(findings)+1L]] <<- list(file=file,line=line,scope=scope,call=call)
  walk <- function(x, scope, line) {
    if (missing(x)) return(invisible(NULL))
    if (is.expression(x) || is.pairlist(x)) { for (value in x) walk(value,scope,line); return(invisible(NULL)) }
    if (!is.call(x)) return(invisible(NULL))
    calls <<- calls+1L
    head <- name(x[[1L]])
    leaf <- sub('^.*::','',head)
    if (head %in% c('quote','expression','substitute','alist')) { templates <<- templates+1L; return(invisible(NULL)) }
    if (head %in% c('<-','=')) {
      if (is.symbol(x[[2L]]) && is.call(x[[3L]]) && identical(x[[3L]][[1L]],as.name('function'))) scope <- as.character(x[[2L]])
    }
    if (bad(head)) add(head,scope,line)
    if (head == '-' && length(x) == 3L && is.call(x[[2L]]) && is.call(x[[3L]]) &&
        all(vapply(x[2:3],function(value)sub('^.*::','',name(value[[1L]])) %in% c('ncol','NCOL','nrow','NROW'),logical(1))))
      add('statistical dimension subtraction',scope,line)
    if (head %in% c('::',':::') && bad(name(x))) add(name(x),scope,line)
    if (leaf == 'do.call' && !(identical(scope,'.rsdeseq2_frozen_binding') && identical(name(x[[2L]]),'base::substitute'))) add('unreviewed do.call',scope,line)
    positions <- c(lapply=2L,sapply=2L,vapply=2L,apply=3L,tapply=3L,mapply=1L,Map=1L,Reduce=1L,Filter=1L,Find=1L,Position=1L,Negate=1L)
    if (leaf %in% names(positions)) {
      values <- as.list(x)[-1L]; index <- if ('FUN' %in% names(values)) match('FUN',names(values)) else positions[[leaf]]
      if (length(values)>=index && bad(name(values[[index]]))) add(paste(head,'callback',name(values[[index]])),scope,line)
    }
    if (head == '<dynamic>' && !identical(scope,'.rsdeseq2_optional_object_accessor')) add('unreviewed dynamic function target',scope,line)
    pieces <- as.list(x)
    if (is.call(pieces[[1L]])) walk(pieces[[1L]],scope,line)
    for (piece in pieces[-1L]) walk(piece,scope,line)
  }
  for (x in expressions) {
    ref <- attr(x,'srcref'); line <- if(is.null(ref))NA_integer_ else as.integer(ref[[1L]])
    walk(x,'<top-level>',line)
  }
  list(findings=findings,quoted_templates=templates,executable_calls=calls)
}
# Real rejection controls distinguish executable calls/function pointers from syntax data.
negative_controls <- c('stats::pnorm(x)','sum(x)','lapply(x,sqrt)','do.call(f,x)','eval(x)','base::lapply(x,sqrt)','base::do.call(f,x)',
 'ncol(full) - ncol(reduced)')
for (source in negative_controls) {
  stopifnot(length(scan(parse(text=source),'<control>')$findings)>0L)
}
stopifnot(length(scan(parse(text='list(expression=quote({ sum(x); stats::pnorm(x) }))'),'<control>')$findings)==0L)
stopifnot(length(scan(parse(text='seq_len(ncol(x)) - 1L'),'<control>')$findings)==0L)
paths <- sort(list.files(args[[1L]],pattern='[.]R$',recursive=TRUE,full.names=TRUE))
if (!length(paths)) stop('no production R sources found')
records <- lapply(paths,function(path) scan(parse(path,keep.source=TRUE),path))
findings <- unlist(lapply(records,`[[`,'findings'),recursive=FALSE)
if (length(findings)) {
  for (f in findings) message(f$file,':',f$line,' ',f$scope,': ',f$call)
  stop('production wrapper contains unreviewed R numerical or dynamic calls')
}
report <- list(status='passed',scope='lexical production-R call guard; runtime identity/firewall tests remain separately required',
  forbidden_names=forbidden,distribution_pattern=distribution,source_files=length(paths),quoted_templates=sum(vapply(records,`[[`,integer(1),'quoted_templates')),
  executable_calls=sum(vapply(records,`[[`,integer(1),'executable_calls')),findings=findings,negative_controls=length(negative_controls),
  permitted_dynamic_target='.rsdeseq2_optional_object_accessor: reviewed object metadata getters',
  permitted_do_call='.rsdeseq2_frozen_binding: base::substitute syntax inspection')
if (length(args)==2L) {
  if(file.exists(args[[2L]]))stop('use a fresh report path')
  jsonlite::write_json(report,args[[2L]],auto_unbox=TRUE,pretty=TRUE,null='null')
}
cat(length(paths),'production R files pass the numerical-boundary lexical guard\n')
