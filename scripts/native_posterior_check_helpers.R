# Shared numerical firewall and primitive encoding for local posterior checks.
decode <- function(v) vapply(v, function(h) {
  h <- sub('^bits:', '', h)
  readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)), 'double', n=1L, size=8L, endian='big')
}, numeric(1L))
encode <- function(v) lapply(as.numeric(v), hex)
encode_optional <- function(v) lapply(as.numeric(v), function(x) if (is.na(x)) NULL else hex(x))
matrix_input <- function(m) if (is.null(m)) NULL else matrix(decode(m$values), nrow=m$rows, ncol=m$columns)
compare <- function(a,b,bound) {
  stopifnot(length(a)==length(b), identical(unname(is.na(a)),unname(is.na(b))))
  finite <- is.finite(a)&is.finite(b)
  stopifnot(all(finite|is.na(a)), all(finite|is.na(b)))
  error <- if (any(finite)) max(abs(a[finite]-b[finite])) else 0
  stopifnot(error<=bound)
  error
}
blocked <- function(expr) {
  functions <- list(base=c('sqrt','exp','expm1','log','log1p','log2','mean','rank','cumsum','^','%*%','lgamma','digamma','set.seed','RNGkind'),
    stats=c('pnorm','dnorm','dnbinom','dbinom','dt','integrate','quantile','optim','optimHess','uniroot','qnorm','rnorm','runif'),
    apeglm=c('apeglm','priorVar','svalue'))
  if(requireNamespace('ashr',quietly=TRUE))functions$ashr<-c('ash','ash.workhorse','mixEM','w_mixEM','mixVBEM','cxxMixSquarem','mixSQP','estimate_mixprop','cdf_post')
  if(requireNamespace('SQUAREM',quietly=TRUE))functions$SQUAREM<-'squarem'
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name, where=asNamespace(package))), add=TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
      tracer=substitute(stop(MESSAGE,call.=FALSE),list(MESSAGE=paste('native regression executed R math:',package,name))))))
  force(expr)
}
guard <- tryCatch(blocked(stats::dnorm(0)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed R math',conditionMessage(guard)))
