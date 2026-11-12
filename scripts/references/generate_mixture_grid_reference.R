#!/usr/bin/env Rscript
# Execute original grid preparation, replacing only its downstream optimizer.
# No candidate calculations or fitted-weight claims enter this grid oracle.
main <- function(args) {
if(length(args)!=1L) stop('usage: generate_mixture_grid_reference.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',
          as.character(packageVersion('ashr'))=='2.2.63',!file.exists(args[[1L]]))
ns <- asNamespace('ashr')
original <- get('estimate_mixprop',ns)
captured <- NULL
capture <- function(data,g,prior,...) {
  captured <<- list(g=g,prior=prior)
  list(g=g,penloglik=NA_real_)
}
unlockBinding('estimate_mixprop',ns);assign('estimate_mixprop',capture,ns);lockBinding('estimate_mixprop',ns)
on.exit({ unlockBinding('estimate_mixprop',ns);assign('estimate_mixprop',original,ns);lockBinding('estimate_mixprop',ns) },add=TRUE)
encode <- function(x) as.list(sprintf('%.17g',as.numeric(x)))
bits <- function(x) lapply(as.numeric(x),function(v) paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
x <- c(-3,-1,-.2,0,.2,1,4,NA_real_)
s <- c(.1,.3,.7,1,2,1,.4,NA_real_)
profiles <- list(
  automatic=list(pointmass=FALSE,prior='uniform'),
  automatic_null=list(pointmass=TRUE,prior='nullbiased',nullweight=7),
  shifted=list(pointmass=FALSE,prior='uniform',mixsd=c(.1,.7,2),mode=.3,grange=c(-2,3)),
  clipped=list(pointmass=TRUE,prior='nullbiased',mixsd=c(.1,.7,2),grange=c(-.25,.25)))
families <- c('normal','uniform','halfuniform','+uniform','-uniform','halfnormal')
cases <- list()
for(profile in names(profiles)) for(family in families) {
  p <- profiles[[profile]];warnings <- character()
  result <- withCallingHandlers(do.call(ashr::ash,c(list(betahat=x,sebetahat=s,
    mixcompdist=family,outputlevel=0),p)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  g <- captured$g
  cases[[length(cases)+1L]] <- list(name=paste(profile,family,sep='_'),family=family,
    estimates=bits(x),standard_errors=bits(s),scales=if(is.null(p$mixsd))NULL else bits(p$mixsd),
    point_mass=p$pointmass,mode=bits(if(is.null(p$mode))0 else p$mode)[[1L]],
    range=bits(if(is.null(p$grange))c(-Inf,Inf)else p$grange),
    nullweight=if(p$prior=='nullbiased')bits(if(is.null(p$nullweight))10 else p$nullweight)[[1L]]else NULL,
    component_class=class(g)[[1L]],initial_weights=encode(g$pi),dirichlet=encode(captured$prior),
    means=if(is.null(g$mean))NULL else encode(g$mean),sds=if(is.null(g$sd))NULL else encode(g$sd),
    lower=if(is.null(g$a))NULL else encode(g$a),upper=if(is.null(g$b))NULL else encode(g$b),
    warnings=as.list(warnings))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),ashr=as.character(packageVersion('ashr')),
  scope='original preparation; downstream optimizer replaced by capture'),cases=cases),
  args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')

}
main(commandArgs(TRUE))
