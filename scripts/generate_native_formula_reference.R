#!/usr/bin/env Rscript
# Independent finite-fitting matrices/offsets from original R, never rsdeseq2.
stopifnot(getRversion()=="4.6.2")
args <- commandArgs(TRUE)
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
paths <- if(length(args))args else file.path(root,"crates/rsdeseq2/tests/data/native_formula_r462.json")
if(any(file.exists(paths)))stop("Use fresh paths; references are never overwritten")
hex <- function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
hexes <- function(v)lapply(unname(v),hex)
d <- data.frame(
  x=c(-10.1,-1.1,-.1,-0,0,.1,1.1,10.1,.3,-.3,1e-10,-1e-10),
  y=c(2.1,-.7,.3,1,-2,3,5,-7,11,13,17,19),
  p=c(.1,.3,.7,1,1.1,2,3,7,10,19,101,1e-10),
  z=c(-0,0,-.1,.1,-.3,.3,-.5,.5,1,2,7,1e-10),
  large=rep(c(2147483648,-2147483648,1e308),4),
  g=factor(rep(c("A","B"),6),levels=c("A","B")),
  h=factor(rep(c("X","X","Y","Y"),3),levels=c("X","Y")))
d$same <- d$x
d$many <- factor(rep(c("L001","L064","L129"),4),levels=sprintf("L%03d",seq_len(129)))
numeric <- lapply(names(d)[vapply(d,is.numeric,TRUE)],function(name)list(name=name,values=hexes(d[[name]])))
factors <- lapply(c("g","h","many"),function(name)list(name=name,sample_levels=as.list(as.character(d[[name]])),levels=as.list(levels(d[[name]])),reference=levels(d[[name]])[[1L]]))
cases <- list()
record <- function(name,formula,native_names=NULL,r_columns=NULL) {
  f <- as.formula(formula)
  value <- tryCatch(suppressWarnings({
    frame <- stats::model.frame(f,d,na.action=stats::na.pass)
    matrix <- stats::model.matrix(f,frame)
    offsets <- stats::model.offset(frame)
    if(is.null(r_columns))r_columns <- colnames(matrix)
    matrix <- matrix[,r_columns,drop=FALSE]
    list(matrix=matrix,offsets=offsets)
  }),error=identity)
  valid <- !inherits(value,"error")&&all(is.finite(value$matrix))&&all(is.finite(value$offsets))
  cases[[length(cases)+1L]] <<- list(name=name,formula=formula,
    native_names=if(valid)as.list(native_names)else NULL,
    r_columns=if(valid)as.list(r_columns)else NULL,
    expected=if(valid)hexes(as.vector(t(value$matrix)))else NULL,
    offsets=if(valid)hexes(value$offsets)else NULL)
}
for(power in c(2L,3L,4L,5L,7L,8L,11L,16L))
  record(paste0("power_",power),paste0("~I(x^",power,")"),c("Intercept",paste0("x_pow_",power)))
for(sign in c("+","-"))record(paste0("unary_",if(sign=="+")"plus"else"minus"),paste0("~I(",sign,"x)"),c("Intercept",if(sign=="+")"x_identity"else"x_neg"))
for(item in list(c("x+2.5","x_plus_2p5"),c("x-2.5","x_minus_2p5"),c("2.5-x","x_rminus_2p5"),c("x*2.5","x_times_2p5"),c("x/2.5","x_div_2p5"),c("2.5/y","y_rdiv_2p5"),c("x+y","x_plus_y"),c("x-y","x_minus_y"),c("x*y","x_times_y"),c("x/y","x_div_y")))
  record(item[[2L]],paste0("~I(",item[[1L]],")"),c("Intercept",item[[2L]]))
for(item in list(c("log(p)","p_log"),c("log(p,base=2)","p_log_base_2"),c("log(p,10)","p_log_base_10"),c("log(p,3)","p_log_base_3"),c("log2(p)","p_log2"),c("log10(p)","p_log10"),c("sqrt(p)","p_sqrt"),c("log1p(z)","z_log1p"),c("as.integer(x)","x_as.integer"),c("as.double(x)","x_as.double"),c("as.numeric(g)","g_as.numeric"),c("as.integer(g)","g_as.integer")))
  record(item[[2L]],paste0("~",item[[1L]]),c("Intercept",item[[2L]]))
record("raw_poly","~poly(x,4,raw=TRUE)",c("Intercept",paste0("x_poly_",1:4)))
record("orthogonal_poly","~poly(x,3)",c("Intercept",paste0("poly(x, 3)",1:3)))
record("scale","~scale(p)",c("Intercept","p_scale"))
record("factor_numeric","~g:x",c("Intercept","gA:x","gB:x"))
record("factor_numeric_no_intercept","~0+g:x",c("gA:x","gB:x"))
for(intercept in c(TRUE,FALSE)) {
  columns <- c("gA:hX","gA:hY","gB:hX","gB:hY")
  record(paste0("factor_product_",intercept),paste0("~",if(intercept)""else"0+","g:h"),c(if(intercept)"Intercept",columns),c(if(intercept)"(Intercept)",columns))
}
columns <- c("gA:hX:x","gA:hY:x","gB:hX:x","gB:hY:x")
record("factor_higher_product","~g:h:x",c("Intercept",columns),c("(Intercept)",columns))
record("numeric_higher_product","~x:y:p",c("Intercept","x:y:p"))
for(item in list(c("single_offset","~x+offset(x)"),c("different_offsets","~x+offset(log(p))+offset(I(x-y))"),c("duplicate_offset","~x+offset(x)+offset(x)"),c("duplicate_offset_backticks","~x+offset( x )+offset(`x`)"),c("duplicate_offset_number","~x+offset(log(p,base=2))+offset(log(p,base=2.0))"),c("different_equal_value_offsets","~x+offset(x)+offset(same)")))
  record(item[[1L]],item[[2L]],c("Intercept","x"))
record("invalid_integer","~as.integer(large)")
record("invalid_power","~I(large^3)")
record("invalid_division","~I(x/0)")
record("invalid_log","~log(x)")
record("invalid_offset","~x+offset(log(x))")
record("many_treatment","~many",c("Intercept",paste0("many_",levels(d$many)[-1L],"_vs_L001")))
record("many_middle_reference","~relevel(many,ref='L064')",c("Intercept",paste0('relevel(many, ref = "L064")_',setdiff(levels(d$many),"L064"),"_vs_L064")))
output <- list(metadata=list(R=as.character(getRversion()),platform=R.version$platform,BLAS=unname(extSoftVersion()["BLAS"]),blas_core=Sys.getenv("OPENBLAS_CORETYPE"),policy="Finite fitting covariates/offsets; original R matrix columns reordered only by explicit names; no fitting/rank threshold changed"),numeric=numeric,factors=factors,native_formula_cases=cases)
for(path in paths)jsonlite::write_json(output,path,auto_unbox=TRUE,null="null",pretty=TRUE,digits=NA)
cat(length(cases),"independent native formula cases generated\n")
