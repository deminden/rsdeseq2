#!/usr/bin/env Rscript
# Independent original-R draws and serialized state, not candidate output.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',!file.exists(args[[1L]]))
hex<-function(v,size)paste(format(writeBin(v,raw(),size=size,endian='big')),collapse='')
state<-function(v)lapply(v,function(x)hex(as.integer(x),4L))
bits<-function(v)lapply(v,function(x)hex(as.double(x),8L))
cases<-list()
record<-function(name,seed=NULL,kind=NULL){
 before<-.Random.seed
 u<-runif(64);after_u<-.Random.seed
 assign('.Random.seed',before,.GlobalEnv)
 n<-rnorm(64,0,.5);after_n<-.Random.seed
 cases[[length(cases)+1L]]<<-list(name=name,seed=seed,kind=kind,state=state(before),uniforms=bits(u),normals=bits(n),
   uniform_state=state(after_u),normal_state=state(after_n))
}
for(kind in c('Mersenne-Twister',"L'Ecuyer-CMRG"))for(seed in c(0L,1L,42L,-42L))for(warm in c(0L,17L,620L)){
 RNGkind(kind,'Inversion','Rejection');set.seed(seed)
 if(warm>0L)invisible(runif(warm))
 record(paste(kind,seed,warm,sep='_'),if(warm==0L)seed else NULL,kind)
}
for(seed in c(-917011752L,-1990828124L)) {
 RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(seed)
 record(paste0('CMRG_seed_range_boundary_',seed),seed,"L'Ecuyer-CMRG")
}
for(kind in c('Mersenne-Twister',"L'Ecuyer-CMRG")){
 suppressWarnings(RNGkind(kind,'Inversion','Rounding'));set.seed(42L)
 record(paste0(kind,'_rounding_header'))
}
RNGkind('Mersenne-Twister','Inversion','Rejection')
.Random.seed<-c(10403L,1L,1L,0L,NA_integer_,rep(1L,621L))
record('MT_zero_word_and_integer_NA_bits')
RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(42L)
.Random.seed<-parallel::nextRNGStream(.Random.seed)
record('CMRG_parallel_stream')
jsonlite::write_json(list(source=list(R=as.character(getRversion()),normal='Inversion',draws_per_profile=64L),
 random_cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
