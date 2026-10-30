# Test-only transport of independently evaluated formula variables.
.model_reference_input <- function(case) {
  variables <- lapply(case$variables,function(v) {
    if(v$type=="numeric")list(name=v$name,kind=0L,
      values=matrix(.row_reference_values(v$values),case$samples,v$columns,byrow=TRUE),
      columnNames=if(is.null(v$column_names))NULL else as.character(unlist(v$column_names)))
    else list(name=v$name,kind=1L,values=as.integer(unlist(v$codes)),levels=as.character(unlist(v$levels)),
      coding=if(v$coding=="matrix")5L else match(v$coding,c("contr.treatment","contr.sum","contr.helmert","contr.poly","contr.SAS"))-1L,
      contrasts=if(is.null(v$contrast_values))NULL else matrix(.row_reference_values(v$contrast_values),length(v$levels),v$contrast_columns,byrow=TRUE),
      columnNames=if(is.null(v$contrast_names))NULL else as.character(unlist(v$contrast_names)))
  })
  list(as.integer(case$samples),variables,
    matrix(as.integer(unlist(case$terms)),length(variables),length(case$terms)),as.logical(case$intercept),as.integer(case$response))
}
.model_reference_frame <- function(case) {
  input <- .model_reference_input(case)
  frame <- lapply(seq_along(input[[2]]),function(i) {
    x<-input[[2]][[i]];v<-case$variables[[i]]
    if(x$kind==0L){value<-x$values;colnames(value)<-x$columnNames;return(value)}
    value<-factor(x$values,levels=seq_along(x$levels),labels=x$levels)
    attr(value,"contrasts")<-if(x$coding==5L){m<-x$contrasts;colnames(m)<-x$columnNames;m}else v$coding
    value
  })
  names(frame)<-vapply(case$variables,`[[`,character(1),"raw_name")
  codes<-input[[3]];rownames(codes)<-vapply(case$variables,`[[`,character(1),"name")
  structure(frame,class="data.frame",row.names=as.character(seq_len(case$samples)),
    terms=structure(list(),class="terms",factors=codes,intercept=as.integer(case$intercept),response=as.integer(case$response)))
}
