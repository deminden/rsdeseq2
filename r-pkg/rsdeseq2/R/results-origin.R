# Formula identifiers are interface metadata; all row eligibility runs in Rust.
.rsdeseq2_formula_origin <- function(formula) {
  list(kind="formula",variables=all.vars(formula))
}

.rsdeseq2_result_design_origin <- function(fit) {
  version <- fit$schemaVersion
  if(!is.null(version) && (!is.numeric(version)||length(version)!=1L||is.na(version)||
      !version %in% c(1,2))) stop("Unsupported fitted schema version; refit",call.=FALSE)
  origin <- fit$designOrigin
  if(is.null(origin)) {
    if(!is.null(version)&&version==2) stop("Fitted formula origin is missing; refit",call.=FALSE)
    return(list(kind="legacy",variables=character()))
  }
  if(!is.list(origin)||!setequal(names(origin),c("kind","variables"))||
      anyDuplicated(names(origin))||!is.character(origin$kind)||length(origin$kind)!=1L||
      is.na(origin$kind)||!origin$kind %in% c("formula","matrix","legacy")||
      !is.character(origin$variables)||anyNA(origin$variables)||
      any(!nzchar(origin$variables))||anyDuplicated(origin$variables)||
      (origin$kind!="formula"&&length(origin$variables))) {
    stop("Fitted design origin is invalid; refit",call.=FALSE)
  }
  if(!is.null(version)&&version==2&&origin$kind=="legacy") {
    stop("New fitted states require explicit formula or matrix origin",call.=FALSE)
  }
  origin
}

.rsdeseq2_cooks_factor_metadata <- function(colData, origin) {
  if(is.null(colData)) return(list())
  fields <- if(identical(origin$kind,"formula")) origin$variables else names(colData)
  factors <- lapply(fields,function(name) {
    value <- colData[[name]]
    if(is.factor(value)) list(codes=as.integer(value),levels=nlevels(value)) else NULL
  })
  names(factors) <- fields
  Filter(Negate(is.null),factors)
}
