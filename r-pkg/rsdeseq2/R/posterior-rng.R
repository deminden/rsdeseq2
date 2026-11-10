# Import stored state only. Draws and state transitions are performed in Rust.
.rsdeseq2_posterior_rng_state <- function(env = .GlobalEnv) {
  if (!exists('.Random.seed',envir=env,inherits=FALSE)) return(NULL)
  value <- .rsdeseq2_frozen_binding('.Random.seed',env,context='random state')
  if (is.object(value) || typeof(value) != 'integer' || !is.null(dim(value))) {
    stop('random state must be a stored primitive integer vector',call.=FALSE)
  }
  if (bindingIsLocked('.Random.seed',env)) stop('random state binding must be writable',call.=FALSE)
  value
}
