#!/usr/bin/env Rscript
# Black-box observations of the EM iterations that mixsqp runs before SQP.
# Only the exported mixsqp() interface is called; no mixsqp source is
# consulted. With convtol.sqp = 1e300 and zero.threshold.solution = 0 the
# first SQP iteration accepts its starting point, so the returned (normalised)
# x is the state after numiter.em + 1 EM updates.
#
# Candidate rule checked against R, with d = 1e-15 and normalised weights w:
#   z <- x + d;  y <- z * t(L) %*% (w / (L %*% z));  x <- y + d * max(y)
# The eps control is varied to show it does not enter the EM update.
#
# Usage: Rscript --vanilla probe_mixsqp_em_floor.R > OUTPUT.txt
stopifnot(getRversion() == "4.6.2")
suppressMessages(library(mixsqp))
cat("mixsqp", as.character(packageVersion("mixsqp")), "\n")
fmt <- function(v) paste(sprintf("%.17g", v), collapse = " ")
run <- function(L, x0, em, w, eps = 1e-8) {
  r <- withCallingHandlers(
    mixsqp(L, w, x0, control = list(
      numiter.em = em, maxiter.sqp = 1, convtol.sqp = 1e300,
      zero.threshold.solution = 0, verbose = FALSE, normalize.rows = FALSE,
      eps = eps)),
    warning = function(w) invokeRestart("muffleWarning"))
  r$x
}
model <- function(L, x0, steps, w, d = 1e-15) {
  w <- w / sum(w); x <- x0 / sum(x0)
  for (s in seq_len(steps)) {
    z <- x + d
    y <- z * drop(crossprod(L, w / drop(L %*% z)))
    x <- y + d * max(y)
  }
  x / sum(x)
}
worst <- 0
report <- function(tag, L, x0, em, w, eps = 1e-8) {
  r <- run(L, x0, em, w, eps)
  m <- model(L, x0, em + 1, w)
  err <- max(abs(r - m) / pmax(abs(m), 1e-300))
  worst <<- max(worst, err)
  cat(tag, "em", em, "eps", eps, "\n  R     ", fmt(r), "\n  model ", fmt(m),
      "\n  max relative difference", sprintf("%.3g", err), "\n")
}
# One row: the EM map of L = (1, 1) is the identity, so only the floor moves x.
for (t in c(1e-20, 1e-15, 1e-10)) for (em in 0:2)
  report(paste("L=(1,1) x0=(t,1-t) t", t), matrix(c(1, 1), 1), c(t, 1 - t), em, 1)
for (Lr in list(c(2, 1), c(1, 2), c(3, 0.5), c(1, 2, 1), c(1, 1, 2), c(2, 1, 1, 1)))
  for (em in 0:1) {
    K <- length(Lr)
    x0 <- c(1e-20, 0.3, rep(0.7 / max(K - 2, 1), max(K - 2, 1)))[seq_len(K)]
    report(paste("L=(", paste(Lr, collapse = ","), ")"), matrix(Lr, 1), x0, em, 1)
  }
# Several weighted rows, with and without identity pseudo-rows (prior counts).
set.seed(3)
Ld <- matrix(runif(15), 5, 3); a <- runif(5) * 3
for (b in list(c(0, 0, 0), c(2, 0, 0), c(1, 0.5, 0)))
  for (x0 in list(c(1e-20, 0.3, 0.7), c(0.2, 0.3, 0.5))) for (em in 0:2) {
    keep <- c(b, a) > 0
    L <- rbind(diag(3), Ld)[keep, , drop = FALSE]
    for (eps in c(0, 1e-6))
      report(paste("rows b", paste(b, collapse = ","), "x0", paste(x0, collapse = ",")),
             L, x0, em, c(b, a)[keep], eps)
  }
cat("worst relative difference", sprintf("%.3g", worst), "\n")
