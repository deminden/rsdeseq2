#!/usr/bin/env Rscript
# Black-box observations of the mixsqp search direction when the reduced
# Hessian of the quadratic subproblem is numerically singular. Only exported
# interfaces of ashr and mixsqp are called; no package source is consulted.
#
# The likelihood matrix is the one ash() builds (outputlevel = 3) for an
# eight-gene input whose standard errors are large next to the default
# normal grid. mixsqp is run for a single SQP iteration with
# maxiter.activeset = 1, ..., 8, so the returned x shows the active-set
# trajectory of the first subproblem (after the line search). Each is
# compared with two candidate rules for the reduced solve, both followed by
# the same Armijo backtracking line search:
#   min-norm: drop eigen-directions with eigenvalue <= n * eps * max|lambda|
#   ridge:    solve (A + r I) p = -g for a few ridges r
#
# Usage: Rscript --vanilla probe_mixsqp_singular_direction.R > OUTPUT.txt
stopifnot(getRversion() == "4.6.2")
suppressMessages({library(ashr); library(mixsqp)})
cat("ashr", as.character(packageVersion("ashr")), "mixsqp",
    as.character(packageVersion("mixsqp")), "\n")
b <- c(0.8938459130692038, -0.30807406620808614, 0.5188649759446969,
       0.6867711762458827, -1.0225709949451243, 0.7285431855079626,
       0.78147041493471, 0.4468318106955514)
s <- c(4.047762889349763, 4.133025421302511, 2.9889064674996098,
       2.839866842552506, 3.4874172938544956, 2.566068424736632,
       4.09496696651059, 2.579278191601985)
fit <- ash(b, s, mixcompdist = "normal", method = "shrink", outputlevel = 3)
cat("ash loglik", sprintf("%.17g", fit$loglik), "pi", fit$fitted_g$pi, "\n")
L <- fit$fit_details$matrix_lik; n <- nrow(L); k <- ncol(L)
w <- rep(1 / n, n); e <- 1e-6
# Starting point: 20 + 1 EM updates with the floor recorded by
# probe_mixsqp_em_floor.R, from the uniform proportions.
x <- rep(1 / k, k)
for (i in 1:21) {
  z <- x + 1e-15; y <- z * drop(crossprod(L, w / drop(L %*% z))); x <- y + 1e-15 * max(y)
}
f <- function(x) sum(x) - sum(w * log(drop(L %*% x) + e))
u <- drop(L %*% x) + e
g <- 1 - drop(crossprod(L, w / u))
H <- crossprod(L, (w / u^2) * L)
cat("objective at start", sprintf("%.10g", f(x) - sum(x)), "(verbose mixsqp prints 1.497916890e-02)\n")
cat("eigenvalues of H", signif(eigen(H, symmetric = TRUE)$values, 4), "\n")
solve_min_norm <- function(A, r) {
  E <- eigen(A, symmetric = TRUE)
  keep <- E$values > nrow(A) * .Machine$double.eps * max(abs(E$values))
  drop(E$vectors[, keep, drop = FALSE] %*% (crossprod(E$vectors[, keep, drop = FALSE], r) / E$values[keep]))
}
solve_ridge <- function(r0) function(A, r) solve(A + r0 * diag(nrow(A)), r)
qp <- function(solver, iters) {
  y <- x; c <- g - drop(H %*% x); free <- y > 0
  for (it in seq_len(iters)) {
    gr <- drop(H %*% y + c); idx <- which(free); p <- rep(0, k)
    if (length(idx)) p[idx] <- solver(H[idx, idx, drop = FALSE], -gr[idx])
    if (max(abs(p)) <= 1e-14) {
      j <- which(!free); if (!length(j)) break
      j <- j[which.min(gr[j])]; if (gr[j] >= -1e-10) break
      free[j] <- TRUE; next
    }
    st <- 1; bl <- NA
    for (i in idx) if (p[i] < 0 && -y[i] / p[i] < st) { st <- -y[i] / p[i]; bl <- i }
    y <- pmax(y + st * p, 0); if (!is.na(bl)) { y[bl] <- 0; free[bl] <- FALSE }
  }
  d <- y - x; slope <- sum(g * d); t <- 1
  while (f(pmax(x + t * d, 0)) > f(x) + 0.01 * t * slope) t <- t * 0.75
  x <- pmax(x + t * d, 0)
  x / sum(x)  # mixsqp returns normalised proportions
}
fmt <- function(v) paste(sprintf("%.6g", v), collapse = " ")
for (iters in 1:8) {
  r <- suppressWarnings(mixsqp(L, x0 = rep(1 / k, k), control = list(
    eps = e, numiter.em = 20, maxiter.sqp = 1, maxiter.activeset = iters,
    verbose = FALSE)))$x
  cat("maxiter.activeset", iters, "\n  mixsqp       ", fmt(r), "\n")
  m <- qp(solve_min_norm, iters)
  cat("  min-norm     ", fmt(m), "  max abs diff", sprintf("%.2g", max(abs(m - r))), "\n")
  for (r0 in c(1e-14, 1e-12, 1e-10)) {
    m <- qp(solve_ridge(r0), iters)
    cat(sprintf("  ridge %-7.0e", r0), fmt(m), "  max abs diff", sprintf("%.2g", max(abs(m - r))), "\n")
  }
}
