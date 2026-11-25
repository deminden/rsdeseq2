#!/usr/bin/env Rscript
# Black-box observations of the compiled unbounded L-BFGS behind the apeglm
# compiled methods, called through the entry point documented in the
# RcppNumerical README, optim_lbfgs(f, x, fx_opt, maxit, eps_f, eps_g).
# Each objective is an R closure; every evaluation point is recorded so that
# trial steps, stopping rules, status codes and warnings can be read off.
# Only inputs and outputs are observed; no optimizer source is consulted.
# Output: docs/provenance/observations/unbounded_lbfgs_probe_r462_macos.txt
# usage: Rscript probe_unbounded_lbfgs.R OUTPUT.txt
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else stdout()
suppressMessages(Rcpp::sourceCpp(code = '
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
class Callback: public MFuncGrad {
public:
  Rcpp::Function fn;
  Callback(Rcpp::Function fn): fn(fn) {}
  double f_grad(Constvec& x, Refvec grad) {
    Rcpp::NumericVector xv(x.data(), x.data() + x.size());
    Rcpp::List r = fn(xv);
    Rcpp::NumericVector g = r[1];
    for (int i = 0; i < x.size(); i++) grad[i] = g[i];
    return Rcpp::as<double>(r[0]);
  }
};
// [[Rcpp::export]]
Rcpp::List probe(Rcpp::Function fn, Rcpp::NumericVector start, int maxit, double eps_f, double eps_g) {
  Eigen::VectorXd x = Rcpp::as<Eigen::VectorXd>(start);
  double fopt = NA_REAL;
  Callback f(fn);
  int status = optim_lbfgs(f, x, fopt, maxit, eps_f, eps_g);
  return Rcpp::List::create(Rcpp::Named("x") = x, Rcpp::Named("f") = fopt,
                            Rcpp::Named("status") = status);
}
'))
lines <- sprintf("# R %s, RcppNumerical %s, %s", getRversion(), packageVersion("RcppNumerical"), R.version$platform)
emit <- function(...) lines <<- c(lines, paste0(...))
run <- function(name, fg, start, maxit = 300, eps_f = 1e-8, eps_g = 1e-8) {
  log <- list()
  wrapped <- function(x) { r <- fg(x); log[[length(log) + 1L]] <<- c(x, r[[1L]], r[[2L]]); r }
  warned <- character()
  res <- tryCatch(withCallingHandlers(probe(wrapped, start, maxit, eps_f, eps_g),
                                      warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") }),
                  error = function(e) list(error = conditionMessage(e)))
  on.exit(for (w in warned) emit("warning: ", w), add = TRUE)
  emit("## ", name, " start=", paste(sprintf("%.17g", start), collapse = ","),
       " maxit=", maxit, " eps_f=", eps_f, " eps_g=", eps_g)
  p <- length(start)
  for (i in seq_along(log)) {
    v <- log[[i]]
    emit(sprintf("eval %d x=%s f=%.17g g=%s", i,
                 paste(sprintf("%.17g", v[1:p]), collapse = ","), v[p + 1L],
                 paste(sprintf("%.17g", v[(p + 2L):(2L * p + 1L)]), collapse = ",")))
  }
  if (!is.null(res$error)) emit("result error: ", res$error) else
    emit(sprintf("result status=%d f=%.17g x=%s", res$status, res$f,
                 paste(sprintf("%.17g", res$x), collapse = ",")))
}
rosen <- function(x) {
  t1 <- x[2] - x[1]^2; t2 <- 1 - x[1]
  list(100 * t1^2 + t2^2, c(-400 * x[1] * t1 - 2 * t2, 200 * t1))
}
quad <- function(d) function(x) list(0.5 * sum(d * x^2), d * x)
quartic <- function(x) list(sum(x^4) + 1, 4 * x^3)
expo <- function(x) list(exp(x[1]) - 2 * x[1] + 10, exp(x[1]) - 2)
nanwall <- function(x) if (x[1] > 3) list(NaN, NaN) else list((x[1] - 2)^2 + 10, 2 * (x[1] - 2))
logcosh <- function(x) list(sum(log(cosh(x))) + 10, tanh(x))
run("rosenbrock_default_tolerances", rosen, c(-1.2, 1), eps_f = 1e-6, eps_g = 1e-5)
run("rosenbrock", rosen, c(-1.2, 1))
run("quadratic_1d", quad(4), 3)
run("quadratic_2d", quad(c(1, 100)), c(1, 1))
run("quadratic_2d_offset10", function(x) { r <- quad(c(1, 100))(x); r[[1]] <- r[[1]] + 10; r }, c(1, 1))
run("quartic_1d", quartic, 1)
run("exp_1d", expo, 0)
run("exp_1d_far", expo, -5)
run("nan_beyond_3", nanwall, 0)
run("nan_beyond_3_close", nanwall, 1.9)
run("logcosh_2d", logcosh, c(2, -1))
run("logcosh_2d_tight", logcosh, c(2, -1), eps_f = 1e-12, eps_g = 1e-10)
run("linear_unbounded", function(x) list(-x[1], -1), 0, maxit = 5)
run("rosenbrock_maxit3", rosen, c(-1.2, 1), maxit = 3)
run("ill_conditioned", quad(c(1, 1e6)), c(1, 1))
# Failure modes: tolerances that cannot be met in double precision.
run("offset_quadratic_no_tolerance", function(x) list(1e3 + (x[1] - 1)^2, 2 * (x[1] - 1)), 0.3, eps_f = 0, eps_g = 0)
run("large_offset_logcosh_no_tolerance", function(x) list(1e8 + sum(log(cosh(x))), tanh(x)), c(1, -0.5), eps_f = 0, eps_g = 0)
run("logcosh_offset10_no_tolerance", function(x) list(10 + sum(log(cosh(x))), tanh(x)), c(2, -1), eps_f = 0, eps_g = 0)
run("rosenbrock_offset10_tight", function(x) { r <- rosen(x); r[[1]] <- r[[1]] + 10; r }, c(-1.2, 1), eps_f = 1e-12, eps_g = 1e-10)
run("inconsistent_gradient", function(x) list((x[1] - 1)^2 + 10, -2 * (x[1] - 1)), 0, eps_f = 0, eps_g = 0)
run("nan_near_start", function(x) if (x[1] > 0.5) list(NaN, NaN) else list((x[1] - 2)^2 + 10, 2 * (x[1] - 2)), 0)
run("inf_after_unit_step", function(x) if (x[1] > 0.5) list(Inf, Inf) else list((x[1] - 2)^2 + 10, 2 * (x[1] - 2)), 0)
# Relative objective test: which value scales the decrease.
run("logcosh_1d_eps_f_0.05", function(x) list(log(cosh(x[1])), tanh(x[1])), 1.5, eps_f = 0.05, eps_g = 0)
run("logcosh_1d_plus1_eps_f_0.05", function(x) list(1 + log(cosh(x[1])), tanh(x[1])), 1.5, eps_f = 0.05, eps_g = 0)
run("logcosh_1d_eps_f_0.5", function(x) list(log(cosh(x[1])), tanh(x[1])), 1.5, eps_f = 0.5, eps_g = 0)
run("logcosh_1d_eps_g_0.1", function(x) list(log(cosh(x[1])), tanh(x[1])), 1.5, eps_f = 0, eps_g = 0.1)
run("logcosh_1d_from3_eps_f_0.6", function(x) list(log(cosh(x[1])), tanh(x[1])), 3, eps_f = 0.6, eps_g = 0)
run("logcosh_1d_from3_eps_f_0.7", function(x) list(log(cosh(x[1])), tanh(x[1])), 3, eps_f = 0.7, eps_g = 0)
run("logcosh_1d_from3_eps_f_0.9", function(x) list(log(cosh(x[1])), tanh(x[1])), 3, eps_f = 0.9, eps_g = 0)
# One-dimensional quadratics k (x - m)^2 started at 0: the unit first step
# lands on x = 1, so m fixes the interpolated or extrapolated trial.
for (m in c(0.002, 0.01, 0.03, 0.05, 0.08, 0.09, 0.095, 0.1, 0.105, 0.11, 0.2, 0.45, 0.51, 0.52, 0.53, 0.9, 0.98, 0.995, 9.9, 10.1))
  run(sprintf("quadratic_minimum_at_%g", m), function(x) list(3 * (x[1] - m)^2 + 10, 6 * (x[1] - m)), 0, eps_f = 0, eps_g = 1e-6)
for (m in c(0.0105, 0.012, 0.014, 0.0155))
  run(sprintf("quadratic_minimum_at_%g", m), function(x) list(3 * (x[1] - m)^2 + 10, 6 * (x[1] - m)), 0, eps_f = 0, eps_g = 1e-6)
# Zoom started from an acceptable but too steep unit step: f(0) = f0,
# g(0) = -1, f(1) = 10, g(1) = 0.95; elsewhere a smooth bowl at 0.99.
for (f0 in c(12, 30, 60, 110, 1000))
  run(sprintf("steep_return_f0_%g", f0), function(x) {
    if (x[1] == 0) list(f0, -1) else if (x[1] == 1) list(10, 0.95) else list(9.5 + (x[1] - 0.99)^2, 2 * (x[1] - 0.99))
  }, 0, maxit = 1, eps_f = 0, eps_g = 0)
# Precision-limited problems with no attainable tolerance.
run("ill_conditioned_offset10_no_tolerance", function(x) { r <- quad(c(1, 1e6))(x); r[[1]] <- r[[1]] + 10; r }, c(1, 1), eps_f = 0, eps_g = 0)
run("rosenbrock_offset10_no_tolerance", function(x) { r <- rosen(x); r[[1]] <- r[[1]] + 10; r }, c(-1.2, 1), eps_f = 0, eps_g = 0)
run("exp_1d_no_tolerance", expo, 0, eps_f = 0, eps_g = 0)
run("quartic_offset10_no_tolerance", function(x) list(10 + sum((x - c(1, -2))^4), 4 * (x - c(1, -2))^3), c(0, 0), eps_f = 0, eps_g = 0)
writeLines(lines, out)
