#!/usr/bin/env Rscript
# Black-box probe of the BLAS reductions behind R's qr() on this host.
#
# R's qr() (LINPACK-style Householder QR) calls the BLAS routines ddot, dnrm2
# and daxpy of the BLAS that R is linked against. This script compiles a tiny
# C shim against that same BLAS ($(BLAS_LIBS), as an R package would) and
# observes the routines on chosen inputs. No BLAS or R source is consulted.
#
# Summation order is revealed with the method of Xie et al. (2024), "FPRev:
# Revealing the order of floating-point summation": with x[i] = 2^60,
# x[j] = -2^60 and every other term 1, the result equals n minus the number of
# leaves under the lowest common ancestor of i and j in the summation tree.
#
# Usage: Rscript --vanilla probe_reduction_order.R OUTPUT_DIR
stopifnot(getRversion() == "4.6.2")
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else "reduction-probe"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

src <- file.path(tempdir(), "probe_blas.c")
writeLines(c(
  "#include <R.h>",
  "#include <Rinternals.h>",
  "#include <R_ext/BLAS.h>",
  "SEXP probe_ddot(SEXP x, SEXP y) {",
  "  int n = LENGTH(x), one = 1;",
  "  return ScalarReal(F77_CALL(ddot)(&n, REAL(x), &one, REAL(y), &one));",
  "}",
  "SEXP probe_dnrm2(SEXP x) {",
  "  int n = LENGTH(x), one = 1;",
  "  return ScalarReal(F77_CALL(dnrm2)(&n, REAL(x), &one));",
  "}",
  "SEXP probe_daxpy(SEXP a, SEXP x, SEXP y) {",
  "  int n = LENGTH(x), one = 1;",
  "  SEXP r = PROTECT(duplicate(y));",
  "  F77_CALL(daxpy)(&n, REAL(a), REAL(x), &one, REAL(r), &one);",
  "  UNPROTECT(1);",
  "  return r;",
  "}"), src)
old <- setwd(tempdir())
Sys.setenv(PKG_LIBS = "$(BLAS_LIBS) $(FLIBS)")
stopifnot(system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "probe_blas.c"),
                  stdout = FALSE, stderr = FALSE) == 0)
setwd(old)
dyn.load(file.path(tempdir(), paste0("probe_blas", .Platform$dynlib.ext)))
ddot <- function(x, y) .Call("probe_ddot", as.double(x), as.double(y))
dnrm2 <- function(x) .Call("probe_dnrm2", as.double(x))
daxpy <- function(a, x, y) .Call("probe_daxpy", as.double(a), as.double(x), as.double(y))
bits <- function(v) paste(format(writeBin(as.double(v), raw(), size = 8, endian = "big")), collapse = "")

info <- c(R = R.version.string, platform = R.version$platform, BLAS = extSoftVersion()[["BLAS"]])
writeLines(paste(names(info), info, sep = "\t"), file.path(out, "environment.tsv"))

# --- Summation tree of ddot with y = 1 ------------------------------------
big <- 2^60
lca_sizes <- function(n) {
  m <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    x <- rep(1, n); x[i] <- big; x[j] <- -big
    m[i, j] <- m[j, i] <- as.integer(n - ddot(x, rep(1, n)))
  }
  m
}
# Rebuild the binary tree bottom-up: the smallest LCA size joins first.
tree_string <- function(m) {
  n <- nrow(m)
  if (n == 1L) return("x1")
  groups <- as.list(seq_len(n))
  label <- paste0("x", seq_len(n))
  repeat {
    if (length(groups) == 1L) return(label[[1L]])
    best <- NULL
    for (a in seq_along(groups)) for (b in seq_along(groups)) if (a < b) {
      s <- m[groups[[a]][[1L]], groups[[b]][[1L]]]
      if (s == length(groups[[a]]) + length(groups[[b]]) &&
          (is.null(best) || s < best[[3L]])) best <- list(a, b, s)
    }
    if (is.null(best)) return(NA_character_)
    a <- best[[1L]]; b <- best[[2L]]
    # Within a node the operand that holds the earlier first term is written first.
    groups[[a]] <- c(groups[[a]], groups[[b]]); groups[[b]] <- NULL
    label[[a]] <- paste0("(", label[[a]], " + ", label[[b]], ")"); label <- label[-b]
  }
}
trees <- vapply(1:80, function(n) tree_string(lca_sizes(n)), "")
writeLines(paste(1:80, trees, sep = "\t"), file.path(out, "ddot_tree.tsv"))

# --- Fused multiply-add in ddot --------------------------------------------
# a*b is inexact for these values; a fused update keeps the low bits.
a <- 1 + 2^-30; b <- 1 - 2^-30   # a*b = 1 - 2^-60, rounds to 1 in binary64
fma_rows <- character()
for (n in c(1:40, 63, 64, 65, 96, 97)) for (pos in unique(c(1L, n))) {
  x <- rep(0, n); y <- rep(0, n)
  x[pos] <- a; y[pos] <- b
  # With the subtraction of 1 elsewhere the result exposes the product's low bits.
  other <- if (pos == 1L) n else 1L
  if (n > 1L) { x[other] <- -1; y[other] <- 1 }
  r <- ddot(x, y)
  fma_rows <- c(fma_rows, paste(n, pos, bits(r), sep = "\t"))
}
writeLines(c("n\tposition\tresult_bits", fma_rows), file.path(out, "ddot_fma.tsv"))
# The same test with the -1 term placed at `previous`, the term that the
# summation tree above adds immediately before `position` in one accumulator.
lane_rows <- character()
for (n in c(32, 33, 48, 64, 65, 96, 128, 130)) for (pos in 17:n) for (gap in c(16, 32)) {
  previous <- pos - gap
  if (previous < 1L) next
  x <- rep(0, n); y <- rep(0, n)
  x[previous] <- -1; y[previous] <- 1
  x[pos] <- a; y[pos] <- b
  lane_rows <- c(lane_rows, paste(n, previous, pos, bits(ddot(x, y)), sep = "\t"))
}
writeLines(c("n\tprevious\tposition\tresult_bits", lane_rows), file.path(out, "ddot_fma_pairs.tsv"))

# --- daxpy: fused or separate multiply and add ----------------------------
ax_rows <- character()
for (n in c(1:40, 64, 65)) {
  x <- rep(a, n); y <- rep(-1, n)
  r <- daxpy(b, x, y)
  ax_rows <- c(ax_rows, paste(n, paste(vapply(r, bits, ""), collapse = ","), sep = "\t"))
}
writeLines(c("n\tresult_bits", ax_rows), file.path(out, "daxpy_fma.tsv"))

# --- dnrm2 ------------------------------------------------------------------
cases <- list(
  tiny = rep(1e-200, 4), huge = rep(1e200, 4), mixed = c(1e-200, 1, 1e200),
  ordinary = c(3, 4), increasing = 1:10, decreasing = 10:1,
  zeros = c(0, 0, 5, 0), alternating = rep(c(1, 1e-8), 8),
  random = seq(0.1, 3.3, length.out = 33))
set.seed(1)
for (k in 1:20) cases[[paste0("random", k)]] <- rnorm(sample(1:70, 1))
nrm_rows <- vapply(names(cases), function(nm) {
  x <- cases[[nm]]
  paste(nm, length(x), paste(vapply(x, bits, ""), collapse = ","), bits(dnrm2(x)), sep = "\t")
}, "")
writeLines(c("case\tn\tinput_bits\tresult_bits", nrm_rows), file.path(out, "dnrm2.tsv"))

# --- Cross-check: qr() of one column returns -sign(x1) * dnrm2(x) ----------
chk <- vapply(cases, function(x) {
  q <- qr(cbind(x))$qr[1L, 1L]
  identical(bits(q), bits(-sign(if (x[[1L]] == 0) 1 else x[[1L]]) * dnrm2(x)))
}, TRUE)
writeLines(paste(names(chk), chk, sep = "\t"), file.path(out, "qr_norm_crosscheck.tsv"))

# --- Random dot products for acceptance ------------------------------------
set.seed(2)
dot_rows <- vapply(1:300, function(k) {
  n <- sample(c(1:100, 200, 257, 1000), 1)
  x <- rnorm(n) * 2^sample(-20:20, n, TRUE); y <- rnorm(n)
  paste(paste(vapply(x, bits, ""), collapse = ","), paste(vapply(y, bits, ""), collapse = ","),
        bits(ddot(x, y)), sep = "\t")
}, "")
writeLines(c("x_bits\ty_bits\tresult_bits", dot_rows), file.path(out, "ddot_random.tsv"))
nrm_random <- vapply(1:300, function(k) {
  n <- sample(c(1:100, 200, 257), 1)
  x <- rnorm(n) * 2^sample(-20:20, n, TRUE)
  paste(paste(vapply(x, bits, ""), collapse = ","), bits(dnrm2(x)), sep = "\t")
}, "")
writeLines(c("x_bits\tresult_bits", nrm_random), file.path(out, "dnrm2_random.tsv"))
