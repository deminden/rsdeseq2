#!/usr/bin/env Rscript
# Black-box observations of R's default qr() (LAPACK = FALSE, tol = 1e-7):
# rank, pivot, compact QR matrix and qraux for chosen matrices, plus
# qr.qy(QR, diag(diag(QR$qr), n, p)) as poly() and contr.poly() bases use it.
# Base R only; no R or LINPACK source is consulted.
#
# Usage: Rscript --vanilla probe_qr_rules.R OUTPUT.json
stopifnot(getRversion() == "4.6.2")
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else "qr_cases_r462_macos.json"
bits <- function(v) vapply(as.double(v), function(x)
  paste(format(writeBin(x, raw(), size = 8, endian = "big")), collapse = ""), "")
set.seed(20261125)
cases <- list()
add <- function(name, x, tol = 1e-7) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  q <- qr(x, tol = tol)
  y <- diag(diag(q$qr), nrow(x), ncol(x))
  qy <- if (nrow(x) > 0L && ncol(x) > 0L) qr.qy(q, y) else y
  cases[[length(cases) + 1L]] <<- list(
    name = name, rows = nrow(x), cols = ncol(x), tolerance = bits(tol),
    x = bits(x), rank = q$rank, pivot = q$pivot,
    qr = bits(q$qr), qraux = bits(q$qraux), qy = bits(qy))
}
add("identity3", diag(3))
add("zero_column_middle", cbind(1:4, 0, c(2, 3, 5, 7)))
add("zero_column_first", cbind(0, 1:4, c(2, 3, 5, 7)))
add("all_zero", matrix(0, 3, 2))
add("one_row", matrix(c(2, 3, 5), 1))
add("one_column", cbind(c(3, 4)))
add("wide_full", matrix(rnorm(6), 2, 3))
add("wide_dependent", cbind(1:2, 2 * (1:2), c(1, -1)))
add("tall_dependent", cbind(1, 1:6, 2 * (1:6)))
add("duplicate_last", cbind(1, 1:5, 1:5))
for (eps in c(1e-5, 1e-6, 1e-7, 1e-8, 1e-11)) {
  base <- cbind(1, seq(0, 1, length.out = 8))
  add(paste0("near_collinear_", format(eps)), cbind(base, base[, 2] + eps * sin(1:8)))
  add(paste0("near_collinear_scaled_", format(eps)), cbind(base, 1e6 * (base[, 2] + eps * cos(1:8))))
}
for (tol in c(1e-3, 1e-12)) add(paste0("tolerance_", format(tol)), cbind(1, 1:5, (1:5) + 1e-5 * (1:5)^2), tol)
add("magnitudes", cbind(1e-8 * (1:6), 1e8 * c(1, 0, 2, 0, 3, 1), c(5, 4, 3, 2, 1, 0)))
for (n in c(3, 5, 8, 12, 20, 33, 40, 70)) {
  x <- seq(-1.3, 2.9, length.out = n)
  add(paste0("vandermonde_", n), outer(x - mean(x), 0:min(n - 1, 6), "^"))
}
for (n in c(3, 4, 7, 10)) {
  center <- (n + 1) / 2
  add(paste0("contr_poly_", n), outer(seq_len(n) - center, 0:(n - 1), "^"))
}
for (k in 1:30) {
  n <- sample(1:40, 1); p <- sample(1:8, 1)
  x <- matrix(rnorm(n * p) * 2^sample(-8:8, n * p, TRUE), n, p)
  if (p > 2 && k %% 3 == 0) x[, p] <- x[, 1] - 2 * x[, 2]
  add(paste0("random_", k), x)
}
json <- function(v) {
  if (is.list(v)) {
    if (!is.null(names(v)))
      return(paste0("{", paste0('"', names(v), '":', vapply(v, json, ""), collapse = ","), "}"))
    return(paste0("[", paste(vapply(v, json, ""), collapse = ","), "]"))
  }
  if (is.character(v)) return(paste0("[", paste0('"', v, '"', collapse = ","), "]"))
  paste0("[", paste(v, collapse = ","), "]")
}
meta <- list(R = R.version.string, platform = R.version$platform,
             BLAS = extSoftVersion()[["BLAS"]],
             note = "x, qr and qy are column-major; values are big-endian binary64 hex")
writeLines(paste0('{"metadata":', json(meta), ',"cases":', json(cases), "}"), out)
