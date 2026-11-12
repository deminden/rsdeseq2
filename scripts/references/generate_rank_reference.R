#!/usr/bin/env Rscript
# Freeze independent rank decisions; never derive expected ranks from rsdeseq2.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, getRversion() == "4.6.2",
          as.character(packageVersion("DESeq2")) == "1.52.0")
x <- seq_len(8) - 4.5
z <- rep(c(-1, 1), 4)
cases <- list(
  independent = cbind(1, x),
  duplicate = cbind(x, x),
  intercept_alias = cbind(1, rep(2, 8)),
  zero_column = cbind(1, rep(0, 8)),
  tiny_covariate = cbind(1, x * 1e-12),
  large_covariate = cbind(1, x * 1e12),
  extreme_finite_units = cbind(1e-200, x * 1e200),
  mixed_units_three_columns = cbind(1, x * 1e-12, z * 1e12),
  near_intercept_1e6 = cbind(1, 1 + x * 1e-6),
  near_intercept_1e9 = cbind(1, 1 + x * 1e-9),
  near_intercept_1e12 = cbind(1, 1 + x * 1e-12),
  near_covariate_1e6 = cbind(1, x, x + z * 1e-6),
  near_covariate_1e9 = cbind(1, x, x + z * 1e-9),
  near_covariate_1e12 = cbind(1, x, x + z * 1e-12)
)
rows <- lapply(names(cases), function(name) {
  m <- cases[[name]]
  scaled <- sweep(m, 2L, apply(abs(m), 2L, max), function(a, b) a / ifelse(b == 0, 1, b))
  data.frame(case = name, n_rows = nrow(m), n_cols = ncol(m),
    r_default_rank = qr(m)$rank, relative_rank = qr(scaled, tol = 1e-10)$rank,
    row_major_values = paste(sprintf("%.17g", as.vector(t(m))), collapse = ","))
})
write.table(do.call(rbind, rows), args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
