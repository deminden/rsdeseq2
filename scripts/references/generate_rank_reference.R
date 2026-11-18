# R qr() rank decisions for near-collinear designs (hardening 6.2).
# Usage: Rscript --vanilla generate_rank_reference.R output.json
# Each case stores the design (column-major, exact %a hex doubles), R's
# qr(x)$rank with the default tol = 1e-7 (LINPACK dqrdc2), and whether
# DESeq2's checkFullRank stops for it.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages(library(DESeq2))
hex <- function(x) sprintf("%a", as.double(x))
cases <- list()
add <- function(name, x) {
  rank <- qr(x)$rank
  rejected <- tryCatch({
    DESeq2:::checkFullRank(x)
    FALSE
  }, error = function(e) TRUE)
  cases[[length(cases) + 1L]] <<- list(
    name = name, rows = nrow(x), cols = ncol(x), values = hex(as.vector(x)),
    rank = rank, deseq2_rejects = rejected
  )
}
set.seed(20261118)
n <- 12
group <- rep(c(0, 1), each = n / 2)
age <- c(31, 44, 52, 38, 61, 47, 35, 58, 42, 50, 39, 55)
noise <- rnorm(n)
for (eps in c(1e-2, 1e-4, 1e-6, 1e-7, 3e-8, 1e-8, 1e-9, 1e-11, 0)) {
  add(sprintf("age_copy_eps_%g", eps), cbind(1, group, age, age + eps * noise))
  add(sprintf("group_copy_eps_%g", eps), cbind(1, group, group + eps * noise))
}
for (scale in c(1e-8, 1e8)) {
  add(sprintf("scaled_columns_%g", scale), cbind(1, group, age * scale))
}
add("exact_sum", cbind(1, group, 1 - group))
add("zero_column", cbind(1, group, 0))
batch <- rep(c(0, 1, 0), length.out = n)
add("interaction_full", cbind(1, group, batch, group * batch))
add("nearly_saturated", diag(n)[, 1:(n - 1)] + cbind(1e-9 * noise, matrix(0, n, n - 2)))
out <- list(
  generator = "scripts/references/generate_rank_reference.R",
  r_version = R.version.string,
  deseq2_version = as.character(packageVersion("DESeq2")),
  tolerance = 1e-7,
  cases = cases
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE), args[[1]])
