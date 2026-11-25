#!/usr/bin/env Rscript
# Independent references for the special functions and the gamma distribution
# (docs/provenance/specs/special_functions.md): lgamma, gamma, digamma,
# pgamma (both tails), qgamma, beta, lbeta, dnorm and qnorm on grids covering
# small, moderate and large arguments, both tails, negative non-integers and
# the edge values. Runs base R only; never loads rsdeseq2.
stopifnot(getRversion() == "4.6.2")
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed = TRUE)), "..", ".."))
path <- if (length(args)) args[[1L]] else
  file.path(root, "crates/rcompat-numeric/tests/data/special_functions_r462.tsv")
if (file.exists(path)) stop("Use a fresh path; references are never overwritten")
hex <- function(v) paste(format(writeBin(as.double(v), raw(), size = 8, endian = "big")), collapse = "")
num <- function(v) sprintf("%.17g", v)
rows <- list()
add <- function(fn, a, b, value) {
  rows[[length(rows) + 1L]] <<- data.frame(
    fn = fn, a = num(a), b = if (is.na(b[1L]) && !is.nan(b[1L])) "-" else num(b),
    value = vapply(value, hex, ""), stringsAsFactors = FALSE)
}
quiet <- function(expr) suppressWarnings(expr)

negatives <- function(shifts) unlist(lapply(shifts, function(k) -seq(0.05, 0.95, by = 0.1) - k))
near_one_two <- c(1 - 2^-seq(1, 52, by = 3), 1 + 2^-seq(1, 52, by = 3),
                  2 - 2^-seq(1, 52, by = 3), 2 + 2^-seq(1, 52, by = 3))
edges <- c(0, -0, -1, -2, -10, Inf, -Inf, NaN)

x <- c(2^-1074, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, seq(0.01, 0.09, by = 0.01),
       seq(0.1, 12, by = 0.1), near_one_two, 12.5, 15, 20, 31.7, 50, 99.9, 171.5,
       500.25, 1e3, 1e4, 1e5, 1e6, 1e10, 1e15, 1e20, 1e50, 1e100, 1e200, 1e300, 1e306,
       negatives(c(0, 1, 2, 3, 4, 5, 9, 20)), -50.3, -100.7, -150.2, -170.4, -170.6,
       -171.5, -200.25, -1000.75, -12345.6789, -1e6 - 0.5, -1e10 - 0.25, -1e15 - 0.5,
       edges, -1e6, -2^52)
add("lgamma", x, NA, quiet(lgamma(x)))

# Range limits of gamma, located by bisection and written in hexadecimal so
# they parse exactly: smallest finite positive argument and its predecessor,
# overflow at about 171.6145, zero below about -170.567, and tiny negatives.
limits <- c(0x1.0292a5d2f227ap-1022, 0x1.0292a5d2f2279p-1022, -0x1.0292a5d2f227ap-1022,
            -0x1.0000000000004p-1024, -0x1p-1024, 0x1.573a9cf985201p+7,
            0x1.573a9cf985202p+7, -0x1.55228f00a54c5p+7, -0x1.55228f00a54c6p+7)
x <- c(limits, 1e-300, 1e-100, 1e-10, 1e-3, seq(0.05, 10, by = 0.1), near_one_two,
       seq(10.5, 171.5, by = 3.5), 1:30, 50, 100, 170, 172, 1e3,
       negatives(c(0, 1, 2, 5, 10, 50, 100, 169)), -170.05, -170.5, -170.55,
       -171.5, -200.5, edges, -50)
add("gamma", x, NA, quiet(gamma(x)))

root_psi <- 1.4616321449683622
x <- c(1e-300, 1e-100, 1e-10, 1e-5, 1e-3, seq(0.05, 12, by = 0.1),
       root_psi + c(-0.5, -0.2, -0.05, -0.01, -1e-3, -1e-5, 1e-5, 1e-3, 0.01, 0.05, 0.2, 0.5),
       15, 20, 50, 100, 1e3, 1e5, 1e10, 1e15, 1e50, 1e300,
       negatives(c(0, 1, 2, 5, 20, 100)), -1000.3, -1e6 - 0.5, edges)
add("digamma", x, NA, quiet(digamma(x)))

shapes <- c(1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 19.5,
            20, 30, 50, 100, 300, 1e3, 1e4, 1e5, 1e6)
probs <- c(1e-300, 1e-100, 1e-30, 1e-10, 1e-4, 0.01, 0.1, 0.3, 0.5)
for (a in shapes) {
  x <- c(quiet(qgamma(probs, a)), quiet(qgamma(probs, a, lower.tail = FALSE)),
         a * c(0.5, 0.9, 1, 1.1, 2), a + 1, a + sqrt(a), max(a - sqrt(a), a / 4))
  x <- sort(unique(x[is.finite(x) & x > 0]))
  add("pgamma_lower", rep(a, length(x)), x, pgamma(x, a))
  add("pgamma_upper", rep(a, length(x)), x, pgamma(x, a, lower.tail = FALSE))
}
a <- c(0.5, 1, 3, 1e3)
add("pgamma_lower", c(a, a), c(rep(0, 4), 1e300 + 0 * a), pgamma(c(rep(0, 4), 1e300 + 0 * a), c(a, a)))
add("pgamma_upper", c(a, a), c(rep(0, 4), 1e300 + 0 * a),
    pgamma(c(rep(0, 4), 1e300 + 0 * a), c(a, a), lower.tail = FALSE))

quantile_shapes <- c(0.01, 0.1, 0.5, 1, 2.5, 7, 20, 50, 250, 1e4)
quantile_probs <- c(1e-300, 1e-50, 1e-10, 1e-3, 0.05, 0.25, 0.5, 0.75, 0.95,
                    1 - 1e-3, 1 - 1e-10, 1 - 2^-52)
for (rate in c(1, 0.5, 3.75)) for (a in quantile_shapes) {
  add(paste0("qgamma_rate_", rate), rep(a, length(quantile_probs)), quantile_probs,
      qgamma(quantile_probs, a, rate))
}
add("qgamma_rate_1", rep(2, 2), c(0, 1), qgamma(c(0, 1), 2))

beta_args <- c(1e-300, 1e-10, 0.001, 0.0625, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3.25,
               5, 7.75, 10, 12.5, 20, 33.5, 50, 100.25, 1e3, 1e4 + 0.5, 1e6, 1e10, 1e15)
pairs <- expand.grid(a = beta_args, b = beta_args)
pairs <- pairs[pairs$a <= pairs$b, ]
pairs <- rbind(pairs, data.frame(a = c(1, 2.5, Inf), b = c(Inf, Inf, Inf)))
add("lbeta", pairs$a, pairs$b, lbeta(pairs$a, pairs$b))
add("beta", pairs$a, pairs$b, beta(pairs$a, pairs$b))

x <- c(0, 1e-300, 1e-8, seq(-40, 40, by = 0.37), 37.5, 38.5, -38.6, 1e10, Inf, -Inf, NaN)
add("dnorm", x, NA, dnorm(x))

p <- c(0, 1, NaN, 5e-324, 1e-310, 10^-(1:300), 1 - 10^-(1:16), seq(0.01, 0.99, by = 0.01),
       0.075, 0.07499999999999999, 0.075000000000000011, 0.925, 0.92499999999999993,
       exp(-25), exp(-25) * (1 - 2^-50), exp(-25) * (1 + 2^-50), 1 - 2^-53)
add("qnorm", p, NA, qnorm(p))

out <- do.call(rbind, rows)
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("#", R.version.string, "on", R.version$platform),
             "function\ta\tb\tvalue_bits"),
           path)
write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE, append = TRUE)
