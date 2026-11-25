#!/usr/bin/env Rscript
# Black-box probe: how does qf(p, df1, df2) behave for very large degrees of
# freedom? Compares qf with the chi-square limits
#   df2 -> Inf: qf(p, df1, df2) -> qchisq(p, df1) / df1
#   df1 -> Inf: qf(p, df1, df2) -> df2 / qchisq(p, df2, lower.tail = FALSE)
# and finds, by bisection on integer degrees of freedom, where qf starts to
# return the limit exactly. Runs base R only.
stopifnot(getRversion() == "4.6.2")
p <- 0.99
bits <- function(v) paste(format(writeBin(as.double(v), raw(), size = 8, endian = "big")), collapse = "")
limit2 <- function(df1) qchisq(p, df1) / df1
limit1 <- function(df2) df2 / qchisq(p, df2, lower.tail = FALSE)
same <- function(a, b) identical(bits(a), bits(b))
first_exact <- function(is_exact, low, high) {
  stopifnot(!is_exact(low), is_exact(high))
  while (high - low > 1) {
    mid <- floor((low + high) / 2)
    if (is_exact(mid)) high <- mid else low <- mid
  }
  high
}
cat("# ", R.version.string, " on ", R.version$platform, "\n", sep = "")
cat("case\tdf1\tthreshold_df_other\tnote\n")
for (df1 in c(1, 2, 5, 20, 200)) {
  t <- first_exact(function(df2) same(qf(p, df1, df2), limit2(df1)), 1000, 1e7)
  cat("df2_limit\t", df1, "\t", t, "\tsmallest integer df2 with qf == qchisq(p, df1)/df1\n", sep = "")
}
for (df2 in c(1, 2, 5, 20, 200)) {
  t <- first_exact(function(df1) same(qf(p, df1, df2), limit1(df2)), 1000, 1e7)
  cat("df1_limit\t", df2, "\t", t, "\tsmallest integer df1 with qf == df2/qchisq(p, df2, lower.tail=FALSE)\n", sep = "")
}
# Which limit applies when both are large?
for (pair in list(c(5e5, 6e5), c(6e5, 5e5), c(5e5, 5e5), c(4e5, 4e5 + 1))) {
  q <- qf(p, pair[[1]], pair[[2]])
  cat("both\t", pair[[1]], "\t", pair[[2]], "\tqf=", bits(q),
      " df2_limit=", same(q, limit2(pair[[1]])),
      " df1_limit=", same(q, limit1(pair[[2]])), "\n", sep = "")
}
