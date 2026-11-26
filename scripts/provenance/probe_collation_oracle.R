#!/usr/bin/env Rscript
# Black-box oracle: ask R to sort or compare strings under one ICU collation.
#
# Only R's own string comparison is consulted (`<`, `order`, `icuSetCollate`,
# `icuGetCollate`). The derivation script build_collation_table.py and the
# verification examples of rcompat-collate write the questions; this script
# writes R's answers.
#
# Usage:
#   Rscript --vanilla probe_collation_oracle.R sort    <locale> <input> <output>
#   Rscript --vanilla probe_collation_oracle.R compare <locale> <input> <output>
# Input lines hold strings as space-separated hexadecimal code points; for
# "compare" two strings separated by a tab. Lines starting with "#" are
# skipped. "sort" writes the distinct strings in R's order (ties broken by
# code points) with the relation to the next string: "<" primary, "<<"
# secondary, "<<<" tertiary difference, "=" tie; the level is reported when
# the level tests below are valid for the collation (they are not under a
# backwards secondary level, as in fr_CA), otherwise "<". "compare" writes
# -1, 0 or 1 per pair. The first output line records the collation ICU
# reported.
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L, args[[1]] %in% c("sort", "compare"))
mode <- args[[1]]
locale <- args[[2]]
icuSetCollate(locale = locale)
actual <- icuGetCollate("actual")
lines <- readLines(args[[3]])
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
dec <- function(h) intToUtf8(strtoi(strsplit(h, " ", fixed = TRUE)[[1]], 16L))
decode <- function(x) vapply(x, dec, "", USE.NAMES = FALSE)
hex <- function(x) vapply(x, function(s) paste(sprintf("%04X", utf8ToInt(s)), collapse = " "), "", USE.NAMES = FALSE)
cmp <- function(a, b) ifelse(a < b, -1L, ifelse(b < a, 1L, 0L))
LO <- "\t"; HI <- "\n"
# Level tests (see probe_collation.R): primary equality, then secondary.
P <- function(x, y) cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
S <- function(x, y) P(x, y) & cmp(paste0(x, "a\u0301"), paste0(y, "a")) > 0 & cmp(paste0(y, "a\u0301"), paste0(x, "a")) > 0
level <- function(x, y) ifelse(!P(x, y), 1L, ifelse(!S(x, y), 2L, 3L))
levels_valid <- S("a", "A") && !S("a", "\u00e1") && !S("ab", "\u00e1b") && !S("a\u00e1", "\u00e1a")
header <- sprintf("# %s on %s; ICU %s; requested %s; collation %s", R.version.string,
  R.version$platform, extSoftVersion()[["ICU"]], locale, actual)
if (mode == "compare") {
  parts <- strsplit(lines, "\t", fixed = TRUE)
  a <- decode(vapply(parts, `[[`, "", 1L))
  b <- decode(vapply(parts, `[[`, "", 2L))
  r <- cmp(a, b)
  out <- c("<", "=", ">")[r + 2L]
  if (levels_valid) {
    lv <- level(a, b)
    out[r != 0] <- strrep(out[r != 0], lv[r != 0])
  }
  writeLines(c(header, out), args[[4]])
} else {
  x <- unique(decode(lines))
  x <- x[order(x, hex(x))]
  n <- length(x)
  c0 <- cmp(x[-n], x[-1])
  stopifnot(all(c0 <= 0))
  rel <- ifelse(c0 == 0, "=", "<")
  if (levels_valid) {
    LO <- "\t"; HI <- "\n"
    P <- function(x, y) cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
    S <- function(x, y) P(x, y) & cmp(paste0(x, "á"), paste0(y, "a")) > 0 & cmp(paste0(y, "á"), paste0(x, "a")) > 0
    i <- which(c0 < 0)
    p <- P(x[i], x[i + 1L])
    s <- p
    s[p] <- S(x[i][p], x[i + 1L][p])
    rel[i] <- ifelse(!p, "<", ifelse(!s, "<<", "<<<"))
  }
  writeLines(c(header, paste(hex(x), c(rel, ""), sep = "\t")), args[[4]], useBytes = TRUE)
}
