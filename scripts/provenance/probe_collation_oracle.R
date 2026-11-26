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
# the level tests below are valid for the collation, otherwise "<". "compare"
# writes -1, 0 or 1 per pair. The first output line records the collation ICU
# reported and the level test used.
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
# Level tests (see probe_collation.R). LO = tab is the lowest primary and
# HI = line feed: x and y have equal primaries exactly when x+LO+HI > y+LO+LO
# and x+LO+LO < y+LO+HI. With equal primaries they have equal secondaries
# when a letter with an accent placed after them decides (x+LO+acc > y+LO+"a"
# and y+LO+acc > x+LO+"a"); under a backwards secondary level (fr_CA) the
# accented letter goes before them instead (acc+LO+x against "a"+LO+y). The
# tab keeps the letter out of contractions with the end of x (Danish "aa");
# acc is the first of a-acute, a-grave, a-breve and a-macron that the
# collation sorts as a secondary variant of "a" (in Icelandic a-acute is a
# letter of its own). The checks below choose the placement or, when neither
# is valid, report no levels.
LO <- "\t"; HI <- "\n"
P <- function(x, y) cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
accs <- c("á", "à", "ă", "ā")
acc <- accs[P("a", accs) & cmp("a", accs) < 0][1]
mkS <- function(before) function(x, y) {
  if (before) {
    P(x, y) & cmp(paste0(acc, LO, x), paste0("a", LO, y)) > 0 & cmp(paste0(acc, LO, y), paste0("a", LO, x)) > 0
  } else {
    P(x, y) & cmp(paste0(x, LO, acc), paste0(y, LO, "a")) > 0 & cmp(paste0(y, LO, acc), paste0(x, LO, "a")) > 0
  }
}
valid <- function(S) S("a", "A") && !S("a", acc) && !S("ab", paste0(acc, "b")) && !S(paste0("a", acc), paste0(acc, "a"))
placement <- "none"
if (!is.na(acc)) {
  for (before in c(FALSE, TRUE)) {
    S <- mkS(before)
    if (valid(S)) {
      placement <- if (before) "before" else "after"
      break
    }
  }
}
levels_valid <- placement != "none"
level <- function(x, y) ifelse(!P(x, y), 1L, ifelse(!S(x, y), 2L, 3L))
test <- if (levels_valid) sprintf("level test U+%04X %s", utf8ToInt(acc), placement) else "no level test"
header <- sprintf("# %s on %s; ICU %s; requested %s; collation %s; %s", R.version.string,
  R.version$platform, extSoftVersion()[["ICU"]], locale, actual, test)
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
    i <- which(c0 < 0)
    p <- P(x[i], x[i + 1L])
    s <- p
    s[p] <- S(x[i][p], x[i + 1L][p])
    rel[i] <- ifelse(!p, "<", ifelse(!s, "<<", "<<<"))
  }
  writeLines(c(header, paste(hex(x), c(rel, ""), sep = "\t")), args[[4]], useBytes = TRUE)
}
