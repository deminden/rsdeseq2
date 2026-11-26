#!/usr/bin/env Rscript
# Black-box probe: R's order of every code point under the ICU root collation
# or under the collation of a locale.
#
# Only R's own string comparison is consulted (`<`, `order`, `icuSetCollate`,
# `icuGetCollate`) plus `grepl(perl = TRUE)` to report which code points
# R's regular-expression engine considers unassigned. No collation table or
# implementation is read. The output is the verification oracle and the
# input of the table derivation (`crates/rcompat-collate/examples/derive.rs`);
# docs/provenance/specs/collation.md describes the method.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation_full.R [output-dir [locale]]
# Writes root_single.tsv (about 12 MB), or <locale>_single.tsv for another
# locale, into /tmp/rcollate by default; the SHA-256 of each output is
# recorded in docs/provenance/observations/collation_r462_macos/oracle_summary.tsv.
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1]] else "/tmp/rcollate"
locale <- if (length(args) > 1L) args[[2]] else "root"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
icuSetCollate(locale = locale)
actual <- icuGetCollate("actual")
stopifnot(locale != "root" || identical(actual, "root"))

cmp <- function(a, b) ifelse(a < b, -1L, ifelse(b < a, 1L, 0L))
LO <- "\t"; HI <- "\n"
# Primary equality (see probe_collation.R): appending LO HI to one side and
# LO LO to the other lets the primary level decide unless the primaries differ.
P <- function(x, y) cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
# Primary and secondary equality: "a" + U+0301 against "a" decides at the
# secondary level unless the secondaries already differ.
S <- function(x, y) P(x, y) & cmp(paste0(x, "á"), paste0(y, "a")) > 0 & cmp(paste0(y, "á"), paste0(x, "a")) > 0
if (locale != "root") {
  # Under a tailoring the accented letter is kept out of contractions with
  # the end of x by a tab, is chosen among letters the collation sorts as a
  # secondary variant of "a", and goes before x under a backwards secondary
  # level, as in probe_collation_oracle.R.
  accs <- c("\u00e1", "\u00e0", "\u0103", "\u0101")
  acc <- accs[P("a", accs) & cmp("a", accs) < 0][1]
  before <- cmp(paste0("a", LO, acc), paste0(acc, LO, "a")) > 0
  S <- function(x, y) {
    if (before) {
      P(x, y) & cmp(paste0(acc, LO, x), paste0("a", LO, y)) > 0 & cmp(paste0(acc, LO, y), paste0("a", LO, x)) > 0
    } else {
      P(x, y) & cmp(paste0(x, LO, acc), paste0(y, LO, "a")) > 0 & cmp(paste0(y, LO, acc), paste0(x, LO, "a")) > 0
    }
  }
  stopifnot(!S("a", acc), !S("ab", paste0(acc, "b")), !S(paste0("a", acc), paste0(acc, "a")))
}
stopifnot(P("a", "A"), !P("a", "b"), S("a", "A"))
if (locale == "root") stopifnot(!S("a", "\u00e1"))

# Every code point except surrogates and NUL.
cps <- c(1:0xD7FF, 0xE000:0x10FFFF)
chars <- intToUtf8(cps, multiple = TRUE)
# R's regular-expression engine: unassigned (general category Cn).
unassigned <- grepl("^\\p{Cn}$", chars, perl = TRUE)
# Completely ignorable: appending the character to "a" changes nothing.
ignorable <- cmp(paste0("a", chars), "a") == 0
live <- which(!ignorable)
o <- live[order(chars[live], cps[live])]
x <- chars[o]
n <- length(x)
c0 <- cmp(x[-n], x[-1])
stopifnot(all(c0 <= 0))
rel <- character(n - 1)
rel[c0 == 0] <- "="
i <- which(c0 < 0)
p <- P(x[i], x[i + 1])
rel[i[!p]] <- "<"
j <- i[p]
s <- S(x[j], x[j + 1])
rel[j[!s]] <- "<<"
rel[j[s]] <- "<<<"
# Primary-ignorable: "a" + c is primary-equal to "a" (LO itself is excluded).
pign <- P(paste0("a", x), "a") & x != LO
header <- if (locale == "root") {
  sprintf("# %s on %s; ICU %s; collation %s", R.version.string, R.version$platform,
    extSoftVersion()[["ICU"]], actual)
} else {
  sprintf("# %s on %s; ICU %s; requested %s; collation %s", R.version.string, R.version$platform,
    extSoftVersion()[["ICU"]], locale, actual)
}
writeLines(c(header,
  "# non-ignorable code points in R's order: code point, relation to the next (< primary, << secondary, <<< tertiary, = tie), primary-ignorable (1/0), unassigned per PCRE (1/0)",
  paste(sprintf("%04X", cps[o]), c(rel, ""), as.integer(pign), as.integer(unassigned[o]), sep = "\t"),
  "# completely ignorable: code point, unassigned per PCRE",
  paste(sprintf("%04X", cps[ignorable]), as.integer(unassigned[ignorable]), sep = "\t")),
  file.path(out, sprintf("%s_single.tsv", gsub("[@=]", "_", locale))), useBytes = TRUE)
cat("code points", length(cps), "ignorable", sum(ignorable), "unassigned(PCRE)", sum(unassigned),
  "primary-ignorable", sum(pign), "relations", paste(names(table(rel)), table(rel), collapse = " "), "\n")
