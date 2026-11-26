#!/usr/bin/env Rscript
# Black-box probe: R's sorted order of chosen and random strings under ICU
# root collation, and which locale identifiers ICU reports as using the root
# collation. The sorted lists check the model derived by probe_collation.R;
# nothing here is used to build the weight table.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation_orders.R [output-dir]
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1]] else "docs/provenance/observations/collation_r462_macos"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
header <- sprintf("# %s on %s; ICU %s", R.version.string, R.version$platform, extSoftVersion()[["ICU"]])

# 1. Locale identifiers: the collation ICU reports as actually used, and
# whether a fixed sample sorts as under root.
sample_words <- c("a", "A", "\u00e5", "\u00e4", "\u00f6", "o", "z", "aa", "ch", "c", "h", "ll", "l",
  "\u00f1", "n", "\u00e7", "\u0131", "i", "I", "\u0130", "\u00fe", "th", "y", "\u00fc", "u", "v", "w",
  "\u0161", "s", "\u017e", "\u0142", "\u00df", "ss", "\u00e6", "\u00f8", "cote", "cot\u00e9", "c\u00f4te",
  "c\u00f4t\u00e9", "\u0451", "\u0435", "\u0439", "\u0438", "\u03ac", "\u03b1", "1", "-", " ")
icuSetCollate(locale = "root")
root_order <- order(sample_words)
locales <- c("root", "", "en", "en_US", "en_GB", "en_AU", "en_CA", "en_IE", "en_NZ", "en_ZA", "en_IN",
  "en_US_POSIX", "de", "de_DE", "de_AT", "de_CH", "fr", "fr_FR", "fr_CA", "fr_BE", "fr_CH", "it", "it_IT",
  "nl", "nl_NL", "pt", "pt_BR", "pt_PT", "es", "es_ES", "es_MX", "ca", "gl", "eu", "ga", "cy", "id", "ms",
  "sw", "af", "sv", "da", "nb", "nn", "fi", "is", "et", "lv", "lt", "pl", "cs", "sk", "sl", "hr", "bs",
  "sr", "hu", "ro", "tr", "az", "sq", "mt", "ru", "uk", "be", "bg", "mk", "el", "he", "ar", "fa", "hi",
  "th", "vi", "ja", "zh", "ko")
loc <- do.call(rbind, lapply(locales, function(l) {
  icuSetCollate(locale = l)
  data.frame(requested = l, actual = icuGetCollate("actual"),
    root_order = identical(order(sample_words), root_order))
}))
icuSetCollate(locale = "root")
writeLines(c(header, "# requested\tactual\tsample sorts as under root",
  paste(loc$requested, loc$actual, loc$root_order, sep = "\t")), file.path(out, "locales.tsv"))

# 2. Sorted sets. Each line: code points of one string, then the relation to
# the next string in R's order ("<" or "=" for a tie).
hex <- function(x) vapply(x, function(s) paste(sprintf("%04X", utf8ToInt(s)), collapse = " "), "")
write_sorted <- function(x, file) {
  x <- unique(x)
  x <- x[order(x, hex(x))]
  n <- length(x)
  rel <- c(ifelse(x[-n] < x[-1], "<", ifelse(x[-1] < x[-n], ">", "=")), "")
  stopifnot(!any(rel == ">"))
  writeLines(c(header, "# code points\trelation to the next string",
    paste(hex(x), rel, sep = "\t")), file.path(out, file), useBytes = TRUE)
}
set.seed(20261126)
rstr <- function(n, alphabet, maxlen) vapply(seq_len(n), function(i)
  paste(sample(alphabet, sample.int(maxlen, 1L), replace = TRUE), collapse = ""), "")
ascii <- intToUtf8(32:126, multiple = TRUE)
letters_digits <- c(letters, LETTERS, 0:9)
write_sorted(c(rstr(3000, ascii, 6), rstr(2000, c(letters_digits, " ", "_", "-", "."), 8),
  rstr(1000, c("a", "A", "b", "B", "1", "2", "_", " ", "."), 5)), "sorted_ascii.tsv")
ranges <- list(c(0x20, 0x7E), c(0xA0, 0x24F), c(0x250, 0x36F), c(0x370, 0x52F), c(0x1E00, 0x1EFF),
  c(0x2000, 0x22FF), c(0xFB00, 0xFB06), c(0xFF01, 0xFF5E))
inventory <- intToUtf8(unlist(lapply(ranges, function(r) r[1]:r[2])), multiple = TRUE)
write_sorted(rstr(6000, inventory, 4), "sorted_inventory.tsv")
# Words: base letters with accents, case and combining marks, close together.
latin <- c("a", "e", "i", "o", "u", "c", "n", "s", "z", "l", "d", "t")
accented <- c("\u00e1", "\u00e0", "\u00e2", "\u00e4", "\u00e9", "\u00e8", "\u00ea", "\u00eb", "\u00ed",
  "\u00f3", "\u00f6", "\u00f8", "\u00fa", "\u00fc", "\u00e7", "\u00f1", "\u0161", "\u017e", "\u0142",
  "\u0111", "\u00f0", "\u00fe", "\u00df", "\u00e6", "\u0153", "\u0140", "\u017f")
marks <- c("\u0300", "\u0301", "\u0302", "\u0303", "\u0308", "\u030a", "\u030c", "\u0327", "\u0328", "\u0335")
greek <- intToUtf8(c(0x3AC:0x3CE), multiple = TRUE)
cyr <- intToUtf8(c(0x430:0x45F), multiple = TRUE)
word <- function(n, base, extra) vapply(seq_len(n), function(i) {
  k <- sample(2:6, 1L)
  w <- sample(c(base, extra), k, replace = TRUE)
  up <- runif(k) < 0.2; w[up] <- toupper(w[up])
  m <- runif(k) < 0.15; w[m] <- paste0(w[m], sample(marks, sum(m), replace = TRUE))
  paste(w, collapse = "")
}, "")
write_sorted(c(word(3000, latin, accented), word(1000, greek, character()), word(1000, cyr, c("\u0438\u0306", "\u0306")),
  word(500, c("l", "L"), c("\u00b7", "\u0387", "j", "a"))), "sorted_words.tsv")
# Other scripts: Hebrew, Arabic, Devanagari, kana, emoji and conjoining jamo
# are in the probed ranges, Hangul syllables are checked through their jamo,
# and CJK ideographs are outside the probed ranges (fallback ordering).
other <- intToUtf8(c(0x5D0:0x5EA, 0x627:0x64A, 0x905:0x939, 0x3041:0x3096, 0x4E00:0x4E40, 0xAC00:0xAC40,
  0x1F600:0x1F620), multiple = TRUE)
write_sorted(c(rstr(1500, c(other, letters, "1", " "), 3)), "sorted_other_scripts.tsv")
# Strings that differ only at the secondary and tertiary levels: characters
# that share a first primary weight in collation_table.tsv (written by
# probe_collation.R) are swapped for one another, so R's order is decided by
# accents, case and variants. The table only chooses the inputs; the
# expected order is R's.
tab <- read.delim(file.path(out, "collation_table.tsv"), comment.char = "#", header = FALSE,
  colClasses = "character")
tab <- tab[tab$V2 != "" & strtoi(tab$V1, 16L) < 0x1100, ]
inv <- intToUtf8(strtoi(tab$V1, 16L), multiple = TRUE)
first <- sub("\\..*", "", tab$V2)
by_primary <- split(inv, first)
by_primary <- by_primary[lengths(by_primary) > 1 & names(by_primary) != "0"]
variants <- vapply(seq_len(6000), function(i) {
  ps <- sample(names(by_primary), sample.int(4L, 1L), replace = TRUE)
  paste(vapply(ps, function(p) sample(by_primary[[p]], 1L), ""), collapse = "")
}, "")
skeletons <- vapply(seq_len(1500), function(i)
  paste(sample(c("a", "s", "l", "d", "o", "e", "i", "n"), sample(2:4, 1L), TRUE), collapse = ""), "")
respelled <- unlist(lapply(skeletons, function(w) {
  ch <- strsplit(w, "")[[1]]
  vapply(1:4, function(j) paste(vapply(ch, function(c) sample(by_primary[[first[match(c, inv)]]], 1L), ""),
    collapse = ""), "")
}))
combining <- inv[first == "0"]
with_marks <- vapply(seq_len(4000), function(i) paste(sample(c(sample(inv, 3L), sample(combining, 2L)),
  sample(2:5, 1L), TRUE), collapse = ""), "")
write_sorted(c(variants, respelled, with_marks), "sorted_variants.tsv")
cat("locales", nrow(loc), "root-like", sum(loc$actual == "root"), "\n")
