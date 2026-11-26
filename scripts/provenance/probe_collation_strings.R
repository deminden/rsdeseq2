#!/usr/bin/env Rscript
# Black-box probe: R's order of random multi-character strings under ICU root
# collation, the verification oracle for rcompat-collate (strings mixing
# scripts, accents, case, compatibility variants, punctuation and digits).
#
# Only R's own string comparison is consulted. The strings are drawn with a
# fixed seed from the characters listed in root_single.tsv
# (probe_collation_full.R); nothing here is used to derive the model.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation_strings.R [dir]
# Reads root_single.tsv and writes root_strings.tsv in /tmp/rcollate by
# default. Each output line: code points of one string, relation to the next
# ("<" primary, "<<" secondary, "<<<" tertiary difference, "=" tie).
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
dir <- if (length(args)) args[[1]] else "/tmp/rcollate"
oracle <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
  "probe_collation_oracle.R")
lines <- readLines(file.path(dir, "root_single.tsv"))
end <- grep("^# completely ignorable", lines)
body <- lines[seq_len(end - 1)]
body <- body[!startsWith(body, "#")]
f <- strsplit(body, "\t", fixed = TRUE)
cp <- strtoi(vapply(f, `[[`, "", 1L), 16L)
rel <- vapply(f, `[[`, "", 2L)
pign <- vapply(f, `[[`, "", 3L) == "1"
ign <- strtoi(sub("\t.*", "", lines[(end + 1):length(lines)]), 16L)
# The head: everything before the trailing implicit run.
n <- length(cp)
inc <- c(cp[-1] > cp[-n] & rel[-n] == "<", FALSE)
head_n <- max(which(!inc[seq_len(n - 3L)]))
keep <- seq_len(head_n)
keep <- keep[!cp[keep] %in% c(0x09, 0x0A)]
marks <- cp[keep][pign[keep]]
letters_ <- keep[!pign[keep]]
gid <- cumsum(c(TRUE, rel[letters_[-length(letters_)]] == "<"))
groups <- split(cp[letters_], gid)
G <- length(groups)
set.seed(20261126)
pick <- function(x) x[sample.int(length(x), 1L)]
member <- function(g) pick(groups[[g]])
# Words over a window of neighbouring groups (one script, roughly), each
# letter replaced by a random member of its group: accents, case and
# compatibility variants of the same letters.
window_word <- function(len, width = 60L) {
  start <- sample.int(G - width, 1L)
  gs <- start + sample.int(width, len, replace = TRUE) - 1L
  gs
}
variants <- function(gs, k) vapply(seq_len(k), function(i) {
  x <- vapply(gs, member, 1L)
  if (runif(1) < 0.3) {
    at <- sample.int(length(x), 1L)
    x <- append(x, pick(marks), after = at)
  }
  intToUtf8(x)
}, "")
out <- character()
out <- c(out, unlist(lapply(1:14000, function(i) variants(window_word(sample(1:5, 1L)), 5L))))
out <- c(out, vapply(1:40000, function(i) intToUtf8(vapply(window_word(sample(1:6, 1L)), member, 1L)), ""))
ascii <- c(32:126)
out <- c(out, vapply(1:30000, function(i) intToUtf8(sample(ascii, sample(1:8, 1L), replace = TRUE)), ""))
anychar <- cp[keep]
out <- c(out, vapply(1:30000, function(i) intToUtf8(sample(c(anychar, marks, ign), sample(1:4, 1L), replace = TRUE)), ""))
latin <- c(97:122, 65:90)
out <- c(out, vapply(1:25000, function(i) {
  b <- sample(c(latin, 0x438, 0x418, 0x627, 0x648, 0x64A, 0x915, 0x995, 0xB95, 0xE01, 0x1000, 0x6C, 0x4C, 0x3B1), 1L)
  m <- sample(c(marks, 0x306, 0x653, 0x654, 0x655, 0x9BE, 0x9D7, 0xBBE, 0xBD7, 0x9C7, 0xBC6, 0xB7, 0x387, ign), sample(1:3, 1L), replace = TRUE)
  intToUtf8(c(b, m, sample(latin, sample(0:2, 1L), replace = TRUE)))
}, ""))
cjk <- c(0x4E00:0x4E80, 0x3400:0x3420, 0x20000:0x20030, 0x9F00:0x9F20, 0xF900:0xF910, 0x2F800:0x2F810,
  0x3041:0x3096, 0x30A1:0x30FA, 0xFF66:0xFF9D, 0xAC00:0xAC40, 0x1100:0x1112, 0x3131:0x3163, 0x2F00:0x2F20, 0x31C0:0x31E3)
out <- c(out, vapply(1:25000, function(i) intToUtf8(sample(c(cjk, latin, 0x30FC, 0x309B, 0x309C, 0x3099, 0x309A), sample(1:5, 1L), replace = TRUE)), ""))
# Digits, punctuation and spaces around words.
punct <- c(32, 33:47, 58:64, 91:96, 123:126, 0xA0, 0x2010:0x2027, 0x3000:0x3003, 0x660:0x669, 0x966:0x96F, 0xFF10:0xFF19)
out <- c(out, vapply(1:30000, function(i) intToUtf8(sample(c(punct, latin, 48:57), sample(1:7, 1L), replace = TRUE)), ""))
out <- unique(out[nzchar(out)])
hex <- vapply(out, function(s) paste(sprintf("%04X", utf8ToInt(s)), collapse = " "), "", USE.NAMES = FALSE)
q <- file.path(dir, "root_strings_input.tsv")
writeLines(hex, q)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(oracle), "sort", "root", shQuote(q),
  shQuote(file.path(dir, "root_strings.tsv"))))
stopifnot(status == 0)
cat("strings", length(out), "\n")
