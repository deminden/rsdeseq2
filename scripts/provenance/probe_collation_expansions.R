#!/usr/bin/env Rscript
# Black-box probe: which characters R sorts as a sequence of several primary
# weights (expansions such as "æ" = a e), and which sequence.
#
# Input: root_single.tsv written by probe_collation_full.R. Only R's string
# comparison is consulted. A primary group (characters that are
# primary-equal) is an expansion when an earlier atom is a primary prefix of
# it: g <=p c <p g + U+FFFF. Its decomposition is found atom by atom with a
# binary search over the atoms in sorted order.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation_expansions.R [dir]
# Reads and writes in /tmp/rcollate by default (root_expansions.tsv).
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
dir <- if (length(args)) args[[1]] else "/tmp/rcollate"
icuSetCollate(locale = "root")
stopifnot(identical(icuGetCollate("actual"), "root"))
cmp <- function(a, b) ifelse(a < b, -1L, ifelse(b < a, 1L, 0L))
LO <- "\t"; HI <- "\n"; TOP <- "￿"
P <- function(x, y) cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
Pcmp <- function(x, y) ifelse(P(x, y), 0L, cmp(x, y))

lines <- readLines(file.path(dir, "root_single.tsv"))
end <- grep("^# completely ignorable", lines)
rows <- strsplit(lines[seq_len(end - 1)][!startsWith(lines[seq_len(end - 1)], "#")], "\t", fixed = TRUE)
cp <- strtoi(vapply(rows, `[[`, "", 1L), 16L)
rel <- vapply(rows, function(r) if (length(r) >= 2) r[[2]] else "", "")
pign <- vapply(rows, function(r) r[[3]] == "1", TRUE)
# The trailing run of code points in increasing order (unassigned code points
# and private use) is ordered by code point and not examined; U+FFFD and
# U+FFFF, which R puts after it, are excluded from that test.
n <- length(cp)
inc <- c(cp[-1] > cp[-n] & rel[-n] == "<", FALSE)
head_n <- max(which(!inc[seq_len(n - 3L)]))
start <- c(TRUE, rel[-n] == "<")
reps <- which(start & !pign & seq_len(n) <= head_n)
ch <- intToUtf8(cp[reps], multiple = TRUE)
# Everything between an atom q and an expansion that starts with q also starts
# with q, so one pass in sorted order suffices: a group is an expansion
# exactly when the first atom of the previous group is a primary prefix of it.
prefix_of <- function(q, c) Pcmp(q, c) <= 0 & Pcmp(c, paste0(q, TOP)) < 0 & !P(q, c)
atom <- rep(TRUE, length(reps))
first <- seq_along(reps)
for (i in seq_along(reps)[-1]) {
  q <- first[i - 1L]
  if (prefix_of(ch[[q]], ch[[i]])) { atom[i] <- FALSE; first[i] <- q }
}
exp_idx <- which(!atom)
A <- ch[atom]
cat("groups", length(reps), "atoms", length(A), "expansions", length(exp_idx), "\n")
# Decompose: W grows by the largest atom a with W + a <=p c.
target <- ch[exp_idx]
W <- rep("", length(target)); parts <- vector("list", length(target)); done <- rep(FALSE, length(target))
for (step in 1:24) {
  act <- which(!done)
  if (!length(act)) break
  lo <- rep(1L, length(act)); hi <- rep(length(A), length(act))
  while (any(lo < hi)) {
    mid <- (lo + hi + 1L) %/% 2L
    ok <- Pcmp(paste0(W[act], A[mid]), target[act]) <= 0
    lo <- ifelse(ok, mid, lo); hi <- ifelse(ok, hi, mid - 1L)
  }
  for (j in seq_along(act)) parts[[act[j]]] <- c(parts[[act[j]]], lo[j])
  W[act] <- paste0(W[act], A[lo])
  done[act] <- P(W[act], target[act])
}
stopifnot(all(done))
atom_cp <- cp[reps][atom]
header <- sprintf("# %s on %s; ICU %s; collation %s", R.version.string, R.version$platform,
  extSoftVersion()[["ICU"]], icuGetCollate("actual"))
writeLines(c(header, "# expansion (first character of its primary group)\tprimary-equal sequence of atoms",
  paste(sprintf("%04X", cp[reps][exp_idx]),
    vapply(parts, function(p) paste(sprintf("%04X", atom_cp[p]), collapse = " "), ""), sep = "\t")),
  file.path(dir, "root_expansions.tsv"), useBytes = TRUE)
