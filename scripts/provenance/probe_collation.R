#!/usr/bin/env Rscript
# Black-box probe: how does R order strings under ICU root collation?
#
# Only R's own string comparison is consulted (`<`, `sort`, `order`,
# `icuSetCollate`, `icuGetCollate`). No collation table or implementation is
# read. Every weight written to collation_table.tsv is derived from the
# answers R gives for strings chosen here; docs/provenance/specs/collation.md
# describes the method and the resulting model.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation.R [output-dir]
# Writes into docs/provenance/observations/collation_r462_macos by default.
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1]] else "docs/provenance/observations/collation_r462_macos"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
icuSetCollate(locale = "root")
stopifnot(identical(icuGetCollate("actual"), "root"))

cps_hex <- function(x) vapply(x, function(s) paste(sprintf("%04X", utf8ToInt(s)), collapse = " "), "")
cmp <- function(a, b) ifelse(a < b, -1L, ifelse(b < a, 1L, 0L))
LO <- "\t"; HI <- "\n"; TOP <- "\uffff"
# Primary equality. Appending LO HI to one side and LO LO to the other makes the
# primary level decide unless the primary sequences already differ. LO (tab) is
# the lowest non-ignorable character observed, so a proper primary prefix is
# never mistaken for equality (strings under test never contain LO).
P <- function(x, y) {
  cmp(paste0(x, LO, HI), paste0(y, LO, LO)) > 0 & cmp(paste0(x, LO, LO), paste0(y, LO, HI)) < 0
}
# Primary and secondary equality: with equal primaries, appending "a" + U+0301
# to one side and "a" to the other decides at the secondary level unless the
# secondary sequences already differ.
ACUTE <- "a\u0301"
S <- function(x, y) P(x, y) & cmp(paste0(x, ACUTE), paste0(y, "a")) > 0 & cmp(paste0(y, ACUTE), paste0(x, "a")) > 0
Pcmp <- function(x, y) ifelse(P(x, y), 0L, cmp(x, y))

# Checks of the assumptions the method rests on.
stopifnot(cmp(LO, HI) < 0, cmp(paste0(LO, "a"), paste0(HI, LO)) < 0, cmp("a", "b") < 0,
  P("a", "A"), !P("a", "b"), S("a", "A"), !S("a", "\u00e1"), cmp("a\u0301", "\u00e1") == 0)

ranges <- list(c(0x01, 0x24F), c(0x250, 0x36F), c(0x370, 0x52F), c(0x1E00, 0x1EFF),
  c(0x2000, 0x22FF), c(0xFB00, 0xFB06), c(0xFF01, 0xFF5E), c(0x590, 0x5FF), c(0x600, 0x6FF),
  c(0x900, 0x97F), c(0x2300, 0x23FF), c(0x2460, 0x27BF), c(0x3040, 0x30FF), c(0x1F300, 0x1F64F),
  c(0x1100, 0x11FF))
inventory <- unlist(lapply(ranges, function(r) r[1]:r[2]))
chars <- intToUtf8(inventory, multiple = TRUE)
names(chars) <- sprintf("%04X", inventory)

# 1. Completely ignorable characters: appending them changes nothing.
ignorable <- cmp(paste0("a", chars), "a") == 0
# 2. Sorted order of the remaining characters.
live <- chars[!ignorable]
live <- live[order(live)]
n <- length(live)
# Level of the first difference between neighbours: 0 tie, 1 primary,
# 2 secondary, 3 tertiary.
lev <- ifelse(cmp(live[-n], live[-1]) == 0, 0L, ifelse(!P(live[-n], live[-1]), 1L,
  ifelse(!S(live[-n], live[-1]), 2L, 3L)))
stopifnot(all(cmp(live[-n], live[-1]) <= 0))
# 3. Primary-ignorable characters (combining marks): "a" + c is primary-equal
# to "a". LO itself is the lowest primary and is excluded from this test.
pign <- P(paste0("a", live), "a") & live != LO
marks <- live[pign]
base <- live[!pign]
nb <- length(base)
gid <- cumsum(c(TRUE, !P(base[-nb], base[-1])))
groups <- split(base, gid)
reps <- vapply(groups, `[[`, "", 1L)
names(reps) <- vapply(groups, function(x) names(x)[1], "")
G <- length(groups)
stopifnot(reps[[1]] == LO)

# 4. Primary decomposition of each group: the shortest chain of groups whose
# concatenation is primary-equal to the group. A group g is a prefix of c
# when g <=p c <p g + TOP.
decomp <- vector("list", G)
decomp[[1]] <- 1L
for (g in seq_len(G)[-1]) {
  c0 <- reps[[g]]; W <- ""; parts <- integer()
  for (step in 1:6) {
    cand <- paste0(W, reps)
    ok <- Pcmp(cand, c0) <= 0 & Pcmp(c0, paste0(cand, TOP)) < 0
    ok[1] <- FALSE
    if (nchar(W)) ok <- ok & !P(cand, W)
    h <- which(ok)
    stopifnot(length(h) > 0)
    parts <- c(parts, h[1]); W <- paste0(W, reps[[h[1]]])
    if (P(W, c0)) break
  }
  stopifnot(P(W, c0))
  decomp[[g]] <- parts
}
atom <- vapply(seq_len(G), function(g) identical(decomp[[g]], g), TRUE)
atom[1] <- TRUE
flat <- function(g) if (atom[g]) g else unlist(lapply(decomp[[g]], flat))
primary <- integer(G); primary[atom] <- seq_len(sum(atom))
prim_of <- lapply(seq_len(G), function(g) primary[flat(g)])

# 5. Canonical ties: members equal to another member plus one combining mark,
# or equal to an earlier member outright.
tie <- setNames(rep("", length(base)), names(base))
for (g in seq_len(G)) {
  mem <- groups[[g]]
  if (length(mem) < 2) next
  cand <- as.vector(outer(mem, marks, paste0))
  for (k in seq_along(mem)) {
    c0 <- mem[[k]]; nm <- names(mem)[k]
    same <- mem[seq_len(k - 1)][cmp(mem[seq_len(k - 1)], c0) == 0]
    if (length(same)) { tie[nm] <- names(same)[1]; next }
    hit <- which(cmp(cand, c0) == 0)
    if (length(hit)) {
      i <- hit[1]; b <- ((i - 1) %% length(mem)) + 1; m <- ((i - 1) %/% length(mem)) + 1
      tie[nm] <- paste(names(mem)[b], names(marks)[m])
    }
  }
}

# Marks equal to a sequence of two other marks.
pairs <- outer(marks, marks, paste0)
for (k in seq_along(marks)) {
  hit <- which(cmp(pairs, marks[[k]]) == 0)
  if (length(hit)) {
    i <- hit[1]; a <- ((i - 1) %% length(marks)) + 1; b <- ((i - 1) %/% length(marks)) + 1
    tie[names(marks)[k]] <- paste(names(marks)[a], names(marks)[b])
  }
}

# 6. Secondary weights. Marks: classes of secondary-equal neighbours in sorted
# order. Members of a group that do not decompose: secondary-equal classes;
# the class of the group's first member is the common secondary.
mclass <- cumsum(c(TRUE, !S(marks[-length(marks)], marks[-1])))
spell <- function(g) paste(reps[flat(g)], collapse = "")
# R (the atoms of group g) with `mark` after element i; i = 0 puts "x" + mark
# in front.
insert_mark <- function(g, i, mark) {
  at <- reps[flat(g)]
  if (i == 0) return(paste0("x", mark, paste(at, collapse = "")))
  paste0(paste(at[seq_len(i)], collapse = ""), mark, paste(at[-seq_len(i)], collapse = ""))
}
# Element that carries a non-common secondary: with the lowest mark m0 after
# element i, v sorts above the spelled string exactly when the secondary
# element lies at or before i (the secondary is below m0, which
# mark_position confirms).
own_element <- function(v, g) {
  k <- length(flat(g))
  if (k == 1) return(1L)
  for (i in seq_len(k)) if (cmp(v, insert_mark(g, i, marks[[1]])) > 0) return(i)
  NA_integer_
}
# Position of a non-common secondary relative to the marks: compare v (with
# "x" in front when it is on the first element) with the spelled string
# carrying each mark just before that element.
mark_position <- function(v, g, j) {
  x <- if (j == 1) paste0("x", v) else v
  r <- cmp(x, vapply(marks, function(m) insert_mark(g, j - 1L, m), ""))
  sum(r > 0)
}
members <- list()
for (g in seq_len(G)) {
  mem <- groups[[g]]
  plain <- mem[tie[names(mem)] == ""]
  cls <- cumsum(c(TRUE, if (length(plain) > 1) !S(plain[-length(plain)], plain[-1]) else logical()))
  # tertiary levels within each class
  for (k in unique(cls)) {
    x <- plain[cls == k]
    tl <- cumsum(c(TRUE, if (length(x) > 1) cmp(x[-length(x)], x[-1]) != 0 else logical()))
    common <- k == 1 && (atom[g] || S(x[[1]], spell(g)))
    # Otherwise: the atoms with one mark after element i (a skeleton) may be
    # secondary-equal to the class.
    skel <- NULL
    if (!common) {
      at <- reps[flat(g)]
      for (i in seq_along(at)) {
        sk <- paste0(paste(at[seq_len(i)], collapse = ""), marks, paste(at[-seq_len(i)], collapse = ""))
        hit <- which(S(sk, x[[1]]))
        if (length(hit)) { skel <- list(pos = i, mark = names(marks)[hit[1]]); break }
      }
    }
    own <- if (common || !is.null(skel)) NA_integer_ else own_element(x[[1]], g)
    pos <- if (is.na(own)) NA_integer_ else mark_position(x[[1]], g, own)
    for (j in seq_along(x)) members[[names(x)[j]]] <- list(group = g, class = k, common = common,
      skeleton = skel, level = tl[j], levels = max(tl), own = own, mark_position = pos,
      code = if (is.na(own)) NA_character_ else sprintf("B%d.%d", g, k))
  }
}
# Non-common secondaries on the same primary are compared with each other:
# for v (atoms R) and w (atoms Q) carrying the secondary on their first
# element, v + tail(Q) and w + tail(R) have the same primaries and differ
# first at that secondary. Secondary-equal pairs share one code.
own_classes <- Filter(function(m) !is.na(m$own), members)
own_classes <- own_classes[!duplicated(vapply(own_classes, `[[`, "", "code"))]
own_char <- vapply(names(own_classes), function(nm) chars[[nm]], "")
own_primary <- vapply(own_classes, function(m) primary[flat(m$group)[m$own]], 1L)
own_rank <- setNames(integer(length(own_classes)), vapply(own_classes, `[[`, "", "code"))
own_alias <- setNames(names(own_rank), names(own_rank))
for (q in unique(own_primary)) {
  idx <- which(own_primary == q)
  if (length(idx) < 2) next
  stopifnot(all(vapply(own_classes[idx], `[[`, 1L, "own") == 1L))
  tail_of <- function(i) paste(reps[flat(own_classes[[i]]$group)][-1], collapse = "")
  less <- function(i, j) cmp(paste0(own_char[[i]], tail_of(j)), paste0(own_char[[j]], tail_of(i)))
  same <- function(i, j) S(paste0(own_char[[i]], tail_of(j)), paste0(own_char[[j]], tail_of(i)))
  # insertion sort with the observed comparison
  ord <- integer()
  for (i in idx) {
    at <- length(ord) + 1L
    for (p in seq_along(ord)) if (less(i, ord[p]) < 0) { at <- p; break }
    ord <- append(ord, i, after = at - 1L)
  }
  for (p in seq_along(ord)) {
    own_rank[names(own_rank)[ord[p]]] <- p
    if (p > 1 && same(ord[p - 1], ord[p])) own_alias[names(own_rank)[ord[p]]] <- own_alias[names(own_rank)[ord[p - 1]]]
  }
}
# A shared secondary is named after a single-character class when there is one.
for (a in unique(own_alias)) {
  set <- names(own_alias)[own_alias == a]
  single_set <- set[vapply(own_classes[match(set, names(own_rank))], function(m) atom[m$group], TRUE)]
  if (length(single_set)) own_alias[set] <- single_set[[1]]
}
# Expansions with a common secondary: fit the tertiary of each element against
# the common classes of its atoms. The first element whose tertiary differs
# from every member decides; later elements are not observable against plain
# strings and take the same side (below all or above all members), offset by
# the character's level in its own class so that members of one class keep
# their observed order.
chain <- function(a) {
  g <- which(primary == a)
  x <- groups[[g]][vapply(names(groups[[g]]), function(nm) {
    m <- members[[nm]]; !is.null(m) && m$class == 1 && isTRUE(m$common)
  }, TRUE)]
  lv <- vapply(names(x), function(nm) members[[nm]]$level, 1L)
  x[!duplicated(lv)]
}
fit_tertiary <- function(c0, atoms, skel = NULL, bump = 0L, ch = lapply(atoms, chain)) {
  k <- length(atoms); t <- integer(k); pre <- character()
  # Spell element choices, with the skeleton's mark after its element.
  build <- function(parts) {
    if (!is.null(skel)) parts[skel$pos] <- paste0(parts[skel$pos], chars[[skel$mark]])
    paste(parts, collapse = "")
  }
  for (i in seq_len(k)) {
    rest <- seq_len(k)[-seq_len(i)]
    rmin <- vapply(ch[rest], function(x) x[[1]], "")
    rmax <- vapply(ch[rest], function(x) x[[length(x)]], "")
    x <- ch[[i]]
    lo <- cmp(c0, vapply(x, function(xx) build(c(pre, xx, rmin)), ""))
    hi <- cmp(c0, vapply(x, function(xx) build(c(pre, xx, rmax)), ""))
    j <- which(lo >= 0 & hi <= 0)
    if (length(j) == 1) { t[i] <- 2L * j; pre <- c(pre, x[[j]]); next }
    below <- sum(hi > 0)
    t[i] <- 2L * below + 1L
    for (r in rest) t[r] <- if (below == 0) 1L else 2L * length(ch[[r]]) + 1L + bump
    break
  }
  t
}

# 7. Assemble collation elements: primary, secondary, tertiary.
# Secondary codes: "C" common, "M<k>" mark class k, "B<g>.<k>.<pos>" a class of
# its own placed after <pos> mark classes. Tertiary: 2 * level within the class,
# odd values between levels.
ce <- list()
ce_of <- function(nm) {
  if (!is.null(ce[[nm]])) return(ce[[nm]])
  if (!is.na(tie[nm]) && tie[[nm]] != "") {
    parts <- strsplit(tie[[nm]], " ")[[1]]
    r <- do.call(c, lapply(parts, ce_of))
  } else if (nm %in% names(marks)) {
    k <- mclass[match(nm, names(marks))]
    x <- marks[mclass == k]; tl <- cumsum(c(TRUE, if (length(x) > 1) cmp(x[-length(x)], x[-1]) != 0 else logical()))
    r <- list(list(p = 0L, s = paste0("M", k), t = 2L * tl[match(nm, names(x))]))
  } else {
    m <- members[[nm]]; g <- m$group; p <- prim_of[[g]]
    if (m$common) {
      t <- if (length(p) == 1) rep(2L * m$level, 1) else fit_tertiary(chars[[nm]], p, bump = m$level)
      r <- lapply(seq_along(p), function(i) list(p = p[i], s = "C", t = t[i]))
    } else if (!is.null(m$skeleton)) {
      t <- fit_tertiary(chars[[nm]], p, m$skeleton, bump = m$level)
      r <- lapply(seq_along(p), function(i) list(p = p[i], s = "C", t = t[i]))
      r <- append(r, ce_of(m$skeleton$mark), after = m$skeleton$pos)
    } else {
      code <- own_alias[[m$code]]
      s1 <- paste0("B", code)
      if (code != m$code) {
        # Secondary-equal to a single character of another class: fit the
        # tertiary against that class at the carrying element.
        target <- own_classes[[match(code, names(own_rank))]]
        stopifnot(length(flat(target$group)) == 1L)
        tm <- names(Filter(function(x) identical(x$code, code), members))
        lv <- vapply(tm, function(z) members[[z]]$level, 1L)
        ch <- lapply(p, chain); ch[[m$own]] <- chars[tm][!duplicated(lv)]
        t <- fit_tertiary(chars[[nm]], p, bump = m$level, ch = ch)
      } else {
        t <- rep(2L * m$level, length(p))
      }
      r <- lapply(seq_along(p), function(i) list(p = p[i], s = if (i == m$own) s1 else "C", t = t[i]))
    }
  }
  ce[[nm]] <<- r
  r
}
for (nm in names(live)) ce_of(nm)

# 8. Contractions: pairs whose comparison changes when a completely ignorable
# character is put between them.
IGN <- chars[ignorable][[1]]
contr <- list()
for (x in live) {
  r <- cmp(paste0(x, live), paste0(x, IGN, live))
  for (y in live[r != 0]) {
    w <- paste0(x, y)
    eq <- live[cmp(live, w) == 0]
    contr[[length(contr) + 1]] <- data.frame(first = cps_hex(x), second = cps_hex(y),
      equals = if (length(eq)) names(eq)[1] else "")
  }
}
contr <- do.call(rbind, contr)
stopifnot(all(contr$equals != ""))

# 9. Write the observations and the derived table.
sec_codes <- unique(unlist(lapply(ce, function(r) vapply(r, `[[`, "", "s"))))
mark_codes <- sprintf("M%d", sort(unique(mclass)))
big <- setdiff(sec_codes, c("C", mark_codes))
big_class <- own_classes[match(sub("^B", "", big), names(own_rank))]
big_pos <- vapply(big_class, `[[`, 1L, "mark_position")
big_q <- own_primary[match(sub("^B", "", big), names(own_rank))]
big_r <- own_rank[sub("^B", "", big)]
big <- big[order(big_pos, big_q, big_r)]
# Order of secondaries: common, then marks with the non-common classes inserted
# after the number of mark classes they exceed. Non-common secondaries are only
# ever compared with others on the same primary, so their order across
# primaries is immaterial and follows the primary.
sec_order <- "C"
for (k in 0:length(mark_codes)) {
  if (k > 0) sec_order <- c(sec_order, mark_codes[k])
  sec_order <- c(sec_order, big[big_pos[order(big_pos, big_q, big_r)] == k])
}
sec_weight <- setNames(seq_along(sec_order), sec_order)
fmt <- function(r) paste(vapply(r, function(e) sprintf("%d.%d.%d", e$p, sec_weight[[e$s]], e$t), ""), collapse = " ")
table <- data.frame(cp = c(names(chars)[ignorable], names(live)),
  ces = c(rep("", sum(ignorable)), vapply(names(live), function(nm) fmt(ce[[nm]]), "")))
table <- table[order(strtoi(table$cp, 16L)), ]
header <- sprintf("# %s on %s; ICU %s; collation %s", R.version.string, R.version$platform,
  extSoftVersion()[["ICU"]], icuGetCollate("actual"))
w <- function(file, lines) writeLines(c(header, lines), file.path(out, file), useBytes = TRUE)
w("collation_table.tsv", c("# code point\tcollation elements primary.secondary.tertiary (empty: ignorable)",
  paste(table$cp, table$ces, sep = "\t")))
w("contractions.tsv", c("# first\tsecond\tequal to", paste(contr$first, contr$second, contr$equals, sep = "\t")))
rel <- c("=", "<", "<<", "<<<")
w("inventory_order.txt", c("# characters in R's sorted order; relation to the next: < primary, << secondary, <<< tertiary, = tie",
  paste(names(live), c(rel[lev + 1L], ""), sep = "\t"),
  "# completely ignorable", names(chars)[ignorable]))
w("decompositions.tsv", c("# group representative\tprimary decomposition (groups)\tatom",
  vapply(seq_len(G), function(g) paste(names(reps)[g], paste(names(reps)[decomp[[g]]], collapse = " "),
    atom[g], sep = "\t"), "")))
w("ties.tsv", c("# character\tequal to (character, or character + mark)",
  paste(names(tie)[tie != ""], tie[tie != ""], sep = "\t")))
# 10. Hangul syllables U+AC00..U+D7A3. The arithmetic of the Unicode
# Standard (section 3.12) spells each syllable as a leading consonant, a vowel
# and an optional trailing consonant from the conjoining jamo probed above;
# check that R ties every syllable with that jamo sequence.
s_index <- 0:11171
syll <- intToUtf8(0xAC00 + s_index, multiple = TRUE)
lead <- 0x1100 + s_index %/% (21 * 28); vowel <- 0x1161 + (s_index %% (21 * 28)) %/% 28
trail <- s_index %% 28
jamo <- vapply(seq_along(s_index), function(i) intToUtf8(c(lead[i], vowel[i], if (trail[i]) 0x11A7 + trail[i])), "")
hangul_ok <- cmp(syll, jamo) == 0
w("hangul.tsv", c("# syllables checked\tequal to their jamo sequence",
  paste(length(syll), sum(hangul_ok), sep = "\t"),
  if (!all(hangul_ok)) paste("# differs:", paste(sprintf("%04X", 0xAC00 + s_index[!hangul_ok]), collapse = " "))))
stopifnot(all(hangul_ok))
cat("inventory", length(chars), "ignorable", sum(ignorable), "marks", length(marks), "groups", G,
  "atoms", sum(atom), "contractions", nrow(contr), "secondaries", length(sec_order), "\n")
