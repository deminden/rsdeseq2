#!/usr/bin/env python3
"""Derive the rule data of rcompat-collate from black-box observations of R.

Inputs (default directory /tmp/rcollate, see docs/provenance/specs/collation.md):

- root_single.tsv: R's order of every code point under the ICU root collation,
  with the level of each difference (scripts/provenance/probe_collation_full.R);
- root_expansions.tsv: characters R sorts as several primary weights and the
  sequence (scripts/provenance/probe_collation_expansions.R);
- decompositions.tsv: canonical and compatibility decompositions and case
  from Rust (`cargo run -p rcompat-collate --example dump_decompositions`).

Questions the single-character order cannot answer (for example whether a
ligature carries a secondary weight of its own) are put to R with
scripts/provenance/probe_collation_oracle.R.

The model (docs/provenance/specs/collation.md) maps a character through rules:
canonical decomposition, a contraction list, explicit exceptions, the atom
order, marks, ignorables, compatibility decomposition with a tertiary variant,
lowercase mapping with the uppercase tertiary, and code point order for the
rest. This script stores only what the rules cannot produce, checks that the
rules reproduce R's order of every code point, and writes
crates/rcompat-collate/src/tables.rs.

Usage: build_collation_table.py [--input DIR] [--check]
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "crates/rcompat-collate/src/tables.rs"
ORACLE = ROOT / "scripts/provenance/probe_collation_oracle.R"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")


def hexs(seq) -> str:
    return " ".join(f"{c:04X}" for c in seq)


def parse_hex(text: str) -> list[int]:
    return [int(x, 16) for x in text.split()] if text else []


class Observations:
    """R's single-code-point order and the Rust decomposition facts."""

    def __init__(self, directory: Path):
        self.rows: list[tuple[int, str, bool]] = []
        self.ignorable: list[int] = []
        section = 0
        for line in (directory / "root_single.tsv").read_text(encoding="utf-8").splitlines():
            if line.startswith("# completely ignorable"):
                section = 1
                continue
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if section == 0:
                self.rows.append((int(f[0], 16), f[1], f[2] == "1"))
            else:
                self.ignorable.append(int(f[0], 16))
        # The trailing run in increasing code point order is the implicit
        # area (unassigned and private use); U+FFFD and U+FFFF follow it.
        n = len(self.rows)
        k = n - 3
        while k > 0 and self.rows[k - 1][0] < self.rows[k][0] and self.rows[k - 1][1] == "<":
            k -= 1
        self.head = k
        self.after_tail = [r[0] for r in self.rows[n - 2 :]]
        self.expansions: dict[int, list[int]] = {}
        for line in (directory / "root_expansions.tsv").read_text(encoding="utf-8").splitlines():
            if line.startswith("#"):
                continue
            a, b = line.split("\t")
            self.expansions[int(a, 16)] = parse_hex(b)
        self.nfd: dict[int, list[int]] = {}
        self.nfkd: dict[int, list[int]] = {}
        self.lower: dict[int, list[int]] = {}
        self.upper_flag: set[int] = set()
        for line in (directory / "decompositions.tsv").read_text(encoding="utf-8").splitlines():
            f = line.split("\t")
            c = int(f[0], 16)
            d, k2, lo, flags = parse_hex(f[1]), parse_hex(f[2]), parse_hex(f[3]), f[4]
            if d != [c]:
                self.nfd[c] = d
            if k2 != [c]:
                self.nfkd[c] = k2
            if lo != [c]:
                self.lower[c] = lo
            if "u" in flags:
                self.upper_flag.add(c)


CACHE: Path | None = None


def oracle(pairs, locale: str = "root") -> list[str]:
    """R's relation for each pair of strings ("<", "<<", "=", ">>>", ...).
    Answers are cached by the SHA-256 of the questions."""
    if not pairs:
        return []
    text = "".join(f"{hexs(x)}\t{hexs(y)}\n" for x, y in pairs)
    digest = hashlib.sha256((locale + "\n" + text).encode()).hexdigest()
    cached = CACHE / f"{digest}.tsv" if CACHE else None
    if cached is None or not cached.exists():
        with tempfile.TemporaryDirectory() as tmp:
            q = Path(tmp) / "q.tsv"
            a = Path(tmp) / "a.tsv"
            q.write_text(text, encoding="utf-8")
            subprocess.run([RSCRIPT, "--vanilla", str(ORACLE), "compare", locale, str(q), str(a)], check=True)
            answer = a.read_text(encoding="utf-8")
        if cached is not None:
            cached.parent.mkdir(parents=True, exist_ok=True)
            cached.write_text(answer, encoding="utf-8")
    else:
        answer = cached.read_text(encoding="utf-8")
    out = [line for line in answer.splitlines() if not line.startswith("#")]
    if len(out) != len(pairs):
        raise RuntimeError("oracle answered a different number of questions")
    return out


class Structure:
    """Groups (primary-equal), classes (secondary-equal) and layers
    (tertiary-equal) of R's single-code-point order."""

    def __init__(self, obs: Observations):
        rows = obs.rows[: obs.head]
        self.marks = [(c, rel) for c, rel, pign in rows if pign]
        letters = [(c, rel) for c, rel, pign in rows if not pign]
        if [r[0] for r in rows[: len(self.marks)]] != [m[0] for m in self.marks]:
            raise RuntimeError("primary-ignorable characters do not sort first")
        self.groups: list[list[tuple[int, str]]] = [[]]
        for c, rel in letters:
            self.groups[-1].append((c, rel))
            if rel == "<":
                self.groups.append([])
        if not self.groups[-1]:
            self.groups.pop()
        self.group_of: dict[int, int] = {}
        for g, members in enumerate(self.groups):
            for c, _ in members:
                self.group_of[c] = g
        # The representative of a group is its first member without a
        # canonical decomposition (R lists ties in code point order, so a
        # compatibility ideograph can come before the ideograph it equals).
        first = [members[0][0] for members in self.groups]
        self.rep = [next((c for c, _ in members if c not in obs.nfd), members[0][0]) for members in self.groups]
        self.atoms = [r for r, f in zip(self.rep, first) if f not in obs.expansions]
        self.atom_index = {c: i for i, c in enumerate(self.atoms)}
        to_rep = dict(zip(first, self.rep))
        self.atom_seq = [[r] if f not in obs.expansions else [to_rep[a] for a in obs.expansions[f]] for r, f in zip(self.rep, first)]

    @staticmethod
    def classes(members):
        """Split a group's members into classes, each a list of layers."""
        out = [[[]]]
        for c, rel in members:
            out[-1][-1].append(c)
            if rel == "<<":
                out.append([[]])
            elif rel == "<<<":
                out[-1].append([])
        return [[layer for layer in cls if layer] for cls in out if any(cls)]


def strongest(relations: list[str]) -> str:
    """The relation between the ends of a chain: its lowest level."""
    order = ["<", "<<", "<<<", "="]
    return min(relations, key=order.index)


class Model:
    """Python mirror of the rules in crates/rcompat-collate/src/elements.rs."""

    COMMON = 1

    def __init__(self, obs: Observations, st: Structure):
        self.obs, self.st = obs, st
        self.ignorable = set(obs.ignorable)
        # Marks: non-decomposing primary-ignorable characters in R's order.
        kept = []
        pending = []
        for c, rel in st.marks:
            if c in obs.nfd:
                pending.append(rel)
                continue
            if kept:
                kept[-1] = (kept[-1][0], strongest([kept[-1][1]] + pending))
            pending = [rel]
            kept.append((c, rel))
        self.mark_list = kept
        self.mark_set = {c for c, _ in kept}
        # Code points R sorts as unassigned although Rust's decomposition or
        # case data knows them (newer Unicode than the ICU of the probe).
        tail = {r[0] for r in obs.rows[obs.head :]}
        self.unassigned = {c for c in tail if c in obs.nfd or c in obs.nfkd or c in obs.lower or c in obs.upper_flag}
        if any(c in obs.nfd for c in self.unassigned):
            raise RuntimeError("an unassigned code point has a canonical decomposition")
        self.explicit: dict[int, list[tuple[int, int, int]]] = {}
        self.compat: dict[int, tuple[int, int]] = {}
        self.contracting: set[int] = set()
        self.extra_contractions: dict[tuple[int, ...], int] = {}
        self.max_own = 0
        self.prevowels: dict[int, list[tuple[int, int]]] = {}
        self.tl, self.tu = 0, 1
        self.mark_tertiary: dict[int, int] = {}

    # -- weights -------------------------------------------------------
    def primary(self, c: int) -> int:
        return (self.st.atom_index[c] + 1) << 8

    def tail_primary(self, c: int) -> int:
        if c == 0xFFFD:
            c = 0x110000
        elif c == 0xFFFF:
            c = 0x110001
        return ((len(self.st.atoms) + 1) << 8) + c

    def mark_weights(self):
        """Secondary code and tertiary layer of each mark."""
        base = self.COMMON + self.max_own + 1
        out = {}
        s, layer = base, 0
        for i, (c, _) in enumerate(self.mark_list):
            if i:
                rel = self.mark_list[i - 1][1]
                if rel in ("<", "<<"):
                    s, layer = s + 1, 0
                elif rel == "<<<":
                    layer += 1
            out[c] = (s, layer)
        return out

    def finish(self):
        self.marks = {c: (s, self.mark_tertiary[c]) for c, (s, _) in self.mark_weights().items()}
        # Primary-ignorable compatibility characters that decompose to
        # several marks sort as those marks.
        for c in self.mark_set:
            if c in self.obs.nfkd:
                parts = [x for x in self.nfd(self.obs.nfkd[c]) if x in self.mark_set]
                if len(parts) >= 2:
                    s0, t0 = self.marks[c]
                    self.explicit[c] = [(0, self.marks[x][0], t0 if i == 0 else self.marks[x][1]) for i, x in enumerate(parts)]
        self.contractions = {tuple(self.obs.nfd[c]): c for c in self.contracting}
        self.contractions.update(self.extra_contractions)
        self.max_key = max((len(k) for k in self.contractions), default=1)
        self.cache: dict[int, list[tuple[int, int, int]] | None] = {}

    # -- elements ------------------------------------------------------
    def char(self, c: int):
        if c in self.cache:
            return self.cache[c]
        self.cache[c] = None  # guards against cycles
        out = self._char(c)
        self.cache[c] = out
        return out

    def _char(self, c: int):
        if c in self.unassigned:
            return [(self.tail_primary(c), self.COMMON, self.tl)]
        if c in self.explicit:
            e = self.explicit[c]
            if isinstance(e, tuple):
                base = self.string(e[1])
                if base is None:
                    return None
                out, done = [], False
                for p, s, t in base:
                    if p and not done:
                        out.append((p, s, e[2]))
                        done = True
                    else:
                        out.append((p, s, t))
                return out
            return e
        if c in self.st.atom_index:
            return [(self.primary(c), self.COMMON, self.tl)]
        if c in self.marks:
            s, t = self.marks[c]
            return [(0, s, t)]
        if c in self.ignorable:
            return []
        if c in self.obs.nfkd:
            els = self.string(self.obs.nfkd[c])
            if els is None or c not in self.compat:
                return None
            lv, uv = self.compat[c]
            out = []
            for p, s, t in els:
                if p == 0:
                    out.append((p, s, t))
                elif t == self.tl:
                    out.append((p, s, lv))
                elif t == self.tu:
                    out.append((p, s, uv))
                else:
                    return None
            return out
        lower = self.obs.lower.get(c)
        if lower is not None and len(lower) == 1:
            els = self.char(lower[0])
            if els is None or any(p and t != self.tl for p, s, t in els):
                return None
            return [(p, s, self.tu if p else t) for p, s, t in els]
        return [(self.tail_primary(c), self.COMMON, self.tl)]

    def nfd(self, seq):
        out = []
        for c in seq:
            out.extend(self.obs.nfd.get(c, [c]))
        return out

    def reorder(self, seq):
        """A prevowel followed by one of its consonants sorts after it."""
        seq = list(seq)
        i = 0
        while i + 1 < len(seq):
            ranges = self.prevowels.get(seq[i])
            if ranges and any(a <= seq[i + 1] <= b for a, b in ranges):
                seq[i], seq[i + 1] = seq[i + 1], seq[i]
                i += 2
            else:
                i += 1
        return seq

    def string(self, seq):
        seq = self.reorder(self.nfd(seq))
        out = []
        i = 0
        while i < len(seq):
            for k in range(min(self.max_key, len(seq) - i), 1, -1):
                target = self.contractions.get(tuple(seq[i : i + k]))
                if target is not None:
                    els = self.char(target)
                    i += k
                    break
            else:
                els = self.char(seq[i])
                i += 1
            if els is None:
                return None
            out.extend(els)
        return out


def key(els):
    return (
        [p for p, s, t in els if p],
        [s for p, s, t in els if s],
        [t for p, s, t in els if t],
    )


def relation(a, b) -> str:
    """Relation of two element lists: "<", "<<", "<<<", "=" or ">..."."""
    ka, kb = key(a), key(b)
    for level in range(3):
        if ka[level] != kb[level]:
            sign = "<" if ka[level] < kb[level] else ">"
            return sign * (level + 1)
    return "="


def primaries_by_rule(obs: Observations, st: Structure, m: Model):
    """Atom sequence each head character gets from the rules alone, or None
    where no rule applies (the character must then be stored)."""
    contr = {tuple(obs.nfd[c]): c for c in m.contracting}
    max_key = max((len(k) for k in contr), default=1)
    memo: dict[int, list[int] | None] = {}

    def char(c):
        if c in memo:
            return memo[c]
        memo[c] = None
        if c in st.atom_index:
            r = [c]
        elif c in m.mark_set or c in m.ignorable:
            r = []
        elif c in obs.nfkd:
            r = seq(obs.nfkd[c])
        elif c in obs.lower and len(obs.lower[c]) == 1:
            r = char(obs.lower[c][0])
        else:
            r = None
        if r is None and c in st.group_of:
            r = st.atom_seq[st.group_of[c]] if c in stored else None
        memo[c] = r
        return r

    def seq(s):
        s = m.nfd(s)
        out = []
        i = 0
        while i < len(s):
            for k in range(min(max_key, len(s) - i), 1, -1):
                t = contr.get(tuple(s[i : i + k]))
                if t is not None:
                    r = char(t)
                    i += k
                    break
            else:
                r = char(s[i])
                i += 1
            if r is None:
                return None
            out.extend(r)
        return out

    stored: set[int] = set()
    return char, seq, stored, memo


PREVOWELS: dict[int, list[tuple[int, int]]] = {}


def find_prevowels(obs: Observations, st: Structure, width: int = 40):
    """Pairs of neighbouring atoms whose order changes when an ignorable
    character is put between them are contractions or reorderings. A
    reordering ties x y with y x; for each such x the following characters
    it swaps with are found among the atoms within 300 of it."""
    ign = obs.ignorable[0]
    # Ideographs (from U+4E00 on in the atom order) are left out.
    atoms = st.atoms[: st.atom_index[0x4E00]]
    pairs = []
    for i, x in enumerate(atoms):
        for d in range(-width, width + 1):
            if d and 0 <= i + d < len(atoms):
                pairs.append(([x, atoms[i + d]], [x, ign, atoms[i + d]]))
    hits = [q for (q, _), rel in zip(pairs, oracle(pairs)) if rel != "="]
    swaps = oracle([(q, q[::-1]) for q in hits])
    pre = sorted({q[0] for q, rel in zip(hits, swaps) if rel == "="})
    pos = {a: i for i, a in enumerate(atoms)}
    out = {}
    for x in pre:
        window = atoms[max(0, pos[x] - 300) : pos[x] + 300]
        rels = oracle([([x, y], [x, ign, y]) for y in window])
        cons = sorted(y for y, rel in zip(window, rels) if rel != "=")
        ranges: list[list[int]] = []
        for y in cons:
            if ranges and y == ranges[-1][1] + 1:
                ranges[-1][1] = y
            else:
                ranges.append([y, y])
        out[x] = [tuple(r) for r in ranges]
    contractions = len(hits) - sum(1 for rel in swaps if rel == "=")
    return out, len(pairs), contractions


def derive(obs: Observations, st: Structure, extra: set[int]) -> Model:
    m = Model(obs, st)
    marks = m.mark_weights()
    head = [c for c, _, pign in obs.rows[: obs.head] if not pign]
    # 1. Primary level. Canonical decompositions that R sorts as one letter
    # become contractions; characters no rule explains are stored.
    for _ in range(4):
        char, seq, stored, memo = primaries_by_rule(obs, st, m)
        changed = False
        for c in head:
            observed = st.atom_seq[st.group_of[c]]
            if c in obs.nfd:
                got = seq([c])
                if got != observed:
                    if c not in m.contracting:
                        m.contracting.add(c)
                        changed = True
            elif c not in st.atom_index:
                memo.clear()
                got = char(c)
                if got != observed:
                    stored.add(c)
        if not changed:
            break
    m.stored = stored | extra | {c for c in m.contracting if c not in st.atom_index}
    m.prevowels = PREVOWELS
    # Stored characters whose compatibility decomposition R ties with them
    # (as "l" + U+00B7 with U+0140) make that sequence a contraction.
    candidates = [c for c in sorted(m.stored) if c in obs.nfkd and c not in obs.nfd]
    answers = oracle([([c], obs.nfkd[c]) for c in candidates])
    for c, rel in zip(candidates, answers):
        if rel == "=":
            m.extra_contractions[tuple(m.nfd(obs.nfkd[c]))] = c
    # Characters whose compatibility or lowercase base carries a tertiary
    # weight other than plain or uppercase are stored too.
    while True:
        label = labels(obs, st, m)
        extra = {c for c in head if label(c)[0] == "?"}
        if not extra:
            break
        m.stored |= extra
    secondaries(obs, st, m)
    tertiaries(obs, st, m)
    return m


def first_char(obs: Observations, m: Model, c: int) -> int:
    """The character whose tertiary weight a canonical decomposition keeps."""
    if c in m.contracting or c not in obs.nfd:
        return c
    return obs.nfd[c][0]


def secondaries(obs: Observations, st: Structure, m: Model):
    """Element templates of the stored characters: (code point, secondary,
    tertiary label). A secondary of "own" is a weight of its own above the
    common one; "like:x" copies the elements of the rule-explained x."""
    m.templates: dict[int, list] = {}
    questions = []
    skeleton_questions = []
    m.skeleton_facts = []
    for g, members in enumerate(st.groups):
        atoms = st.atom_seq[g]
        own_rank = 0
        for k, cls in enumerate(Structure.classes(members)):
            flat = [c for layer in cls for c in layer]
            stored = [c for c in flat if c in m.stored]
            if not stored:
                continue
            ruled = [c for c in flat if c not in m.stored and c not in m.contracting]
            if k == 0 and len(atoms) == 1:
                kind = ("common", None)
            elif k == 0:
                # Expansion with no rule-explained member: ask R whether it is
                # secondary-equal to its spelled-out atoms. Rank 1 is kept
                # for it in case it is not.
                kind = ("ask", None)
                own_rank = 1
                for c in stored:
                    questions.append((c, atoms))
            elif ruled:
                kind = ("like", [ruled[0]])
            else:
                # A weight of its own, unless the class is secondary-equal to
                # the atoms with one mark (as "đ" is to "d" + U+0335).
                own_rank += 1
                kind = ("own", own_rank)
                skeleton_questions.append((stored[0], atoms, g, k))
            for c in stored:
                m.templates[c] = (kind, atoms)
    # Skeletons: the stored character against its atoms with each mark
    # inserted after the first atom.
    marks = [c for c, _ in m.mark_list]
    pairs = [([c], atoms[:1] + [mk] + atoms[1:]) for c, atoms, _, _ in skeleton_questions for mk in marks]
    answers = oracle(pairs)
    for i, (c, atoms, g, k) in enumerate(skeleton_questions):
        rels = answers[i * len(marks) : (i + 1) * len(marks)]
        hit = [mk for mk, rel in zip(marks, rels) if rel == "="]
        hit += [mk for mk, rel in zip(marks, rels) if rel in ("<<<", ">>>")]
        if hit:
            # The skeleton has a plain first letter: the class member R ties
            # with it carries the plain tertiary weight, one above it more.
            rel = rels[marks.index(hit[0])]
            m.skeleton_facts.append((c, rel))
        if hit:
            cls = Structure.classes(st.groups[g])[k]
            for layer in cls:
                for x in layer:
                    if x in m.stored:
                        m.templates[x] = (("like", atoms[:1] + [hit[0]] + atoms[1:]), atoms)
    answers = oracle([([c], atoms) for c, atoms in questions])
    for (c, atoms), rel in zip(questions, answers):
        if rel in ("<<<", ">>>", "="):
            m.templates[c] = (("common", None), atoms)
        elif rel in ("<<", ">>"):
            m.templates[c] = (("own", 1), atoms)
        else:
            raise RuntimeError(f"{c:04X} is not primary-equal to its atoms")
    m.max_own = max([k[1] for k, _ in m.templates.values() if k[0] == "own"], default=0)


def labels(obs: Observations, st: Structure, m: Model):
    """Tertiary label of every head character: "L" plain, "U" uppercase,
    "K<cp>" a compatibility variant (with the case of its base), "V<cp>" a
    stored character, "M<cp>" a mark variant."""
    contr = {tuple(obs.nfd[c]): c for c in m.contracting}
    memo: dict[int, tuple[str, str]] = {}

    def label(c):
        # returns (label, slot) where slot is the case of the base
        if c in memo:
            return memo[c]
        memo[c] = ("L", "L")
        if c in m.stored:
            r = (f"V{c:X}", "V")
        elif c in st.atom_index:
            r = ("L", "L")
        elif c in obs.nfd and c not in m.contracting:
            s = obs.nfd[c]
            for k in range(len(s), 1, -1):
                if tuple(s[:k]) in contr:
                    r = label(contr[tuple(s[:k])])
                    break
            else:
                r = label(s[0])
        elif c in obs.nfkd:
            base = [x for x in m.nfd(obs.nfkd[c]) if x not in m.mark_set]
            slot = label(base[0])[1] if base else "L"
            r = (f"K{c:X}", slot) if slot in ("L", "U") else ("?", "?")
        elif c in obs.lower and len(obs.lower[c]) == 1:
            r = ("U", "U") if label(obs.lower[c][0])[1] == "L" else ("?", "?")
        else:
            r = ("L", "L")
        memo[c] = r
        return r

    return label


def fit_expansions(obs: Observations, st: Structure, m: Model, label):
    """Place the tertiary weight of characters that sort as several atoms
    (ligatures, squared and parenthesised forms, stored expansions) against
    the variants of their first letter. With x_j the first letter replaced
    by the class member of layer j and the rest unchanged, R puts the
    character after x_j and before x_(j+1) when its weight is that of layer
    j or lies between j and j+1; comparing with x_j followed by the highest
    variant of the second letter tells the two apart."""
    items = []
    for c in sorted(set(obs.nfkd) | m.stored):
        if c in obs.nfd or c not in st.group_of:
            continue
        lab = label(c)[0]
        if not (lab.startswith("K") or lab.startswith("V")):
            continue
        spelled = m.nfd(obs.nfkd[c]) if c in obs.nfkd and c not in m.stored else list(st.atom_seq[st.group_of[c]])
        bases = [i for i, x in enumerate(spelled) if x not in m.mark_set and x not in m.ignorable]
        if len(bases) < 2:
            continue
        first, second = spelled[bases[0]], spelled[bases[1]]
        if first not in st.group_of or second not in st.group_of:
            continue
        cls = Structure.classes(st.groups[st.group_of[first]])[0]
        reps = [layer[0] for layer in cls]
        top2 = Structure.classes(st.groups[st.group_of[second]])[0][-1][0]
        items.append((c, lab, spelled, bases, reps, top2))
    pairs = []
    for c, lab, spelled, bases, reps, top2 in items:
        for r in reps:
            x = list(spelled)
            x[bases[0]] = r
            pairs.append(([c], x))
            y = list(x)
            y[bases[1]] = top2
            pairs.append(([c], y))
    answers = iter(oracle(pairs))
    out = []
    fitted = set()
    for c, lab, spelled, bases, reps, top2 in items:
        fitted.add(c)
        rel = [(next(answers), next(answers)) for _ in reps]
        if any(a.strip("<>=") or a in ("<", ">", "<<", ">>") for a, _ in rel):
            continue  # not secondary-equal: nothing to fit at the tertiary level
        above = [i for i, (a, _) in enumerate(rel) if a in (">>>", "=")]
        if not above:
            out.append([[lab], [label(reps[0])[0]]])
            continue
        j = above[-1]
        lo = label(reps[j])[0]
        if rel[j][1] in ("<<<", "<", "<<"):
            out.append(("=", lo, lab))  # the weight of layer j
        elif j + 1 < len(reps):
            out.append([[lo], [lab], [label(reps[j + 1])[0]]])
        else:
            out.append([[lo], [lab]])
    return out, fitted


def tertiaries(obs: Observations, st: Structure, m: Model):
    label = labels(obs, st, m)
    parent: dict[str, str] = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    # Characters that sort as several atoms are placed by fitting their first
    # element; their order among themselves follows from the later elements.
    facts, fitted = fit_expansions(obs, st, m, label)
    sequences = []
    for members in st.groups:
        for cls in Structure.classes(members):
            seq = [[label(c)[0] for c in layer if c not in fitted or c in m.stored] for layer in cls]
            sequences.append([layer for layer in seq if layer])
    mark_layers = [[[]]]
    for i, (c, rel) in enumerate(m.mark_list):
        mark_layers[-1][-1].append(c)
        if rel == "<<" or rel == "<":
            mark_layers.append([[]])
        elif rel == "<<<":
            mark_layers[-1].append([])
    m.mark_label = {}
    for cls in mark_layers:
        cls = [layer for layer in cls if layer]
        for li, layer in enumerate(cls):
            for c in layer:
                m.mark_label[c] = "L" if li == 0 else f"M{c:X}"
        if cls:
            sequences.append([[m.mark_label[c] for c in layer] for layer in cls])
    for seq in sequences:
        for layer in seq:
            for x in layer[1:]:
                union(layer[0], x)
    for c, rel in m.skeleton_facts:
        if rel == "=":
            union("L", label(c)[0])
        elif rel == ">>>":
            sequences.append([["L"], [label(c)[0]]])
        else:
            sequences.append([[label(c)[0]], ["L"]])
    for fact in facts:
        if os.environ.get("DEBUG_FIT") and any(x in str(fact) for x in os.environ["DEBUG_FIT"].split(",")):
            print("fit", fact, file=sys.stderr)
        if isinstance(fact, tuple):
            union(fact[1], fact[2])
        else:
            sequences.append(fact)
    edges: dict[str, set[str]] = defaultdict(set)
    nodes = set()
    for seq in sequences:
        for a, b in zip(seq, seq[1:]):
            edges[find(a[0])].add(find(b[0]))
        nodes.update(find(layer[0]) for layer in seq)
    # Longest-path layering (Kahn's algorithm); a cycle is an error.
    indeg = defaultdict(int)
    for a in edges:
        for b in edges[a]:
            indeg[b] += 1
    value = {n: 1 for n in nodes}
    ready = [n for n in nodes if indeg[n] == 0]
    done = 0
    while ready:
        n = ready.pop()
        done += 1
        for b in edges.get(n, ()):
            value[b] = max(value[b], value[n] + 1)
            indeg[b] -= 1
            if indeg[b] == 0:
                ready.append(b)
    if done != len(nodes):
        cyc = [n for n in nodes if indeg[n] > 0]
        raise RuntimeError(f"tertiary labels form a cycle: {sorted(cyc)[:20]}")
    val = lambda lab: value.get(find(lab), 1)
    m.tl, m.tu = val("L"), val("U")
    m.label, m.val = label, val
    m.tertiary_levels = max(value.values())


def assemble(obs: Observations, st: Structure, m: Model):
    """Turn templates and labels into weights and check every adjacent pair
    of R's single-code-point order."""
    m.mark_tertiary = {c: m.val(m.mark_label[c]) for c, _ in m.mark_list}
    for c, (kind, atoms) in m.templates.items():
        t = m.val(m.label(c)[0])
        if kind[0] == "like":
            m.explicit[c] = ("like", kind[1], t)
        else:
            s0 = m.COMMON + (kind[1] if kind[0] == "own" else 0)
            m.explicit[c] = [(m.primary(a), s0 if i == 0 else m.COMMON, t) for i, a in enumerate(atoms)]
    known: dict[int, dict[str, int]] = {}
    for c in obs.nfkd:
        if c in m.stored or c in obs.nfd:
            continue
        lab, slot = m.label(c)
        if lab.startswith("K"):
            known[c] = {slot: m.val(lab)}
    # The case slot a character does not use is borrowed from the nearest
    # compatibility character that uses it, so that runs stay long.
    keys = sorted(known)
    for i, c in enumerate(keys):
        for slot in ("L", "U"):
            if slot in known[c]:
                continue
            for d in range(1, 64):
                for j in (i - d, i + d):
                    if 0 <= j < len(keys) and slot in known[keys[j]]:
                        known[c][slot] = known[keys[j]][slot]
                        break
                if slot in known[c]:
                    break
            known[c].setdefault(slot, m.tu if slot == "U" else m.tl)
    m.compat = {c: (v["L"], v["U"]) for c, v in known.items()}
    m.finish()


def longest_increasing(seq):
    """One longest strictly increasing subsequence."""
    import bisect

    tails, tails_idx, prev = [], [], [-1] * len(seq)
    for i, x in enumerate(seq):
        k = bisect.bisect_left(tails, x)
        if k == len(tails):
            tails.append(x)
            tails_idx.append(i)
        else:
            tails[k] = x
            tails_idx[k] = i
        prev[i] = tails_idx[k - 1] if k else -1
    out, i = [], tails_idx[-1] if tails_idx else -1
    while i >= 0:
        out.append(seq[i])
        i = prev[i]
    return out[::-1]


def runs(values):
    """Maximal runs of consecutive code points: [(first, length)]."""
    out = []
    for c in values:
        if out and c == out[-1][0] + out[-1][1]:
            out[-1][1] += 1
        else:
            out.append([c, 1])
    return out


def rust_array(name: str, typ: str, items: list[str], doc: str, per_line: int = 8) -> str:
    lines = [f"/// {line}" for line in doc.splitlines()]
    lines.append(f"pub(crate) static {name}: [{typ}; {len(items)}] = [")
    for i in range(0, len(items), per_line):
        lines.append("    " + ", ".join(items[i : i + per_line]) + ",")
    lines.append("];")
    return "\n".join(lines) + "\n"


def emit(obs: Observations, st: Structure, m: Model) -> tuple[str, dict]:
    """The Rust tables and their sizes in bytes."""
    sizes = {}
    parts = []
    han_at = st.atom_index[0x4E00]
    plain = runs(st.atoms[:han_at])
    atom_runs = []
    for first, n in plain:
        while n > 0:
            k = min(n, 2048)
            atom_runs.append(f"0x{first << 11 | (k - 1):08X}")
            first += k
            n -= k
    parts.append(rust_array("ATOM_RUNS", "u32", atom_runs,
        "Atoms (characters with a primary weight of their own) before the ideographs, in R's\nprimary order, as runs of consecutive code points: `first << 11 | (length - 1)`."))
    sizes["ATOM_RUNS"] = 4 * len(atom_runs)
    # Ideographs: merged from blocks of consecutive code points.
    han = st.atoms[han_at:]
    blocks = runs(sorted(han))
    block_of = {}
    for b, (first, n) in enumerate(blocks):
        for c in range(first, first + n):
            block_of[c] = b
    if len(blocks) > 31:
        raise RuntimeError("too many ideograph blocks for the encoding")
    # Within each block R's order is mostly code point order: the characters
    # outside the longest increasing subsequence are listed where they occur
    # ("moved"); the rest are runs that continue where their block stopped.
    moved = set()
    for b, (first, n) in enumerate(blocks):
        seq = [c for c in han if block_of[c] == b]
        moved |= set(seq) - set(longest_increasing(seq))
    han_runs, moved_list = [], []
    cursor = [first for first, _ in blocks]
    prev_block, length = None, 0
    def flush():
        if prev_block is not None and length:
            han_runs.append(f"0x{prev_block << 10 | (length - 1):04X}")
    for c in han:
        if c in moved:
            flush()
            prev_block, length = None, 0
            han_runs.append(f"0x{0x8000 | block_of[c]:04X}")
            moved_list.append(f"{c - blocks[block_of[c]][0]}")
            continue
        b = block_of[c]
        while cursor[b] in moved:
            cursor[b] += 1
        if cursor[b] != c:
            raise RuntimeError("ideograph order is not a merge of increasing block sequences")
        cursor[b] += 1
        if b == prev_block and length < 1024:
            length += 1
        else:
            flush()
            prev_block, length = b, 1
    flush()
    parts.append(rust_array("HAN_BLOCKS", "(u32, u32)", [f"(0x{f:05X}, {n})" for f, n in blocks],
        "Blocks of ideographs that R orders by radical and stroke: (first code point, length).", 4))
    parts.append(rust_array("HAN_RUNS", "u16", han_runs,
        "Ideographs in R's primary order after the atoms of `ATOM_RUNS`: an entry takes the next\n`(entry & 0x3FF) + 1` code points of block `entry >> 10` that are not moved; an entry\n`0x8000 | block` takes the next offset of `HAN_MOVED` in that block.", 12))
    parts.append(rust_array("HAN_MOVED", "u16", moved_list,
        "Ideographs R does not order by code point within their block, in R's order, as offsets\nfrom the start of their block.", 16))
    sizes["HAN"] = 8 * len(blocks) + 2 * len(han_runs) + 2 * len(moved_list)
    # Marks.
    marks = []
    for i, (c, _) in enumerate(m.mark_list):
        new = i > 0 and m.mark_list[i - 1][1] in ("<", "<<")
        t = m.marks[c][1]
        marks.append(f"0x{c << 8 | int(new) << 7 | t:08X}")
    parts.append(rust_array("MARKS", "u32", marks,
        "Characters without a primary weight in R's order: `code point << 8 | new << 7 | tertiary`,\nwhere `new` starts the next secondary weight."))
    sizes["MARKS"] = 4 * len(marks)
    ign = runs(sorted(obs.ignorable))
    parts.append(rust_array("IGNORABLE", "(u32, u32)", [f"(0x{f:04X}, 0x{f + n - 1:04X})" for f, n in ign],
        "Completely ignorable code points: (first, last).", 4))
    sizes["IGNORABLE"] = 8 * len(ign)
    un = runs(sorted(m.unassigned))
    parts.append(rust_array("UNASSIGNED", "(u32, u32)", [f"(0x{f:04X}, 0x{f + n - 1:04X})" for f, n in un],
        "Code points R sorts as unassigned although Rust's Unicode data knows them: (first, last).", 4))
    sizes["UNASSIGNED"] = 8 * len(un)
    parts.append(rust_array("CONTRACTING", "u32", [f"0x{c:04X}" for c in sorted(m.contracting)],
        "Characters whose canonical decomposition R sorts as the single character."))
    extra = sorted(m.extra_contractions.items())
    if any(len(k) != 2 for k, _ in extra):
        raise RuntimeError("only two-character extra contractions are encoded")
    parts.append(rust_array("PAIR_CONTRACTIONS", "(u32, u32, u32)", [f"(0x{k[0]:04X}, 0x{k[1]:04X}, 0x{v:04X})" for k, v in extra],
        "Two-character sequences R sorts as one character: (first, second, character).", 4))
    pre = [f"(0x{x:04X}, 0x{a:04X}, 0x{b:04X})" for x, ranges in sorted(m.prevowels.items()) for a, b in ranges]
    parts.append(rust_array("PREVOWELS", "(u32, u32, u32)", pre,
        "Prevowels and the consonants they sort after: (prevowel, first consonant, last consonant).", 4))
    sizes["CONTRACTIONS"] = 4 * len(m.contracting) + 12 * len(extra) + 12 * len(pre)
    # Stored characters.
    atom_of_primary = {m.primary(a): a for a in st.atoms}
    mark_of_secondary = {}
    for c, _ in m.mark_list:
        mark_of_secondary.setdefault(m.marks[c][0], c)

    def encode(e):
        p, s, t = e
        if p:
            ref, sec = atom_of_primary[p], s - m.COMMON
        else:
            ref, sec = mark_of_secondary[s], 0
        if not (0 <= sec < 32 and 0 < t < 64):
            raise RuntimeError(f"weight out of range: {e}")
        return ref << 11 | sec << 6 | t

    resolved = {}
    for c in sorted(m.explicit):
        els = m.char(c)
        if els is None:
            raise RuntimeError(f"{c:04X} has no elements")
        resolved[c] = [encode(e) for e in els]
    singles = sorted(c for c, v in resolved.items() if len(v) == 1 and v[0] >> 11 in m.st.atom_index)
    variant_runs = []
    for c in singles:
        v = resolved[c][0]
        if variant_runs:
            first, n, v0 = variant_runs[-1]
            if c == first + n and v == v0 + (n << 11):
                variant_runs[-1][1] += 1
                continue
        variant_runs.append([c, 1, v])
    in_runs = set(singles)
    parts.append(rust_array("VARIANTS", "(u32, u32)", [f"(0x{f << 11 | (n - 1):08X}, 0x{v:08X})" for f, n, v in variant_runs],
        "Stored characters with one element, as runs: (`first << 11 | (length - 1)`, element of\nthe first); the atom of each further character is the next code point.", 4))
    multi = [c for c in sorted(resolved) if c not in in_runs]
    offsets, elements = [], []
    for c in multi:
        offsets.append(f"(0x{c:05X}, {len(elements)})")
        elements.extend(f"0x{v:08X}" for v in resolved[c])
    offsets.append(f"(0x{0x110000:05X}, {len(elements)})")
    parts.append(rust_array("EXPLICIT", "(u32, u16)", offsets,
        "Other stored characters: (code point, index of the first element in `ELEMENTS`); a\nfinal entry closes the list.", 6))
    parts.append(rust_array("ELEMENTS", "u32", elements,
        "Collation elements of stored characters: `reference << 11 | secondary << 6 | tertiary`.\nThe reference is an atom (primary weight; secondary 0 is the common weight, others are\nweights of its own) or a mark (secondary weight of the mark)."))
    sizes["STORED"] = 8 * len(variant_runs) + 8 * len(offsets) + 4 * len(elements)
    # Compatibility tertiary weights as runs over code points.
    comp = []
    for c in sorted(m.compat):
        v = m.compat[c]
        if not comp or comp[-1][1] != v:
            comp.append((c, v))
        if max(v) >= 32:
            raise RuntimeError("tertiary weight out of range")
    parts.append(rust_array("COMPAT", "u32", [f"0x{c << 10 | v[0] << 5 | v[1]:08X}" for c, v in comp],
        "Tertiary weights of compatibility decompositions: `first << 10 | lowercase << 5 |\nuppercase`, for the characters from `first` up to the next entry."))
    sizes["COMPAT"] = 4 * len(comp)
    head = (
        "// @generated by scripts/provenance/build_collation_table.py from black-box observations\n"
        "// of R 4.6.2 with ICU 78.3; see docs/provenance/specs/collation.md. Do not edit.\n\n"
        "//! Data the rules of the root collation cannot produce.\n\n"
        f"/// Number of atoms (`ATOM_RUNS` and the ideographs).\npub(crate) const ATOM_COUNT: u32 = {len(st.atoms)};\n"
        f"/// Tertiary weight of plain (lowercase or uncased) letters.\npub(crate) const PLAIN: u16 = {m.tl};\n"
        f"/// Tertiary weight of uppercase letters.\npub(crate) const UPPER: u16 = {m.tu};\n"
        f"/// Number of secondary weights of their own between the common one and the marks.\npub(crate) const OWN_SECONDARIES: u16 = {m.max_own};\n\n"
    )
    return head + "\n".join(parts), sizes


def check_single(obs: Observations, m: Model):
    """Adjacent pairs of R's order where the model differs."""
    rows = [(c, rel) for c, rel, _ in obs.rows]
    bad = []
    prev = None
    for c, rel in rows:
        els = m.string([c])
        if prev is not None:
            pc, prel, pels = prev
            want = "=" if prel == "=" else prel
            got = "?" if els is None or pels is None else relation(pels, els)
            if got != want:
                bad.append((pc, c, want, got))
        prev = (c, rel, els)
    return bad


def check_sorted(m: Model, path: Path):
    """Adjacent pairs of an R-sorted list of strings where the model differs."""
    bad, total = [], 0
    prev = None
    text = path.read_text(encoding="utf-8")
    levels = "\t<<" in text
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        h, rel = line.split("\t")
        els = m.string(parse_hex(h))
        if prev is not None:
            total += 1
            got = "?" if els is None or prev[1] is None else relation(prev[1], els)
            want = prev[2]
            if not levels and got.startswith("<"):
                got = "<"
            if got != want:
                bad.append((prev[0], h, want, got))
        prev = (h, els, rel)
    return total, bad


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=Path("/tmp/rcollate"))
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--verify", type=Path, action="append", default=[])
    args = parser.parse_args()
    global CACHE
    CACHE = args.input / "oracle_cache"
    obs = Observations(args.input)
    st = Structure(obs)
    PREVOWELS.update(find_prevowels(obs, st)[0])
    print(f"head {obs.head} groups {len(st.groups)} atoms {len(st.atoms)} marks {len(st.marks)}", file=sys.stderr)
    # Characters a rule explains at the primary level but not at the others
    # are stored as well, until R's order of every code point is reproduced.
    extra: set[int] = set()
    for _ in range(8):
        m = derive(obs, st, extra)
        assemble(obs, st, m)
        bad = check_single(obs, m)
        new = {
            c
            for a, b, _, _ in bad
            for c in (a, b)
            if c not in m.stored and c not in obs.nfd and c not in st.atom_index and c not in m.mark_set
        }
        print(f"stored {len(m.stored)} mismatches {len(bad)} newly stored {len(new)}", file=sys.stderr)
        if not new:
            break
        extra |= new
    print(f"contracting {len(m.contracting)} stored {len(m.stored)} own {m.max_own} tertiary levels {m.tertiary_levels} tl {m.tl} tu {m.tu}", file=sys.stderr)
    print(f"single-code-point mismatches {len(bad)}", file=sys.stderr)
    for a, b, want, got in bad[:60]:
        print(f"  {a:04X} {b:04X} R {want} model {got}", file=sys.stderr)
    text, sizes = emit(obs, st, m)
    with tempfile.TemporaryDirectory() as tmp:
        formatted = Path(tmp) / "tables.rs"
        formatted.write_text(text, encoding="utf-8")
        subprocess.run(["rustfmt", "--edition", "2024", str(formatted)], check=True)
        text = formatted.read_text(encoding="utf-8")
    print("data bytes " + ", ".join(f"{k} {v}" for k, v in sizes.items()) + f"; total {sum(sizes.values())}", file=sys.stderr)
    if args.check:
        if OUTPUT.read_text(encoding="utf-8") != text:
            print(f"{OUTPUT} is not current", file=sys.stderr)
            return 1
    elif not bad:
        OUTPUT.write_text(text, encoding="utf-8")
    for path in args.verify:
        total, bad2 = check_sorted(m, path)
        print(f"{path.name}: {total} adjacent pairs, {len(bad2)} differ", file=sys.stderr)
        for a, b, want, got in bad2[:40]:
            print(f"  [{a}] [{b}] R {want} model {got}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
