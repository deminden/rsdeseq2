# Specification: R's string collation (C and ICU root)

Status: written from allowed inputs only: R's help pages (`?Comparison`,
`?sort`, `?icuSetCollate`), the general description of multi-level collation
(base letters, then accents, then case and variants), the Unicode Standard's
arithmetic for spelling Hangul syllables with conjoining jamo (section 3.12),
and black-box runs of R 4.6.2 with ICU 78.3 on macOS arm64. No collation
implementation, table or locale data was consulted: no ICU or ICU4X source or
data, no CLDR files, no DUCET (`allkeys.txt`), no glibc locale files and no R
source. Every weight in `crates/rcompat-collate/src/table.rs` is derived from
R's answers by the probe scripts below.

## Purpose

R orders character data with the collation of the session. Formula
evaluation needs that order for factor levels (`factor()` sorts the unique
values) and for string comparisons (`<`, `>`, `<=`, `>=`). R reports the
collation in use with `icuGetCollate("actual")`; the rsdeseq2 R wrapper passes
`"C"` when ICU is not in use and the reported ICU locale otherwise. For
English and most western European locales ICU reports `"root"`.

## Probes

- `scripts/provenance/probe_collation.R` derives the weight table.
- `scripts/provenance/probe_collation_orders.R` records R's sorted order of
  about 32,000 chosen and random strings, and the collation ICU reports for 77
  locale identifiers. These orders check the model; they are not used to
  derive it (the near-duplicate set uses the table only to choose inputs).
- `scripts/provenance/build_collation_table.py` turns the derived table into
  `crates/rcompat-collate/src/table.rs` (`--check` verifies it is current).

Outputs: `docs/provenance/observations/collation_r462_macos/`. The sorted sets
are copied to `crates/rcompat-collate/tests/data/` for the crate's tests.

## Observed facts

Probed code points: U+0001..U+036F (ASCII, Latin-1, Latin Extended-A/B, IPA,
spacing modifiers, combining diacritics), U+0370..U+052F (Greek, Cyrillic),
U+0590..U+06FF (Hebrew, Arabic), U+0900..U+097F (Devanagari), U+1100..U+11FF
(conjoining jamo), U+1E00..U+1EFF (Latin Extended Additional), U+2000..U+22FF
(punctuation, super- and subscripts, currency, letterlike symbols, number
forms, arrows, mathematical operators), U+2300..U+23FF, U+2460..U+27BF
(enclosed alphanumerics, box drawing, shapes, symbols, dingbats),
U+3040..U+30FF (kana), U+FB00..U+FB06 (Latin ligatures), U+FF01..U+FF5E
(fullwidth ASCII), U+1F300..U+1F64F (pictographs, emoticons): 5,364 code
points.

1. Ties are reported as ties. `x < y` and `y < x` are both `FALSE` for, for
   example, `"a\u0301"` and `"\u00e1"`, or `"a\u00adb"` and `"ab"`; `order()`
   with a second key then decides by that key.
   `icuSetCollate(strength = ...)` has no observable effect in this build,
   so levels are separated with the string tests below rather than with
   strength settings.
2. 162 probed code points are completely ignorable: appending them to `"a"`
   gives a string that ties with `"a"` (C0 and C1 controls other than
   tab..carriage return and U+0085, soft hyphen, zero-width and directional
   format characters, U+034F, and others).
3. Tab (U+0009) is the lowest non-ignorable character, then the other
   whitespace controls, then spaces, punctuation, symbols, currency signs,
   digits and letters. Punctuation and spaces are not ignored: `"a c"` <
   `"a-d"` < `"ab"`.
4. Level tests. With `LO` = tab and `HI` = line feed, two strings x and y
   have equal base (primary) weights exactly when `x+LO+HI > y+LO+LO` and
   `x+LO+LO < y+LO+HI` (tab is never part of the strings tested). With equal
   primaries, they also have equal accent (secondary) weights exactly when
   `x+"a\u0301" > y+"a"` and `y+"a\u0301" > x+"a"`. Checked on `a`/`A`
   (tertiary), `a`/`á` (secondary) and `a`/`b` (primary).
5. 184 probed characters carry no primary weight (`"a"+c` is primary-equal
   to `"a"`): the combining marks and a few others. Their secondary order is
   their sorted order; all of them sort above the plain secondary of letters.
6. The remaining characters fall into 3,540 primary groups. 3,368 groups are
   single primary weights; 172 are expansions, primary-equal to a sequence
   of other groups found by prefix search (`g <=p c <p g + U+FFFF`): for
   example `ß` = s s, `æ` = a e, `œ` = o e, `½` = 1 ⁄ 2, `…` = . . .,
   `ﬃ` = f f i, `ǳ` = d z, `№` = n o, `™` = t m (`decompositions.tsv`).
7. 816 characters tie with another character, or with another character of
   their group followed by one combining mark: `á` = `a`+U+0301,
   `ǖ` = `ü`+U+0304, `ł` = `l`+U+0335, `ø` = `o`+U+0338, `ǆ` = `ǳ`+U+030C;
   the mark U+0344 ties with U+0308 U+0301 (`ties.tsv`).
8. Within a group, the characters that do not decompose form classes of
   equal secondary weight (neighbours in sorted order that pass the
   secondary test). The first class holds the group's first character and
   has the plain secondary weight; tertiary levels within a class are the
   sorted order of its members (`a` < `ａ` < U+0363 < `A` < `Ａ` < `ª` < `ₐ`).
9. A further class is either secondary-equal to the group's primaries with
   one combining mark inserted (a skeleton, for example `ǆ` against `dz` +
   U+030C), or has a secondary weight of its own. The weights of its own
   (`ð`, `ſ`, `ŀ`, `ß`, `æ`, `œ`, `ґ`, the curly quotes, ...) all sort
   between the plain secondary and every combining mark: with the lowest
   mark m0, `"x"+v` < `"x"+m+R` for every mark m, where R spells the
   primaries of v. The element that carries it is found with m0 inserted
   after each element in turn. Weights of their own on the same primary are
   ordered by comparing `v + tail(w)` with `w + tail(v)` (for example
   `ß` < `ſ`), and a class secondary-equal to another (`ﬅ` and `ſt`) shares
   its weight.
10. Contractions: inserting a completely ignorable character between two
    characters changes the comparison for 11 pairs, each equal to one
    character: `l·`, `L·` (U+00B7 or U+0387) = `ŀ`, `Ŀ`; `и`+U+0306 = `й`,
    `И`+U+0306 = `Й`; Arabic alef, waw and yeh with U+0653..U+0655 = the
    precomposed letters. No other probed pair contracts.
11. Every Hangul syllable U+AC00..U+D7A3 ties with its leading consonant,
    vowel and optional trailing consonant from U+1100..U+11FF
    (`hangul.tsv`: 11,172 of 11,172).
12. Locales: ICU reports `root` for `""`, `root`, `en` and every probed
    English region, `de` (DE, AT, CH), `fr` (FR, BE, CH), `it`, `nl`, `pt`
    (BR, PT), `ca`, `eu`, `ga`, `id`, `ms` and `sw`. It reports a tailoring of
    its own for `fr_CA`, `en_US_POSIX`, `es`, `gl`, `cy`, `af`, the Nordic,
    Baltic, Slavic and other languages probed (`locales.tsv`).

## Model

Each code point maps to a list of collation elements (primary, secondary,
tertiary); a weight of zero means "no weight at this level".

- Primary: 1..3,368, one per single-primary group in sorted order; 0 for
  marks; expansions list the primaries of their decomposition.
- Secondary: 1 for the plain weight; then the weights of their own (finding
  9); then the mark classes in sorted order.
- Tertiary: 2 × the level within the element's class; odd values lie between
  levels. For expansions the tertiary of each element is fitted against the
  classes of its atoms: if the expansion lies between `prefix+x+min(rest)`
  and `prefix+x+max(rest)` the element takes x's level, otherwise it lies
  between levels and decides; later elements are then not observable
  against plain strings and are placed above all levels (or below all when
  the first lies below all), offset by the expansion's level within its own
  class so that, for example, `ǋ` < `Ǌ`.
- Contractions replace the pair by the elements of the equal character.
  Hangul syllables are replaced by their jamo.

Comparison: the sequence of non-zero primaries is compared first
(lexicographically, a proper prefix sorting first), then the non-zero
secondaries, then the non-zero tertiaries. Equal at all three levels is a
tie, as in R. The sort key is the three sequences joined with a 0 separator.

`C` and `POSIX` compare UTF-8 bytes (code point order), as R does in the C
locale.

## Fallback and limits

- Code points outside the probed ranges get one element with primary
  `3,368 + 1 + code point`, plain secondary and tertiary: after every probed
  character, in code point order. This matches R for CJK ideographs relative
  to the probed scripts but not among themselves (R orders ideographs
  differently from code point order: 5 of 1,289 adjacent pairs in
  `sorted_other_scripts.tsv`, all between strings with ideographs), and it
  places scripts that were not probed (Thai, Ethiopic, Georgian, ...) after
  the probed ones rather than at their ICU position.
- Contractions of three or more characters, and discontiguous matching of
  combining marks, were not probed.
- Tailored locales are rejected with an error rather than approximated.

## Verification

`crates/rcompat-collate/tests/r_orders.rs` checks every adjacent pair of the
recorded orders: ASCII strings (4,947 pairs), random strings over the probed
characters (5,622), accented Latin, Greek and Cyrillic words with combining
marks and contractions (5,414), near-duplicates and mark sequences (15,082):
all agree with R, including the 165 ties. `tests/properties.rs` checks
antisymmetry, transitivity, agreement of sort keys with `compare`, and the C
order on random strings.
