# Specification: R's string collation (C, ICU root and tailorings)

Status: written from allowed inputs only: R's help pages (`?Comparison`,
`?sort`, `?icuSetCollate`), the general description of multi-level collation
(base letters, then accents, then case and variants), the public API of the
`unicode-normalization` crate (canonical and compatibility decomposition,
combining classes) and of Rust's `char` (case), and black-box runs of R 4.6.2
with ICU 78.3 on macOS arm64. No collation implementation, table or locale
data was consulted: no ICU or ICU4X source or data, no CLDR files, no DUCET,
no `UnicodeData.txt` or Unihan, no glibc locale files and no R source. Every
fact below is an observation of R on inputs chosen here.

## Purpose

R orders character data with the collation of the session; formula
evaluation needs that order for factor levels and string comparisons. R
reports the collation with `icuGetCollate("actual")`; the rsdeseq2 R wrapper
passes `"C"` when ICU is not in use and the reported name otherwise.

## Probes (all outputs default to `/tmp/rcollate`)

| Script | Output | Use |
|---|---|---|
| `probe_collation_full.R` | `root_single.tsv`: all 1,112,063 code points (no surrogates) in R's order with the level of each difference | derivation and oracle |
| `probe_collation_expansions.R` | `root_expansions.tsv`: groups that are several primary weights, and which | derivation |
| `probe_collation_oracle.R` | answers `compare`/`sort` questions under any locale | derivation questions, oracles |
| `probe_collation_strings.R` | `root_strings.tsv`: 200,778 random strings in R's order | oracle |
| `probe_collation_tailorings.R` | `tailoring_<locale>.tsv` (115,727 strings each), `tailoring_locales.tsv` | derivation of tailorings, oracle |
| `probe_collation_orders.R` | the committed `tests/data/sorted_*.tsv`, `locales.tsv` | crate tests |
| `probe_collation.R` | the first probe's weight table for 5,364 code points (committed observations); `probe_collation_orders.R` uses it to choose near-duplicates | historical, not used by the model |

`cargo run -p rcompat-collate --example dump_decompositions` writes
`decompositions.tsv` (NFD, NFKD, lowercase and case flags from Rust).
`build_collation_table.py` derives `crates/rcompat-collate/src/tables.rs`;
`sample_collation_orders.py` writes `oracle_summary.tsv` (counts and SHA-256
of every oracle, in `observations/collation_r462_macos`) and the committed
sample `tests/data/sorted_tailorings.tsv`. The oracles (about 30 MB) are not
committed; the scripts regenerate them.

Level tests: with LO = tab (the lowest primary) and HI = line feed, x and y
have equal primaries exactly when `x+LO+HI > y+LO+LO` and
`x+LO+LO < y+LO+HI`; with equal primaries, equal secondaries exactly when
`x+"á" > y+"a"` and `y+"á" > x+"a"`.

## Observed facts and the rules derived from them

1. 961 code points are completely ignorable (appending them changes
   nothing); 678 sort before everything with no primary weight (marks).
2. R does not reorder combining marks: `"a\u301\u325"` and `"a\u325\u301"`
   differ. Each character is replaced by its canonical decomposition, as the
   ties of every precomposed character with its decomposition show.
3. Primary groups (characters that are primary-equal): 150,381. A group is
   an expansion when an earlier atom is a primary prefix of it
   (`q <=p c <p q + U+FFFF`); everything between an atom and an expansion
   starting with it starts with it too, so one pass in sorted order finds all
   12,126 expansions, and a binary search over the atoms spells each one.
   The other 138,255 groups are atoms, each one primary weight.
4. Accents: canonical decompositions reproduce R except for 82 characters
   whose decomposition R sorts as the single character (й, Arabic letters
   with hamza or madda, Indic two-part vowels, Tibetan, Balinese and others):
   these are contractions. R also ties `l`+U+00B7 with ŀ (and `L`).
5. Compatibility decompositions (ligatures, fullwidth, super- and subscripts,
   circled, squared, fractions, Arabic presentation forms) reproduce R's
   primaries except for spacing accents and a few letters, which are stored.
   Their tertiary weight is a variant weight per character range for the
   lowercase (or uncased) and uppercase letters they contain, read off each
   character's place among the variants of its letter; for characters of
   several letters the first letter is fitted against the variants of that
   letter (with x_j the first letter replaced by the variant of layer j, R
   puts the character between x_j and x_(j+1); comparing with the highest
   variant of the second letter separates "equal to layer j" from "between").
6. Case: an uppercase letter sorts as its lowercase letter (Rust's
   `to_lowercase`) with the uppercase tertiary weight, lowercase first.
7. Variants R ties or ranks with a letter that no rule explains (digits of
   all scripts with ASCII digits, kana, ð, ł, ø, đ, ſ, æ, ß, œ, ...) are
   stored: 2,344 characters, mostly as runs. A class of a group that has no
   rule-explained member is either the atom plus one mark (đ = d + U+0335,
   found by asking R against every mark) or a secondary weight of its own
   above the common one (ask R whether the class is secondary-equal to its
   spelled atoms); tertiary weights are a longest-path layering of the
   observed order of classes.
8. Ideographs are not in code point order: R orders them by radical and
   stroke, with URO, Extension A, the unified compatibility ideographs and
   Extensions B to J interleaved. Within each block the order is code point
   order apart from 9,600 characters; the data are runs per block plus those
   characters (64 kB). Tangut, Khitan and Nüshu are atoms in code point order
   before the ideographs. Hangul syllables tie with their jamo (canonical
   decomposition).
9. Unassigned code points and private use sort after everything in code point
   order, then U+FFFD and U+FFFF. Code points that the Unicode data of the
   pinned `unicode-normalization` 0.1.25 or of Rust's `char` know but that
   are unassigned for this ICU would be listed in `UNASSIGNED`; there are
   none (with 0.1.26 there were 68, and 0.1.26 also adds the Unicode-3.0
   licence, so the dependency is pinned).
10. Prevowels: Thai U+0E40..0E44, Lao U+0EC0..0EC4, Tai Viet and New Tai Lue
    prevowels followed by one of their consonants tie with the swapped pair
    (found by the ignorable-insertion test over neighbouring atoms, 2.26 M
    pairs, which found no other contractions).

## Model

Element = (primary, secondary, tertiary); a zero weight is absent at that
level. Primary of the i-th atom = (i + 1) · 256 (gaps for tailorings);
secondary 32 = common, then weights of their own, then mark classes in R's
order; tertiary codes · 32. A code point maps by the first rule that applies:
unassigned list → stored → atom → mark → ignorable → compatibility
decomposition with variant tertiary → lowercase with uppercase tertiary →
code point order after all atoms. Strings compare by primaries, then
secondaries, then tertiaries (a proper prefix first); the sort key joins the
three sequences with 0. `C`/`POSIX` compare UTF-8 bytes.

## Tailorings

`icuSetCollate(locale = ...)` changes the order (`tailoring_locales.tsv` lists
what ICU reports; `ca`, `de`, `en`, `fr`, `it`, `nl`, `pt` report `root`,
`nb`/`nn` report `no`). Each tailoring is a rule list on top of root
(`src/tailoring/rules.rs`): `&x < y` (new letter after x), `<<`, `<<<`, `=`,
`&[before 1]`, `&[before 2]`, contractions (Czech `ch`, Hungarian `cs`, `dzs`,
...), `x/y` (Hungarian `ccs` = `cs` + `cs`, Norwegian `þ` = `t` + `h`),
script moves (Cyrillic first for ru, uk, be, bg, sr, mk; Cyrillic after Latin
for hr, bs, az) and switching off a root contraction (й in sr, mk). They were
read off the differences between each locale's order and root, and a tailored
base+mark also matches with marks of lower combining class in between
(observed: `ç̌` is `č` with a cedilla in Czech). Rules were then corrected
until R's order of the 115,727 probe strings was reproduced exactly.

## Verification (`cargo run --release -p rcompat-collate --example verify`)

| Oracle | Adjacent pairs | Differ |
|---|---|---|
| every code point (`root_single.tsv`), with levels | 1,111,101 | 0 |
| random strings (`root_strings.tsv`), with levels | 200,777 | 0 |
| committed `tests/data/sorted_*.tsv` (old probe) | 32,354 | 0 |
| each of sv, fi, nb, nn, es, gl, pl, cs, sk, hu, hr, bs, sl, ro, lt, lv, et, tr, az, ru, uk, be, bg, sr, mk, de@collation=phonebook | 115,726 each | 0 |

Not supported (error): da (31 pairs still differed: uppercase-first details),
is (secondary order of marks changes), el (34: secondary order of some
marks), fr_CA (backwards accents and ŀ), and every other tailored locale.

## Limits

- Contractions of root were searched among neighbouring atoms and pairs with
  marks; contractions between distant characters were not probed.
- Discontiguous matching was observed and implemented for tailored letters;
  root contractions match contiguously, which reproduced every oracle.
- The decomposition and case facts come from `unicode-normalization` and
  Rust's `char`; a newer Unicode version there can assign code points ICU 78
  treats as unassigned (rerun the derivation and the verification).
