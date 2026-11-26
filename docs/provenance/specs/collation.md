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
| `probe_collation_full.R` | `root_single.tsv`: all 1,112,063 code points (no surrogates) in R's order with the level of each difference; with a locale, `<locale>_single.tsv` under that collation | derivation and oracle |
| `probe_collation_expansions.R` | `root_expansions.tsv`: groups that are several primary weights, and which | derivation |
| `probe_collation_oracle.R` | answers `compare`/`sort` questions under any locale | derivation questions, oracles |
| `probe_collation_strings.R` | `root_strings.tsv`: 200,778 random strings in R's order | oracle |
| `probe_collation_tailorings.R` | `tailoring_<locale>.tsv` (115,727 strings each), `tailoring_locales.tsv` | derivation of tailorings, oracle |
| `collation_expansion_strings.py` | `expansion_input.tsv`: 25,172 strings, every character R sorts as several atoms with its spelling in lowercase, uppercase and capitalised atoms and its compatibility decomposition | oracle input |
| `probe_collation_supported.R` | for every supported tailoring `<locale>_single.tsv` (every code point, through `probe_collation_full.R`) and `expansion_<locale>.tsv`, and `expansion_root.tsv` | derivation of tailorings, oracle |
| `probe_collation_orders.R` | the committed `tests/data/sorted_*.tsv`, `locales.tsv` | crate tests |
| `probe_collation.R` | the first probe's weight table for 5,364 code points (committed observations); `probe_collation_orders.R` uses it to choose near-duplicates | historical, not used by the model |

`cargo run -p rcompat-collate --example dump_decompositions` writes
`decompositions.tsv` (NFD, NFKD, lowercase and case flags from Rust).
`build_collation_table.py` derives `crates/rcompat-collate/src/tables.rs`;
`sample_collation_orders.py` writes `oracle_summary.tsv` (counts and SHA-256
of every oracle, in `observations/collation_r462_macos`) and the committed
sample `tests/data/sorted_tailorings.tsv`. The oracles (about 450 MB) are
not committed; the scripts regenerate them.

Level tests: with LO = tab (the lowest primary) and HI = line feed, x and y
have equal primaries exactly when `x+LO+HI > y+LO+LO` and
`x+LO+LO < y+LO+HI`; with equal primaries, equal secondaries exactly when
`x+"á" > y+"a"` and `y+"á" > x+"a"` (root, `probe_collation_full.R` and
`probe_collation_strings.R`). Under a tailoring that test can mislead: in
Danish and Norwegian `x` ending in `a` forms the contraction `aa` with the
appended letter, and in Icelandic `á` is a letter of its own. The tailored
oracles therefore keep the accented letter apart with a tab
(`x+LO+acc > y+LO+"a"` both ways), take as `acc` the first of á, à, ă, ā
that the collation sorts as a secondary variant of `a`, and put it before
`x` (`acc+LO+x` against `"a"+LO+y`) when the collation compares accents from
the end, as `fr_CA` does; the header of each oracle records the test used.
The earlier test had labelled `ª` a secondary variant in Norwegian (it is a
tertiary one, as in root), which a rule had followed; that rule is gone.

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
   characters ("moved"; 7,930 once the longest increasing subsequence of
   each block is kept in place). The moved characters of a block are
   themselves mostly consecutive (6,287 follow the previous moved character
   of their block), so the order is stored as a byte stream of 27,025 bytes
   (27 kB with the block list, 64 kB before): each entry names a block in its
   high nibble and takes the next 1 to 269 characters of the block's main
   sequence, the next moved character, or moved characters at a signed
   distance from where the block's moved sequence stopped
   (`tables::HAN_ORDER` documents the format; it is decoded once). The
   characters are those of the earlier encoding in the same order, so every
   comparison is unchanged. Tangut, Khitan and Nüshu are atoms in code point order
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
11. Expansions against their spelling (`expansion_root.tsv`): the squared
    katakana words with voiced sound marks (42 characters from U+3300) are
    secondary-equal to their compatibility decomposition (R reports a
    tertiary difference), not to the atoms alone; the derivation asks R for
    the decomposition (when it contains marks) and then for the atoms with
    one mark after each of them, and stores such a character as that
    sequence.
12. Some differences in root have two representations that every root
    order agrees with; tailored orders tell them apart. Under the backwards
    accents of `fr_CA`, U+0140 (l with middle dot) and U+013F sort after `l`
    followed by any mark: they are `l` followed by a secondary weight above
    every mark (reference 0 in `ELEMENTS`), where the other characters with
    a weight of their own (ð, ...) have it below the marks; the derivation
    asks R under `fr_CA` for each such character. Under Danish, U+1F18C
    (squared "PA") sorts before U+33A9 (squared "Pa") although root puts it
    after: their first letters have the same tertiary weight and only the
    second differs (`SAME_FIRST_TERTIARY` in the derivation).
13. Case, which only uppercase-first collations use (`da_single.tsv`,
    `expansion_da.tsv`): an element is uppercase when its character is an
    uppercase letter or has a lowercase form (Rust's `char`), and the
    elements of a compatibility decomposition keep the case of their
    letters (Ⓐ, ᴬ, ㏇ uppercase; ª, ᵃ lowercase). Exceptions observed in R:
    normal-size kana (hiragana, katakana and their halfwidth forms) and the
    halfwidth forms U+FF61..FFDC and U+FFE8..FFEE are uppercase, the small
    kana (28 ranges listed in `root.rs`) are
    lowercase; U+20A7 (peseta sign), U+214D and the squared words
    U+1F18A..1F19A and U+1F1A6..1F1AD are uppercase; the superscript signs
    U+2120, U+2122 and U+1F16A..1F16C are lowercase although they decompose
    to capitals. A tailored item gets the case of its characters: all
    uppercase, all lowercase, or mixed (Danish `Aa`, sorted between `AA` and
    `å`).

## Model

Element = (primary, secondary, tertiary, case); a zero weight is absent at
that level, and the case (lowercase, mixed, uppercase) counts only under
uppercase-first collations. Primary of the i-th atom = (i + 1) · 256 (gaps for tailorings);
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
for hr, bs, az; Greek, from α up to Coptic, before Latin for el) and
switching off a root contraction (й in sr, mk). They were read off the
differences between each locale's order and root, and a tailored base+mark
also matches with marks of lower combining class in between (observed:
`ç̌` is `č` with a cedilla in Czech; in Icelandic `a\u301\u325` ties with
`a\u325\u301`). Three options complete the rules, each observed in R:

- Uppercase first (`da`): at the tertiary level each element compares by
  case first, uppercase before mixed before lowercase, then by weight.
  Observed: `A < ᴬ < a < ª`, `TH < Th < Þ < th < þ`, `Å < AA < Aa < å < aa`,
  `あ < ア < ｱ < ぁ < ァ < ｧ` (root: `a < A`, `ぁ < あ < ァ < ｧ < ア < ｱ`);
  the case of each character is fact 13.
- Backwards accents (`fr_CA`): secondary weights compare from the end of
  the string. Observed: `cote < côte < coté < côté` (root:
  `cote < coté < côte < côté`) and `a\u325\u301 < a\u301\u325` (root: the
  reverse).
- Canonical order of marks (`el`): combining marks are put in canonical
  order before collation, so `a\u301\u325` ties with `a\u325\u301` and
  `\u323\u327` with `\u327\u323` (root and the other locales keep them
  apart).

The rules of the four locales added last, and corrections to earlier ones,
came from the orders of every code point (`<locale>_single.tsv`):

- `da`: like Norwegian, but `ę` stays an `e` and `œ` is `oe` with an accent
  (`oé < œ < o̲e`); `aa` is a tertiary variant of `å` (`å <<< aa`); `ª`
  stays a tertiary variant of `a`; uppercase first.
- `is`: á, é, í, ó, ú, ý are letters after the last form of a, e, i, o, u, y
  (`&[before 1] b < á`, ...); `đ` sorts after `d` with any mark (`&d << đ`,
  `d͡ < đ`); `ð` is a letter after `d`; æ (with ä), ö (with ø) and å follow
  the Latin letters.
- `el`: Greek before Latin and canonical order of marks; no rules.
- `fr_CA`: backwards accents; no rules (`ŀ` is fact 12).
- `sv`: `ð` is an accent variant of `d` after `đ` (`ð << ᷙ`, root
  `ð = ᷙ`).
- `fi`: `ŋ`, `ǥ`, `ŧ`, `ʒ` sort right after their letter with U+0335
  (`n̵ < ŋ < n̅`).
- `et`: `õ` follows every form of `w` (`ʍ < ꭩ < õ`), not `w` itself.
- `uk`: `ї` follows Cyrillic iota (`ꙇ < ї`), with U+A676 as its tertiary
  variant (`ї <<< ꙶ <<< Ї`).

Rules were corrected until R's order of the 115,727 probe strings, of every
code point and of the 25,172 expansion strings was reproduced exactly.

## Verification (`cargo run --release -p rcompat-collate --example verify -- /tmp/rcollate`)

| Oracle | Adjacent pairs | Differ |
|---|---|---|
| every code point (`root_single.tsv`), with levels | 1,111,101 | 0 |
| random strings (`root_strings.tsv`), with levels | 200,777 | 0 |
| expansions under root (`expansion_root.tsv`), with levels | 25,171 | 0 |
| committed `tests/data/sorted_*.tsv` (old probe) | 32,354 | 0 |
| each of sv, da, nb, nn, fi, is, es, gl, pl, cs, sk, hu, hr, bs, sl, ro, lt, lv, et, tr, az, ru, uk, be, bg, sr, mk, el, fr_CA, de@collation=phonebook: probe strings, every code point, expansions, with levels | 115,726 + 1,111,101 + 25,171 each | 0 |

Not supported (error): the tailored locales whose rules were not derived
(`sq`, `mt`, `cy`, `af`, `fo`, `kl`, `se`, `es@collation=traditional`,
`sr_Latn`, `en_US_POSIX`, `vi`, `ha`, `yo`, `eo` among those probed).

## Limits

- Contractions of root were searched among neighbouring atoms and pairs with
  marks; contractions between distant characters were not probed.
- Discontiguous matching was observed and implemented for tailored letters;
  root contractions match contiguously, which reproduced every oracle.
- Case (fact 13) was observed under Danish, the only uppercase-first
  collation supported, through every code point and the expansions; a
  character whose case no such order shows keeps the rule's case.
- The decomposition and case facts come from `unicode-normalization` and
  Rust's `char`; a newer Unicode version there can assign code points ICU 78
  treats as unassigned (rerun the derivation and the verification).
