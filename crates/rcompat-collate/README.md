# rcompat-collate

Compare and sort strings in the order R uses. R sorts character vectors and
factor levels by byte order in the `C` locale and, on builds with ICU, by
the Unicode root collation in English and most western European locales.
Locales with rules of their own, such as Swedish or Polish, tailor that
order. This crate reproduces these in pure Rust, with no data files or
native libraries.

## Install

Requires Rust 1.96 or newer.

```sh
cargo add rcompat-collate@0.1.0
```

## Example

Sort factor levels as `sort()` does in an `en_US` R session, and as it does
under `LC_COLLATE=C`.

```rust
use rcompat_collate::Collator;

fn main() -> Result<(), rcompat_collate::Error> {
    let levels = ["treated", "Control", "control", "\u{e9}tat", "etat", "_ref", "10", "9"];

    let root = Collator::new("en_US")?;
    let mut sorted = levels.to_vec();
    sorted.sort_by(|a, b| root.compare(a, b));
    assert_eq!(sorted, ["_ref", "10", "9", "control", "Control", "etat", "\u{e9}tat", "treated"]);

    let c = Collator::new("C")?;
    sorted.sort_by_key(|s| c.sort_key(s));
    assert_eq!(sorted, ["10", "9", "Control", "_ref", "control", "etat", "treated", "\u{e9}tat"]);
    Ok(())
}
```

## What it provides

- `Collator::new` for the names R reports: `"C"`/`"POSIX"` (byte order),
  `"root"`, the locales that use the root collation (`en`, `de`, `fr`, `it`,
  `nl`, `pt` and others, with any region), and 29 tailorings listed in
  `SUPPORTED_LOCALES`: `sv`, `da`, `nb`/`nn`/`no`, `fi`, `is`, `es`, `gl`,
  `pl`, `cs`, `sk`, `hu`, `hr`, `bs`, `sl`, `ro`, `lt`, `lv`, `et`, `tr`,
  `az`, `ru`, `uk`, `be`, `bg`, `sr`, `mk`, `el`, `fr_CA` and
  `de@collation=phonebook`, including Danish uppercase-first order, Canadian
  French accents compared from the end of the string and Greek before Latin.
  Other locales with rules of their own, whose rules were not derived
  (Welsh, Esperanto, Vietnamese, traditional Spanish, ...), return an error
  instead of a wrong order.
- `compare` and `compare_with_strength`, with R's ties: strings that differ
  only in ignorable characters or in how an accent is encoded compare equal.
- `sort_key`, for sorting many strings with precomputed keys.

The root collation is a set of rules (canonical and compatibility
decomposition, case, code point ranges, prevowels) plus about 62 kB of
data for what the rules cannot produce, 27 kB of it the radical-and-stroke
order of the ideographs. Everything was derived from black-box runs of R
4.6.2 with ICU 78.3 and checked against R: under the root collation and
under each tailoring, the order of every one of the 1,112,063 code points
and of 25,172 strings with ligatures and other expansions; under root also
200,778 random strings, and under each tailoring 115,727 strings of
accented letters, with no difference. Comparing two short ASCII strings
takes about 15 ns. The method and its limits are documented in the
repository under `docs/provenance/specs/collation.md`.

[API documentation](https://docs.rs/rcompat-collate) ·
[MIT license](LICENSE)
