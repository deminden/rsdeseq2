# rcompat-collate

Compare and sort strings in the order R uses. R sorts character vectors and
factor levels by byte order in the `C` locale and, on builds with ICU, by
the Unicode root collation in English and most western European locales.
This crate reproduces both in pure Rust, with no data files or native
libraries.

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
  `"root"`, and locales that use the root collation (`en`, `de`, `fr`, `it`,
  `nl`, `pt` and others, with any region). Locales with their own rules, such
  as Swedish or Spanish, return an error instead of a wrong order.
- `compare`, with R's ties: strings that differ only in ignorable characters
  or in how an accent is encoded compare equal.
- `sort_key`, for sorting many strings with precomputed keys.

The weights were derived from black-box runs of R 4.6.2 with ICU 78.3 and
checked against R's order of about 32,000 strings. They cover Latin, Greek,
Cyrillic, Hebrew, Arabic, Devanagari, Hangul, kana, punctuation, symbols and
emoji; other code points sort after these in code point order. The method
and its limits are documented in the repository under
`docs/provenance/specs/collation.md`.

[API documentation](https://docs.rs/rcompat-collate) ·
[MIT license](LICENSE)
