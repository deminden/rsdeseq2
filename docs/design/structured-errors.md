# Design note 2.1 — structured error variants

Status: **implemented in step 2.1b** (2026-11-18), on Denis's request
"do error redesign".

## Problem

`DeseqError` has 13 variants. Five of them (`InvalidOptions`,
`InvalidDispersion`, `InvalidSizeFactors`, `InvalidCounts`,
`UnsupportedFeature`) carry one free-text `reason`. There are about 700 such
messages, and 417 of them are `InvalidOptions`. A caller can only tell
"bad parameter value", "bad formula", "unknown coefficient" and "internal bug"
apart by parsing English. The crate does this itself in one place:
`contrasts::aliases` matches on a message prefix to decide whether a
coefficient search failed.

## Constraints

- **Messages stay byte-identical.** 190 R `expect_error` checks and 158 Rust
  `.contains(...)` checks match message text, and so may users' scripts.
  Every new variant formats exactly the text its sites produced before.
  This makes the change IDENTICAL for numeric output and for messages.
- **`ErrorKind` stays the coarse, stable classification.** Each new variant
  maps to the kind its sites had, except internal invariants, which get their
  own kind. `ErrorKind` itself is unreleased, so adding kinds costs nothing.
- `DeseqError` is `#[non_exhaustive]`, so adding variants is not a breaking
  change. Moving existing errors out of `InvalidOptions` is a behaviour
  change for code that matches that variant; the CHANGELOG says so.

## Variants added

| Variant | Replaces | Display (unchanged) | `kind()` |
| --- | --- | --- | --- |
| `InvalidParameter { domain, parameter, requirement }` | "*name* must be finite / positive / non-negative / finite and positive / finite and non-negative / non-empty / finite and within [0, 1]" in `InvalidOptions`, `InvalidDispersion`, `InvalidSizeFactors`, `InvalidCounts` | `"{domain prefix}: {parameter} must be {requirement}"` | `domain` |
| `InvalidFormula { reason }` | `InvalidOptions` whose message starts with `formula` | `"invalid options: {reason}"` | `Formula` (new) |
| `CoefficientNotFound { reason, candidates }` | the "not present in coefficient names" family | `"invalid options: {reason}"` | `Options` |
| `Internal { what }` | `invariant(..)`, which returned `InvalidOptions` | `"invalid options: internal invariant violated: {what}"` | `Internal` (new) |

`Requirement` is a public `#[non_exhaustive]` enum whose `Display` is the
exact phrase. `domain` is an `ErrorKind` limited to `Options`, `Dispersion`,
`Normalization` and `Counts`; its prefix is the old variant's prefix.

## Not in this step

- The ~150 "*x* requires *y*" messages: their wording varies too much to
  format from fields without changing text. A later step can add
  `MissingInput { what, needed_for }` together with message normalisation,
  which would need Denis's sign-off because it changes messages.
- `scoped_invalid_options` sites keep their `"{scope}: {reason}"` text.
- The free-text variants stay; nothing is deprecated.
