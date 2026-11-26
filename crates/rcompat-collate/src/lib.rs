//! String collation that orders text the way R does.
//!
//! Part of the rcompat series. R compares and sorts character data, and
//! orders factor levels, by the collation of the session: byte order in the
//! `C` locale, and otherwise (on builds with ICU) the Unicode root collation
//! for English and most other western locales. [`Collator`] reproduces both:
//!
//! - `"C"` and `"POSIX"` compare UTF-8 bytes, which is code point order.
//! - `"root"`, `"und"`, `""` and the locales that ICU maps to the root
//!   collation (see [`Collator::new`]) compare with three levels: base
//!   letters and symbols first, then accents, then case and letter variants.
//!
//! The weights were derived from black-box runs of R 4.6.2 with ICU 78.3; the
//! method, the observations and the limits of the model are documented in the
//! repository under `docs/provenance/specs/collation.md`.
//!
//! ```
//! use rcompat_collate::Collator;
//!
//! let root = Collator::new("root")?;
//! let mut words = vec!["banana", "Apple", "apple", "\u{e9}clair", "eclair", "_x", "10", "9"];
//! words.sort_by(|a, b| root.compare(a, b));
//! assert_eq!(words, ["_x", "10", "9", "apple", "Apple", "banana", "eclair", "\u{e9}clair"]);
//! # Ok::<(), rcompat_collate::Error>(())
//! ```

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code)]

mod locale;
mod table;

// Compile and run the README example as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

use std::cmp::Ordering;
use std::fmt;

/// Error returned for a collation name that is not supported.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// The locale has collation rules of its own (a tailoring), or is not a
    /// locale identifier this crate recognises.
    UnsupportedLocale(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedLocale(name) => write!(
                f,
                "collation '{name}' is not supported; use \"C\", \"root\" or a locale that uses \
                 the root collation (such as \"en\", \"en_US\", \"de\" or \"fr\")"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// Compares strings in R's order for one collation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Collator {
    kind: Kind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Bytes,
    Root,
}

/// A precomputed comparison key: comparing keys gives the same result as
/// [`Collator::compare`] on the strings they were made from.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SortKey(Vec<u32>);

impl Collator {
    /// Collator for a collation name as R or ICU report it.
    ///
    /// Accepted names, without regard to letter case and with `_` or `-`
    /// between subtags:
    ///
    /// - `"C"` and `"POSIX"`: byte order;
    /// - `"root"`, `"und"` and `""`: the root collation;
    /// - a language, optionally with a region, for which ICU uses the root
    ///   collation: `en`, `de`, `fr`, `it`, `nl`, `pt`, `ca`, `eu`, `ga`,
    ///   `id`, `ms` and `sw` (for example `"en_US"` or `"de-AT"`), except
    ///   `fr_CA`.
    ///
    /// Other locales tailor the root order (Swedish, Danish, Spanish, Polish,
    /// Russian and many more) and return [`Error::UnsupportedLocale`], as do
    /// identifiers with variants or extensions.
    pub fn new(locale: &str) -> Result<Self, Error> {
        match locale::classify(locale) {
            Some(locale::Collation::Bytes) => Ok(Self::c()),
            Some(locale::Collation::Root) => Ok(Self::root()),
            None => Err(Error::UnsupportedLocale(locale.to_owned())),
        }
    }

    /// The root collation.
    pub fn root() -> Self {
        Self { kind: Kind::Root }
    }

    /// Byte (code point) order, as in R's `C` locale.
    pub fn c() -> Self {
        Self { kind: Kind::Bytes }
    }

    /// Whether this collator compares bytes.
    pub fn is_bytewise(&self) -> bool {
        self.kind == Kind::Bytes
    }

    /// Compare two strings. `Equal` means R treats them as tied, which under
    /// the root collation also happens for distinct strings that differ only
    /// in ignorable characters or in the encoding of an accent (precomposed
    /// or combining).
    pub fn compare(&self, left: &str, right: &str) -> Ordering {
        match self.kind {
            Kind::Bytes => left.as_bytes().cmp(right.as_bytes()),
            Kind::Root => {
                let a = elements(left);
                let b = elements(right);
                for level in [Level::Primary, Level::Secondary, Level::Tertiary] {
                    let order = weights(&a, level).cmp(weights(&b, level));
                    if order != Ordering::Equal {
                        return order;
                    }
                }
                Ordering::Equal
            }
        }
    }

    /// Key whose order equals [`Collator::compare`], for sorting many strings.
    pub fn sort_key(&self, text: &str) -> SortKey {
        match self.kind {
            Kind::Bytes => SortKey(text.bytes().map(u32::from).collect()),
            Kind::Root => {
                let e = elements(text);
                let mut key = Vec::with_capacity(3 * e.len() + 2);
                for (i, level) in [Level::Primary, Level::Secondary, Level::Tertiary]
                    .into_iter()
                    .enumerate()
                {
                    if i > 0 {
                        // Weights are positive, so the separator ends the
                        // shorter sequence first.
                        key.push(0);
                    }
                    key.extend(weights(&e, level));
                }
                SortKey(key)
            }
        }
    }
}

impl Default for Collator {
    /// The root collation.
    fn default() -> Self {
        Self::root()
    }
}

#[derive(Clone, Copy)]
enum Level {
    Primary,
    Secondary,
    Tertiary,
}

/// One collation element: primary, secondary and tertiary weight. Zero means
/// the element has no weight at that level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Element {
    primary: u32,
    secondary: u8,
    tertiary: u8,
}

const COMMON_SECONDARY: u8 = 1;
const PLAIN_TERTIARY: u8 = 2;

fn weights(elements: &[Element], level: Level) -> impl Iterator<Item = u32> + '_ {
    elements
        .iter()
        .map(move |e| match level {
            Level::Primary => e.primary,
            Level::Secondary => u32::from(e.secondary),
            Level::Tertiary => u32::from(e.tertiary),
        })
        .filter(|w| *w != 0)
}

fn unpack(value: u32) -> Element {
    Element {
        primary: value >> 16,
        secondary: ((value >> 8) & 0xFF) as u8,
        tertiary: (value & 0xFF) as u8,
    }
}

/// Push the elements of one code point.
fn push_code_point(code_point: u32, out: &mut Vec<Element>) {
    match table::CODE_POINTS.binary_search(&code_point) {
        Ok(index) => {
            let start = usize::from(table::OFFSETS[index]);
            let end = usize::from(table::OFFSETS[index + 1]);
            out.extend(table::ELEMENTS[start..end].iter().map(|v| unpack(*v)));
        }
        // Code points outside the observed set sort after every observed
        // primary, in code point order, with plain secondary and tertiary.
        Err(_) => out.push(Element {
            primary: table::PRIMARY_COUNT + 1 + code_point,
            secondary: COMMON_SECONDARY,
            tertiary: PLAIN_TERTIARY,
        }),
    }
}

fn contraction(first: u32, second: u32) -> Option<u32> {
    table::CONTRACTIONS
        .iter()
        .find(|(a, b, _)| *a == first && *b == second)
        .map(|(_, _, target)| *target)
}

/// First Hangul syllable and the bases of the conjoining jamo it is spelled
/// with (Unicode Standard, section 3.12).
const SYLLABLE_BASE: u32 = 0xAC00;
const SYLLABLE_COUNT: u32 = 11_172;
const LEADING_BASE: u32 = 0x1100;
const VOWEL_BASE: u32 = 0x1161;
const TRAILING_BASE: u32 = 0x11A7;
const VOWEL_COUNT: u32 = 21;
const TRAILING_COUNT: u32 = 28;

fn elements(text: &str) -> Vec<Element> {
    let mut out = Vec::with_capacity(text.len());
    let mut chars = text.chars().map(u32::from).peekable();
    while let Some(c) = chars.next() {
        if let Some(target) = chars.peek().and_then(|next| contraction(c, *next)) {
            chars.next();
            push_code_point(target, &mut out);
        } else if (SYLLABLE_BASE..SYLLABLE_BASE + SYLLABLE_COUNT).contains(&c) {
            // R ties every syllable with its jamo sequence.
            let index = c - SYLLABLE_BASE;
            push_code_point(
                LEADING_BASE + index / (VOWEL_COUNT * TRAILING_COUNT),
                &mut out,
            );
            push_code_point(
                VOWEL_BASE + (index % (VOWEL_COUNT * TRAILING_COUNT)) / TRAILING_COUNT,
                &mut out,
            );
            if !index.is_multiple_of(TRAILING_COUNT) {
                push_code_point(TRAILING_BASE + index % TRAILING_COUNT, &mut out);
            }
        } else {
            push_code_point(c, &mut out);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_consistent() {
        assert_eq!(table::OFFSETS.len(), table::CODE_POINTS.len() + 1);
        assert!(table::CODE_POINTS.windows(2).all(|w| w[0] < w[1]));
        assert!(table::OFFSETS.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(
            usize::from(table::OFFSETS[table::OFFSETS.len() - 1]),
            table::ELEMENTS.len()
        );
        for value in table::ELEMENTS {
            let e = unpack(value);
            assert!(e.primary <= table::PRIMARY_COUNT);
            assert!(e.secondary >= COMMON_SECONDARY && e.tertiary > 0);
        }
    }

    #[test]
    fn fallback_follows_every_observed_primary() {
        let root = Collator::root();
        assert_eq!(root.compare("\u{4e00}", "\u{ff5e}"), Ordering::Greater);
        assert_eq!(root.compare("\u{4e00}", "\u{4e01}"), Ordering::Less);
    }
}
