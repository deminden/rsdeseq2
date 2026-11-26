//! String collation that orders text the way R does.
//!
//! Part of the rcompat series. R compares and sorts character data, and
//! orders factor levels, by the collation of the session: byte order in the
//! `C` locale and, on builds with ICU, the Unicode root collation or a
//! language's tailoring of it. [`Collator`] reproduces these:
//!
//! - `"C"` and `"POSIX"` compare UTF-8 bytes, which is code point order.
//! - `"root"`, `"und"`, `""` and the locales that ICU maps to the root
//!   collation compare with three levels: base letters and symbols first,
//!   then accents, then case and letter variants.
//! - Swedish, Spanish, Polish, Russian, Turkish and the other locales listed
//!   in [`SUPPORTED_LOCALES`] apply their rules on top of the root collation.
//!
//! The root collation is a set of rules (canonical and compatibility
//! decomposition, case, code point ranges) plus about 100 kB of data for
//! what the rules cannot produce, chiefly the radical-and-stroke order of the
//! ideographs. Everything was derived from black-box runs of R 4.6.2 with ICU
//! 78.3; the method, the observations and the verification are documented in
//! the repository under `docs/provenance/specs/collation.md`.
//!
//! ```
//! use rcompat_collate::Collator;
//!
//! let root = Collator::new("root")?;
//! let mut words = vec!["banana", "Apple", "apple", "\u{e9}clair", "eclair", "_x", "10", "9"];
//! words.sort_by(|a, b| root.compare(a, b));
//! assert_eq!(words, ["_x", "10", "9", "apple", "Apple", "banana", "eclair", "\u{e9}clair"]);
//!
//! let swedish = Collator::new("sv_SE")?;
//! let mut names = vec!["\u{f6}rn", "zebra", "\u{e5}sa", "anna"];
//! names.sort_by(|a, b| swedish.compare(a, b));
//! assert_eq!(names, ["anna", "zebra", "\u{e5}sa", "\u{f6}rn"]);
//! # Ok::<(), rcompat_collate::Error>(())
//! ```

#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]
#![deny(unsafe_code)]

mod locale;
mod root;
mod tables;
mod tailoring;

// Compile and run the README example as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

use std::cmp::Ordering;
use std::fmt;

use root::{Element, root};
use tailoring::Tailoring;

pub use locale::SUPPORTED_LOCALES;

/// Error returned for a collation name that is not supported.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// The name is not a locale identifier this crate recognises, or the
    /// locale has collation rules this crate does not implement.
    UnsupportedLocale(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedLocale(name) => write!(
                f,
                "collation '{name}' is not supported; use \"C\", \"root\" or one of the locales \
                 in rcompat_collate::SUPPORTED_LOCALES"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// How many levels of difference a comparison looks at.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strength {
    /// Base letters, digits and symbols only.
    Primary,
    /// Also accents.
    Secondary,
    /// Also case and letter variants (the strength R uses).
    Tertiary,
}

/// Compares strings in R's order for one collation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Collator {
    kind: Kind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Bytes,
    Rules(Option<Tailoring>),
}

/// A precomputed comparison key: comparing keys gives the same result as
/// [`Collator::compare`] on the strings they were made from.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SortKey(Vec<u32>);

impl Collator {
    /// Collator for a collation name as R or ICU report it
    /// (`icuGetCollate("actual")`), or for a locale identifier.
    ///
    /// Names are matched without regard to letter case, with `_` or `-`
    /// between subtags and an optional region (`"sv_SE"`, `"pt-BR"`):
    ///
    /// - `"C"` and `"POSIX"`: byte order;
    /// - `"root"`, `"und"` and `""`: the root collation;
    /// - languages for which ICU uses the root collation: `en`, `de`, `fr`,
    ///   `it`, `nl`, `pt`, `ca`, `eu`, `ga`, `id`, `ms` and `sw`;
    /// - languages with rules of their own, listed in [`SUPPORTED_LOCALES`].
    ///
    /// Any other name returns [`Error::UnsupportedLocale`] rather than a
    /// silently different order.
    pub fn new(locale: &str) -> Result<Self, Error> {
        match locale::classify(locale) {
            Some(locale::Collation::Bytes) => Ok(Self::c()),
            Some(locale::Collation::Root) => Ok(Self::root()),
            Some(locale::Collation::Tailored(t)) => Ok(Self {
                kind: Kind::Rules(Some(t)),
            }),
            None => Err(Error::UnsupportedLocale(locale.to_owned())),
        }
    }

    /// The root collation.
    pub fn root() -> Self {
        Self {
            kind: Kind::Rules(None),
        }
    }

    /// Byte (code point) order, as in R's `C` locale.
    pub fn c() -> Self {
        Self { kind: Kind::Bytes }
    }

    /// Whether this collator compares bytes.
    pub fn is_bytewise(&self) -> bool {
        self.kind == Kind::Bytes
    }

    /// The name of the collation as ICU reports it: `"C"`, `"root"`, or the
    /// tailoring (`"sv"`, `"no"`, `"de@collation=phonebook"`, ...).
    pub fn name(&self) -> &'static str {
        match self.kind {
            Kind::Bytes => "C",
            Kind::Rules(None) => "root",
            Kind::Rules(Some(t)) => t.name(),
        }
    }

    /// Compare two strings. `Equal` means R treats them as tied, which under
    /// the root collation also happens for distinct strings that differ only
    /// in ignorable characters or in the encoding of an accent (precomposed
    /// or combining).
    pub fn compare(&self, left: &str, right: &str) -> Ordering {
        self.compare_with_strength(left, right, Strength::Tertiary)
    }

    /// Compare two strings looking only at the first `strength` levels. Under
    /// `"C"` the strength is ignored.
    pub fn compare_with_strength(&self, left: &str, right: &str, strength: Strength) -> Ordering {
        let tailoring = match self.kind {
            Kind::Bytes => return left.as_bytes().cmp(right.as_bytes()),
            Kind::Rules(t) => t,
        };
        if left.is_ascii() && right.is_ascii() && tailoring.is_none_or(|t| !t.touches_ascii()) {
            return compare_ascii(left.as_bytes(), right.as_bytes(), strength);
        }
        let a = elements(tailoring, left);
        let b = elements(tailoring, right);
        compare_elements(&a, &b, strength)
    }

    /// Key whose order equals [`Collator::compare`], for sorting many strings.
    pub fn sort_key(&self, text: &str) -> SortKey {
        let tailoring = match self.kind {
            Kind::Bytes => return SortKey(text.bytes().map(u32::from).collect()),
            Kind::Rules(t) => t,
        };
        let e = elements(tailoring, text);
        let mut key = Vec::with_capacity(3 * e.len() + 2);
        key.extend(e.iter().map(|e| e.primary).filter(|w| *w != 0));
        // Weights are positive, so the separator ends the shorter sequence
        // first.
        key.push(0);
        key.extend(secondaries(&e).map(u32::from));
        key.push(0);
        key.extend(tertiaries(&e).map(u32::from));
        SortKey(key)
    }
}

impl Default for Collator {
    /// The root collation.
    fn default() -> Self {
        Self::root()
    }
}

fn elements(tailoring: Option<Tailoring>, text: &str) -> Vec<Element> {
    let seq = root::decompose(text);
    let mut out = Vec::with_capacity(seq.len());
    match tailoring {
        None => root().push_decomposed(&seq, &mut out),
        Some(t) => t.push_decomposed(&seq, &mut out),
    }
    out
}

fn primaries(e: &[Element]) -> impl Iterator<Item = u32> + '_ {
    e.iter().map(|e| e.primary).filter(|w| *w != 0)
}

fn secondaries(e: &[Element]) -> impl Iterator<Item = u16> + '_ {
    e.iter().map(|e| e.secondary).filter(|w| *w != 0)
}

fn tertiaries(e: &[Element]) -> impl Iterator<Item = u16> + '_ {
    e.iter().map(|e| e.tertiary).filter(|w| *w != 0)
}

fn compare_elements(a: &[Element], b: &[Element], strength: Strength) -> Ordering {
    let order = primaries(a).cmp(primaries(b));
    if order != Ordering::Equal || strength == Strength::Primary {
        return order;
    }
    let order = secondaries(a).cmp(secondaries(b));
    if order != Ordering::Equal || strength == Strength::Secondary {
        return order;
    }
    tertiaries(a).cmp(tertiaries(b))
}

/// Comparison of two ASCII strings without allocating. Every ASCII character
/// stands alone and has the common secondary weight, so equal primaries mean
/// equal secondaries.
fn compare_ascii(a: &[u8], b: &[u8], strength: Strength) -> Ordering {
    let order = ascii_weights(a, |e| e.primary).cmp(ascii_weights(b, |e| e.primary));
    if order != Ordering::Equal || strength != Strength::Tertiary {
        return order;
    }
    let tertiary: fn(Element) -> u32 = |e| u32::from(e.tertiary);
    ascii_weights(a, tertiary).cmp(ascii_weights(b, tertiary))
}

fn ascii_weights(s: &[u8], weight: fn(Element) -> u32) -> impl Iterator<Item = u32> + '_ {
    let root = root();
    s.iter()
        .filter_map(move |&c| root.ascii(c))
        .map(weight)
        .filter(|w| *w != 0)
}
