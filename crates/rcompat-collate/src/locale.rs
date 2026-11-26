//! Collation names accepted by [`crate::Collator::new`].
//!
//! The root languages are those for which `icuGetCollate("actual")` reported
//! `root` after `icuSetCollate(locale = ...)` in R 4.6.2 with ICU 78.3
//! (`docs/provenance/observations/collation_r462_macos/locales.tsv`).

use crate::tailoring::Tailoring;

/// What a collation name selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Collation {
    Bytes,
    Root,
    Tailored(Tailoring),
}

/// Languages whose collation ICU reported as the root collation.
const ROOT_LANGUAGES: [&str; 12] = [
    "ca", "de", "en", "eu", "fr", "ga", "id", "it", "ms", "nl", "pt", "sw",
];

/// Locale names with a tailoring of their own that [`crate::Collator::new`]
/// accepts (a region may be added, as in `"sv_SE"`; `nb` and `nn` are the
/// tailoring ICU calls `no`).
pub const SUPPORTED_LOCALES: [&str; 27] = [
    "az",
    "be",
    "bg",
    "bs",
    "cs",
    "de@collation=phonebook",
    "es",
    "et",
    "fi",
    "gl",
    "hr",
    "hu",
    "lt",
    "lv",
    "mk",
    "nb",
    "nn",
    "no",
    "pl",
    "ro",
    "ru",
    "sk",
    "sl",
    "sr",
    "sv",
    "tr",
    "uk",
];

pub(crate) fn classify(name: &str) -> Option<Collation> {
    if name == "C" || name == "POSIX" {
        return Some(Collation::Bytes);
    }
    let lower = name.to_ascii_lowercase();
    if matches!(lower.as_str(), "" | "root" | "und") {
        return Some(Collation::Root);
    }
    let (tag, keyword) = match lower.split_once('@') {
        Some((tag, keyword)) => (tag, Some(keyword)),
        None => (lower.as_str(), None),
    };
    let mut subtags = tag.split(['-', '_']);
    let language = subtags.next()?;
    let region = subtags.next();
    if subtags.next().is_some() || region.is_some_and(|r| !is_region(r)) {
        return None;
    }
    match (language, region, keyword) {
        ("de", _, Some("collation=phonebook")) => {
            Tailoring::find("de@collation=phonebook").map(Collation::Tailored)
        }
        (_, _, Some(_)) => None,
        ("fr", Some("ca"), None) => None,
        ("nb" | "nn" | "no", _, None) => Tailoring::find("no").map(Collation::Tailored),
        (language, _, None) if ROOT_LANGUAGES.contains(&language) => Some(Collation::Root),
        (language, _, None) => Tailoring::find(language).map(Collation::Tailored),
    }
}

fn is_region(subtag: &str) -> bool {
    (subtag.len() == 2 && subtag.bytes().all(|b| b.is_ascii_alphabetic()))
        || (subtag.len() == 3 && subtag.bytes().all(|b| b.is_ascii_digit()))
}
