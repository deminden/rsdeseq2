//! Collation names accepted by [`crate::Collator::new`].
//!
//! The root languages are those for which `icuGetCollate("actual")` reported
//! `root` after `icuSetCollate(locale = ...)` in R 4.6.2 with ICU 78.3
//! (`docs/provenance/observations/collation_r462_macos/locales.tsv`).

/// What a collation name selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Collation {
    Bytes,
    Root,
}

/// Languages whose collation ICU reported as the root collation.
const ROOT_LANGUAGES: [&str; 12] = [
    "ca", "de", "en", "eu", "fr", "ga", "id", "it", "ms", "nl", "pt", "sw",
];

/// Language-region pairs of those languages that ICU reported as tailored.
const TAILORED: [(&str, &str); 1] = [("fr", "ca")];

pub(crate) fn classify(name: &str) -> Option<Collation> {
    if name == "C" || name == "POSIX" {
        return Some(Collation::Bytes);
    }
    let lower = name.to_ascii_lowercase();
    if matches!(lower.as_str(), "" | "root" | "und") {
        return Some(Collation::Root);
    }
    let mut subtags = lower.split(['-', '_']);
    let language = subtags.next()?;
    let region = subtags.next();
    if subtags.next().is_some() || !ROOT_LANGUAGES.contains(&language) {
        return None;
    }
    match region {
        None => Some(Collation::Root),
        Some(region) if is_region(region) && !TAILORED.contains(&(language, region)) => {
            Some(Collation::Root)
        }
        Some(_) => None,
    }
}

fn is_region(subtag: &str) -> bool {
    (subtag.len() == 2 && subtag.bytes().all(|b| b.is_ascii_alphabetic()))
        || (subtag.len() == 3 && subtag.bytes().all(|b| b.is_ascii_digit()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_from_the_locale_probe() {
        for name in [
            "root", "", "und", "en", "en_US", "en-GB", "EN_us", "de_AT", "fr_BE", "pt_BR",
        ] {
            assert_eq!(classify(name), Some(Collation::Root), "{name}");
        }
        for name in ["C", "POSIX"] {
            assert_eq!(classify(name), Some(Collation::Bytes), "{name}");
        }
        for name in [
            "fr_CA",
            "en_US_POSIX",
            "sv",
            "da",
            "es",
            "ru",
            "tr",
            "zh",
            "c",
            "de@collation=phonebook",
            "en-u-co-phonebk",
            "english",
        ] {
            assert_eq!(classify(name), None, "{name}");
        }
    }
}
