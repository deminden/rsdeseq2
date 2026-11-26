//! Order properties of the collators on random strings.

use proptest::prelude::*;
use proptest::sample::select;
use rcompat_collate::Collator;
use std::cmp::Ordering;

/// Characters from the probed ranges (letters, case and accent variants,
/// combining marks, ignorables, punctuation, contractions, expansions) and
/// a few outside them (fallback order).
fn alphabet() -> Vec<char> {
    let mut chars: Vec<char> = "aAbBeEilLoOsSzZ019 _-.,'\u{e1}\u{c1}\u{e9}\u{df}\u{e6}\u{c6}\u{f0}\u{142}\u{17f}\
        \u{b7}\u{301}\u{300}\u{308}\u{327}\u{ad}\u{200b}\u{1}\u{3b1}\u{3ac}\u{391}\u{438}\u{439}\u{418}\u{306}\
        \u{2153}\u{fb01}\u{ff41}\u{1c4}\u{1c5}\u{1c6}\u{5d0}\u{627}\u{653}\u{3042}\u{ac00}\u{4e00}\u{1f600}\u{10ffff}"
        .chars()
        .collect();
    chars.sort_unstable();
    chars.dedup();
    chars
}

fn text() -> impl Strategy<Value = String> {
    proptest::collection::vec(select(alphabet()), 0..7).prop_map(|v| v.into_iter().collect())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    #[test]
    fn root_order_is_antisymmetric_and_matches_sort_keys(a in text(), b in text()) {
        let root = Collator::root();
        let ab = root.compare(&a, &b);
        prop_assert_eq!(ab, root.compare(&b, &a).reverse());
        prop_assert_eq!(ab, root.sort_key(&a).cmp(&root.sort_key(&b)));
        prop_assert_eq!(root.compare(&a, &a), Ordering::Equal);
    }

    #[test]
    fn root_order_is_transitive(a in text(), b in text(), c in text()) {
        let root = Collator::root();
        let mut v = [a, b, c];
        v.sort_by(|x, y| root.compare(x, y));
        prop_assert_ne!(root.compare(&v[0], &v[1]), Ordering::Greater);
        prop_assert_ne!(root.compare(&v[1], &v[2]), Ordering::Greater);
        prop_assert_ne!(root.compare(&v[0], &v[2]), Ordering::Greater);
    }

    #[test]
    fn c_order_is_byte_order(a in text(), b in text()) {
        let c = Collator::c();
        prop_assert_eq!(c.compare(&a, &b), a.as_bytes().cmp(b.as_bytes()));
        prop_assert_eq!(c.compare(&a, &b), c.sort_key(&a).cmp(&c.sort_key(&b)));
    }
}

#[test]
fn c_and_root_orders_differ_as_in_r() {
    let c = Collator::new("C").unwrap();
    let root = Collator::new("root").unwrap();
    // Upper case before lower case in C, letters before case in root.
    assert_eq!(c.compare("B", "a"), Ordering::Less);
    assert_eq!(root.compare("B", "a"), Ordering::Greater);
    assert_eq!(root.compare("a", "A"), Ordering::Less);
    // Accents after every unaccented letter in C, as a tie-breaker in root.
    assert_eq!(c.compare("\u{e9}", "f"), Ordering::Greater);
    assert_eq!(root.compare("\u{e9}", "f"), Ordering::Less);
    // Punctuation and digits before letters in root, by code point in C.
    assert_eq!(c.compare("a", "_"), Ordering::Greater);
    assert_eq!(c.compare("_", "1"), Ordering::Greater);
    assert_eq!(root.compare("_", "1"), Ordering::Less);
    // Precomposed and combining accents tie in root, not in C.
    assert_eq!(root.compare("e\u{301}", "\u{e9}"), Ordering::Equal);
    assert_ne!(c.compare("e\u{301}", "\u{e9}"), Ordering::Equal);
    // Ignorable characters do not count in root.
    assert_eq!(root.compare("a\u{ad}b", "ab"), Ordering::Equal);
}

#[test]
fn locale_names() {
    assert!(Collator::new("en_US").unwrap() == Collator::root());
    assert!(Collator::new("ca").unwrap() == Collator::root());
    assert!(Collator::new("POSIX").unwrap().is_bytewise());
    assert_eq!(Collator::new("sv_SE").unwrap().name(), "sv");
    assert_eq!(Collator::new("nb-NO").unwrap().name(), "no");
    assert_eq!(Collator::new("nn").unwrap(), Collator::new("no").unwrap());
    assert_eq!(
        Collator::new("de_DE@collation=phonebook").unwrap().name(),
        "de@collation=phonebook"
    );
    assert_eq!(Collator::new("da_DK").unwrap().name(), "da");
    assert_eq!(Collator::new("is").unwrap().name(), "is");
    assert_eq!(Collator::new("el-GR").unwrap().name(), "el");
    assert_eq!(Collator::new("fr-CA").unwrap().name(), "fr_CA");
    assert!(Collator::new("fr_FR").unwrap() == Collator::root());
    for name in rcompat_collate::SUPPORTED_LOCALES {
        assert!(Collator::new(name).is_ok(), "{name}");
    }
    for name in [
        "cy",
        "eo",
        "vi",
        "es@collation=traditional",
        "en_US_POSIX",
        "sr_Latn",
        "zh",
        "ja",
        "c",
        "english",
    ] {
        let error = Collator::new(name).unwrap_err();
        assert!(error.to_string().contains(name), "{name}");
    }
}

fn sorted<'a>(locale: &str, words: &[&'a str]) -> Vec<&'a str> {
    let collator = Collator::new(locale).unwrap();
    let mut out = words.to_vec();
    out.sort_by(|a, b| collator.compare(a, b));
    out
}

#[test]
fn tailorings_as_in_r() {
    let w = ["\u{f6}rn", "zebra", "\u{e5}sa", "anna", "\u{e4}ng"];
    assert_eq!(
        sorted("sv", &w),
        ["anna", "zebra", "\u{e5}sa", "\u{e4}ng", "\u{f6}rn"]
    );
    assert_eq!(
        sorted("root", &w),
        ["\u{e4}ng", "anna", "\u{e5}sa", "\u{f6}rn", "zebra"]
    );
    assert_eq!(
        sorted("cs", &["chata", "hrad", "cesta", "ia"]),
        ["cesta", "hrad", "chata", "ia"]
    );
    assert_eq!(
        sorted("es", &["\u{f1}u", "nz", "o"]),
        ["nz", "\u{f1}u", "o"]
    );
    assert_eq!(
        sorted("tr", &["ia", "\u{131}z", "h"]),
        ["h", "\u{131}z", "ia"]
    );
    assert_eq!(sorted("ru", &["b", "\u{431}"]), ["\u{431}", "b"]);
    assert_eq!(sorted("root", &["b", "\u{431}"]), ["b", "\u{431}"]);
    assert_eq!(
        sorted("de@collation=phonebook", &["af", "\u{e4}", "ae"]),
        ["ae", "\u{e4}", "af"]
    );
    assert_eq!(sorted("hu", &["csa", "cz", "da"]), ["cz", "csa", "da"]);
}

#[test]
fn whole_unicode_range_is_ordered() {
    // Ideographs in radical-and-stroke order, unassigned code points last.
    let root = Collator::root();
    assert_eq!(root.compare("\u{4e00}", "\u{20000}"), Ordering::Less);
    assert_eq!(root.compare("\u{9fa5}", "\u{20000}"), Ordering::Greater);
    assert_eq!(root.compare("\u{ff5e}", "\u{4e00}"), Ordering::Less);
    assert_eq!(root.compare("\u{4e00}", "\u{e000}"), Ordering::Less);
    assert_eq!(root.compare("\u{10ffff}", "\u{fffd}"), Ordering::Less);
    // A Hangul syllable ties with its jamo.
    assert_eq!(
        root.compare("\u{ac00}", "\u{1100}\u{1161}"),
        Ordering::Equal
    );
}
