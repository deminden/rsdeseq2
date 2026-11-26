//! Orders observed in R 4.6.2 (ICU root collation), written by
//! `scripts/provenance/probe_collation_orders.R`. Each line holds the code
//! points of one string and its relation to the next string in R's order.

use rcompat_collate::Collator;
use std::cmp::Ordering;

fn parse(data: &str) -> Vec<(String, String)> {
    data.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|line| {
            let (points, relation) = line.split_once('\t').unwrap();
            let text = points
                .split(' ')
                .map(|p| char::from_u32(u32::from_str_radix(p, 16).unwrap()).unwrap())
                .collect();
            (text, relation.to_owned())
        })
        .collect()
}

fn mismatches(data: &str) -> Vec<String> {
    let root = Collator::root();
    let rows = parse(data);
    let mut out = Vec::new();
    for pair in rows.windows(2) {
        let ((a, relation), (b, _)) = (&pair[0], &pair[1]);
        let expected = if relation == "=" {
            Ordering::Equal
        } else {
            Ordering::Less
        };
        let actual = root.compare(a, b);
        let by_key = root.sort_key(a).cmp(&root.sort_key(b));
        if actual != expected || by_key != actual {
            out.push(format!(
                "{a:?} vs {b:?}: R {expected:?}, compare {actual:?}, key {by_key:?}"
            ));
        }
    }
    out
}

fn check(name: &str, data: &str) {
    let bad = mismatches(data);
    assert!(
        bad.is_empty(),
        "{name}: {} adjacent pairs differ from R, first: {:#?}",
        bad.len(),
        &bad[..bad.len().min(20)]
    );
}

#[test]
fn ascii_strings_follow_r() {
    check("ascii", include_str!("data/sorted_ascii.tsv"));
}

#[test]
fn accented_words_greek_and_cyrillic_follow_r() {
    check("words", include_str!("data/sorted_words.tsv"));
}

#[test]
fn random_strings_over_the_probed_characters_follow_r() {
    check("inventory", include_str!("data/sorted_inventory.tsv"));
}

#[test]
fn near_duplicates_and_combining_marks_follow_r() {
    check("variants", include_str!("data/sorted_variants.tsv"));
}

#[test]
fn other_scripts_and_ideographs_follow_r() {
    check(
        "other scripts",
        include_str!("data/sorted_other_scripts.tsv"),
    );
}

/// Samples of R's order under each supported tailoring
/// (`scripts/provenance/sample_collation_orders.py`).
#[test]
fn tailorings_follow_r() {
    let text = include_str!("data/sorted_tailorings.tsv");
    let mut sections = 0;
    for section in text.split("\n@ ").skip(1) {
        let (locale, body) = section.split_once('\n').unwrap();
        let collator = Collator::new(locale).unwrap();
        let rows = parse(body);
        for pair in rows.windows(2) {
            let ((a, relation), (b, _)) = (&pair[0], &pair[1]);
            let expected = if relation == "=" {
                Ordering::Equal
            } else {
                Ordering::Less
            };
            assert_eq!(collator.compare(a, b), expected, "{locale}: {a:?} vs {b:?}");
            assert_eq!(
                collator.sort_key(a).cmp(&collator.sort_key(b)),
                expected,
                "{locale}"
            );
        }
        sections += 1;
    }
    assert_eq!(sections, 26);
}
