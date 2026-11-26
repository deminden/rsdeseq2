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

/// Other scripts follow R except where a CJK ideograph decides the order:
/// ideographs are outside the probed ranges and take the code point
/// fallback, while R orders them differently.
#[test]
fn other_scripts_follow_r_apart_from_ideographs() {
    let root = Collator::root();
    let rows = parse(include_str!("data/sorted_other_scripts.tsv"));
    let ideograph = |s: &str| s.chars().any(|c| ('\u{4e00}'..='\u{9fff}').contains(&c));
    let mut differ = 0;
    for pair in rows.windows(2) {
        let ((a, relation), (b, _)) = (&pair[0], &pair[1]);
        let expected = if relation == "=" {
            Ordering::Equal
        } else {
            Ordering::Less
        };
        if root.compare(a, b) != expected {
            differ += 1;
            assert!(ideograph(a) && ideograph(b), "{a:?} vs {b:?}");
        }
    }
    assert!(differ < rows.len() / 100, "{differ} pairs differ");
}
