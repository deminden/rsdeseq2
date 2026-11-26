//! Check the collator against orders recorded from R.
//!
//! Each file lists strings in R's order, one per line as hexadecimal code
//! points, with the relation to the next string: `<` (primary), `<<`
//! (secondary), `<<<` (tertiary), `=` (tie), or `<` alone when the level was
//! not recorded. The collation is taken from the header written by
//! `scripts/provenance/probe_collation_oracle.R` ("requested xx;"), and is
//! the root collation otherwise. Files written by
//! `scripts/provenance/probe_collation_full.R` are read up to their list of
//! ignorable code points.
//!
//! ```sh
//! cargo run --release -p rcompat-collate --example verify -- /tmp/rcollate/root_single.tsv
//! ```

use std::cmp::Ordering;
use std::process::ExitCode;

use rcompat_collate::{Collator, Strength};

fn decode(hex: &str) -> Option<String> {
    hex.split(' ')
        .map(|p| u32::from_str_radix(p, 16).ok().and_then(char::from_u32))
        .collect()
}

fn relation(c: &Collator, a: &str, b: &str) -> &'static str {
    let order = c.compare(a, b);
    let sign = match order {
        Ordering::Less => "<",
        Ordering::Equal => return "=",
        Ordering::Greater => ">",
    };
    let level = if c.compare_with_strength(a, b, Strength::Primary) != Ordering::Equal {
        1
    } else if c.compare_with_strength(a, b, Strength::Secondary) != Ordering::Equal {
        2
    } else {
        3
    };
    match (sign, level) {
        ("<", 1) => "<",
        ("<", 2) => "<<",
        ("<", _) => "<<<",
        (_, 1) => ">",
        (_, 2) => ">>",
        _ => ">>>",
    }
}

fn check(path: &str) -> Result<(usize, Vec<String>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let locale = text
        .lines()
        .take_while(|l| l.starts_with('#'))
        .find_map(|l| {
            l.split("requested ")
                .nth(1)
                .and_then(|r| r.split(';').next())
        })
        .unwrap_or("root")
        .to_owned();
    let collator = Collator::new(&locale).map_err(|e| e.to_string())?;
    let levels = text.contains("\t<<");
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with("# completely ignorable") {
            break;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let hex = fields.next().unwrap_or("");
        let rel = fields.next().unwrap_or("");
        let s = decode(hex).ok_or_else(|| format!("{path}: bad line {line}"))?;
        rows.push((hex.to_owned(), s, rel.to_owned()));
    }
    let mut bad = Vec::new();
    for pair in rows.windows(2) {
        let (ha, a, want) = &pair[0];
        let (hb, b, _) = &pair[1];
        let mut got = relation(&collator, a, b);
        if !levels && got.starts_with('<') {
            got = "<";
        }
        let by_key = collator.sort_key(a).cmp(&collator.sort_key(b));
        if got != want || by_key != collator.compare(a, b) {
            bad.push(format!("[{ha}] [{hb}] R {want} got {got} (key {by_key:?})"));
        }
    }
    Ok((rows.len().saturating_sub(1), bad))
}

fn main() -> ExitCode {
    let mut failed = false;
    for path in std::env::args().skip(1) {
        match check(&path) {
            Ok((pairs, bad)) => {
                println!("{path}: {pairs} adjacent pairs, {} differ", bad.len());
                for line in bad.iter().take(25) {
                    println!("  {line}");
                }
                failed |= !bad.is_empty();
            }
            Err(e) => {
                println!("{e}");
                failed = true;
            }
        }
    }
    if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}
