//! Check the collator against orders recorded from R.
//!
//! Each file lists strings in R's order, one per line as hexadecimal code
//! points, with the relation to the next string: `<` (primary), `<<`
//! (secondary), `<<<` (tertiary), `=` (tie), or `<` alone when the level was
//! not recorded. The collation is taken from the header written by
//! `scripts/provenance/probe_collation_oracle.R` ("requested xx;"), and is
//! the root collation otherwise. Files written by
//! `scripts/provenance/probe_collation_full.R` are read up to their list of
//! ignorable code points. A directory argument checks every oracle listed in
//! `docs/provenance/observations/collation_r462_macos/oracle_summary.tsv`
//! (root, and for each supported tailoring the probe strings, every code
//! point and the expansions) and their pair counts.
//!
//! ```sh
//! cargo run --release -p rcompat-collate --example verify -- /tmp/rcollate/root_single.tsv
//! cargo run --release -p rcompat-collate --example verify -- /tmp/rcollate
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

/// The oracle files of the summary with their numbers of adjacent pairs.
fn summary() -> Result<Vec<(String, usize)>, String> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../docs/provenance/observations/collation_r462_macos/oracle_summary.tsv"
    );
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let pairs = f.get(2).and_then(|p| p.parse().ok());
            match (f.first(), pairs) {
                (Some(name), Some(pairs)) => Ok(((*name).to_owned(), pairs)),
                _ => Err(format!("{path}: bad line {l}")),
            }
        })
        .collect()
}

fn main() -> ExitCode {
    let mut failed = false;
    let mut paths = Vec::new();
    for arg in std::env::args().skip(1) {
        if std::path::Path::new(&arg).is_dir() {
            match summary() {
                Ok(files) => paths.extend(
                    files
                        .into_iter()
                        .map(|(f, n)| (format!("{arg}/{f}"), Some(n))),
                ),
                Err(e) => {
                    println!("{e}");
                    failed = true;
                }
            }
        } else {
            paths.push((arg, None));
        }
    }
    // Check the files in parallel and report them in order.
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    let mut results = Vec::with_capacity(paths.len());
    for chunk in paths.chunks(threads) {
        std::thread::scope(|scope| {
            let handles: Vec<_> = chunk
                .iter()
                .map(|(path, _)| scope.spawn(|| check(path)))
                .collect();
            for h in handles {
                results.push(
                    h.join()
                        .unwrap_or_else(|_| Err("verification panicked".to_owned())),
                );
            }
        });
    }
    for ((path, expected), result) in paths.into_iter().zip(results) {
        match result {
            Ok((pairs, bad)) => {
                println!("{path}: {pairs} adjacent pairs, {} differ", bad.len());
                if let Some(n) = expected.filter(|&n| n != pairs) {
                    println!("  the summary records {n} adjacent pairs");
                    failed = true;
                }
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
