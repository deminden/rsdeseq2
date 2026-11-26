//! Write the facts the table derivation takes from Rust rather than from R:
//! canonical and compatibility decompositions (`unicode-normalization`) and
//! case (`char::to_lowercase`, `char::is_uppercase`, `char::is_lowercase`).
//!
//! One line per code point whose decompositions or case differ from the code
//! point itself: code point, NFD, NFKD, lowercase mapping, flags (`l`
//! lowercase, `u` uppercase). Used by
//! `scripts/provenance/build_collation_table.py`.
//!
//! ```sh
//! cargo run -p rcompat-collate --example dump_decompositions > /tmp/rcollate/decompositions.tsv
//! ```

use std::io::{self, BufWriter, Write};
use unicode_normalization::UnicodeNormalization;

fn hex(text: &str) -> String {
    text.chars()
        .map(|c| format!("{:04X}", u32::from(c)))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() -> io::Result<()> {
    let mut out = BufWriter::new(io::stdout().lock());
    for c in (1..=0x10_FFFF).filter_map(char::from_u32) {
        let text = c.to_string();
        let nfd: String = text.nfd().collect();
        let nfkd: String = text.nfkd().collect();
        let lower: String = c.to_lowercase().collect();
        let flags = format!(
            "{}{}",
            if c.is_lowercase() { "l" } else { "" },
            if c.is_uppercase() { "u" } else { "" }
        );
        if nfd == text && nfkd == text && lower == text && flags.is_empty() {
            continue;
        }
        writeln!(
            out,
            "{:04X}\t{}\t{}\t{}\t{}",
            u32::from(c),
            hex(&nfd),
            hex(&nfkd),
            hex(&lower),
            flags
        )?;
    }
    out.flush()
}
