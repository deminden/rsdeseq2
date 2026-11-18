//! Design rank decisions against R 4.6.2 `qr(x)$rank` and DESeq2 1.52.0
//! `checkFullRank` (reference: tests/data/rank_reference_r462.json, generated
//! by scripts/references/generate_rank_reference.R).

use rsdeseq2::prelude::*;
use serde_json::Value;

/// Parse C99 `%a` hexadecimal floats as R's `sprintf("%a")` writes them.
fn hex_float(text: &str) -> f64 {
    let (sign, rest) = match text.strip_prefix('-') {
        Some(rest) => (-1.0, rest),
        None => (1.0, text),
    };
    let rest = rest.strip_prefix("0x").expect("hex float prefix");
    let (mantissa, exponent) = rest.split_once('p').expect("hex float exponent");
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let mut value = u64::from_str_radix(whole, 16).expect("hex digits") as f64;
    let mut scale = 1.0 / 16.0;
    for digit in fraction.chars() {
        value += digit.to_digit(16).expect("hex digit") as f64 * scale;
        scale /= 16.0;
    }
    sign * value * 2f64.powi(exponent.parse::<i32>().expect("exponent"))
}

#[test]
fn design_rank_matches_r_qr_and_deseq2_full_rank_check() {
    let reference: Value =
        serde_json::from_str(include_str!("data/rank_reference_r462.json")).unwrap();
    assert_eq!(reference["tolerance"].as_f64(), Some(1e-7));
    let cases = reference["cases"].as_array().unwrap();
    assert!(cases.len() >= 20);
    let mut mismatches = Vec::new();
    for case in cases {
        let name = case["name"].as_str().unwrap();
        let rows = case["rows"].as_u64().unwrap() as usize;
        let cols = case["cols"].as_u64().unwrap() as usize;
        let column_major = case["values"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| hex_float(v.as_str().unwrap()))
            .collect::<Vec<_>>();
        let row_major = (0..rows * cols)
            .map(|i| column_major[(i % cols) * rows + i / cols])
            .collect();
        let design = DesignMatrix::from_row_major(rows, cols, row_major, None).unwrap();
        let rank = design.rank().unwrap();
        let rejected = design.validate_full_rank("reference").is_err();
        if rank as u64 != case["rank"].as_u64().unwrap()
            || rejected != case["deseq2_rejects"].as_bool().unwrap()
        {
            mismatches.push(format!(
                "{name}: rank {rank} (R {}), rejected {rejected} (DESeq2 {})",
                case["rank"], case["deseq2_rejects"]
            ));
        }
    }
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}

#[test]
fn hex_floats_parse_exactly() {
    assert_eq!(hex_float("0x1.8p+1"), 3.0);
    assert_eq!(hex_float("-0x1p-2"), -0.25);
    assert_eq!(hex_float("0x0p+0"), 0.0);
    assert_eq!(hex_float("0x1.999999999999ap-4"), 0.1);
}
