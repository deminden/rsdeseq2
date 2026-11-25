//! Compiled prefit and refined posterior of the two gtex_compact genes whose
//! refined intercepts are most sensitive to last-bit prefit differences.
use crate::support as common;
use common::hex_bits as bits;
use std::error::Error;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/heavy_tail.rs"]
mod heavy_tail;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/random.rs"]
mod random;
fn parse_numeric_input(s: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(s) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(s, 16)?)
    } else {
        s.parse()?
    })
}
fn parse_numeric_inputs(v: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    v.iter().map(|s| parse_numeric_input(s)).collect()
}
fn expected(c: &serde_json::Value, field: &str) -> Vec<f64> {
    c[field]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap())
        .collect()
}
/// Native row-major output in the column-major order of the reference.
fn column_major(
    m: &rsdeseq2::matrix::RowMajorMatrix<f64>,
    rows: usize,
    columns: usize,
) -> Vec<f64> {
    (0..columns)
        .flat_map(|j| (0..rows).map(move |i| *m.get(i, j).unwrap()))
        .collect()
}

#[test]
fn sensitive_gtex_prefits_and_refined_posteriors_match_original_apeglm() {
    let _ = bits;
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/gtex_apeglm_prefit_r462.json")).unwrap();
    let c = &reference["cases"][0];
    let mut input: heavy_tail::Case = serde_json::from_value(c.clone()).unwrap();
    let (g, p) = (input.counts.rows, input.design.columns);
    let mut failures = vec![];

    // Compiled prefit (apeglm method nbinomC), with the public posterior
    // bound of 0.01 log2 units; the observed differences are below 2e-6.
    let bound = 0.01 * std::f64::consts::LN_2;
    input.map_only = true;
    let prefit = heavy_tail::run(&input).unwrap();
    for (i, (a, b)) in column_major(&prefit.mode, g, p)
        .into_iter()
        .zip(expected(c, "prefit_mode"))
        .enumerate()
    {
        if (a - b).abs() > bound {
            failures.push(format!("prefit[{i}]: {a:.17e} != {b:.17e}"));
        }
    }

    // Refined posterior (nbinomCR). Slopes, standard deviations and the
    // reported coefficient's sign probabilities use the public posterior
    // bound of 0.01 log2 units (0.02 for probabilities). The refined
    // intercept is where the original's numerical refinement turns prefit
    // differences of about 1e-12 into differences of about 1e-2: it is
    // checked against a tenth of its posterior standard deviation and the
    // observed difference is recorded in the report.
    input.map_only = false;
    let refined = heavy_tail::run(&input).unwrap();
    let mode = column_major(&refined.mode, g, p);
    let sd = column_major(refined.standard_deviation.as_ref().unwrap(), g, p);
    let reference_sd = expected(c, "sd");
    for (i, (a, b)) in mode.iter().zip(expected(c, "mode")).enumerate() {
        let intercept = i < g;
        let tolerance = if intercept {
            0.1 * reference_sd[i]
        } else {
            bound
        };
        if (a - b).abs() > tolerance {
            failures.push(format!("mode[{i}]: {a:.17e} != {b:.17e}"));
        }
    }
    for (i, (a, b)) in sd.iter().zip(&reference_sd).enumerate() {
        if (a - b).abs() > bound {
            failures.push(format!("sd[{i}]: {a:.17e} != {b:.17e}"));
        }
    }
    for (field, values) in [
        ("fsr", refined.false_sign_rate.as_ref().unwrap()),
        ("sign_svalues", refined.sign_svalues.as_ref().unwrap()),
    ] {
        for (i, (a, b)) in values.iter().zip(expected(c, field)).enumerate() {
            let a = a.unwrap();
            if (a - b).abs() > 0.02 {
                failures.push(format!("{field}[{i}]: {a:.17e} != {b:.17e}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
