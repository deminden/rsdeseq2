//! Shared independent-posterior transports and comparison gates.
use std::error::Error;
#[allow(dead_code)]
#[path = "../../examples/evaluate_r_interface/formula.rs"]
pub(crate) mod formula;
#[allow(dead_code)]
#[path = "../../examples/evaluate_r_interface/heavy_tail.rs"]
pub(crate) mod heavy_tail;
#[allow(dead_code)]
#[path = "../../examples/evaluate_r_interface/likelihood.rs"]
pub(crate) mod likelihood;
#[allow(dead_code)]
#[path = "../../examples/evaluate_r_interface/random.rs"]
pub(crate) mod random;
#[allow(dead_code)]
#[path = "../../examples/evaluate_r_interface/regression.rs"]
pub(crate) mod regression;
pub(crate) fn bits(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}
pub(crate) fn parse_numeric_input(s: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(s) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(s, 16)?)
    } else if s == "NA" {
        f64::from_bits(0x7ff0_0000_0000_07a2)
    } else {
        s.parse()?
    })
}
pub(crate) fn parse_numeric_inputs(v: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    v.iter().map(|s| parse_numeric_input(s)).collect()
}
pub(crate) fn compare(
    a: &[f64],
    b: &serde_json::Value,
    absolute: f64,
    relative: f64,
    context: &str,
    failures: &mut Vec<String>,
) {
    let b = b.as_array().unwrap();
    assert_eq!(a.len(), b.len(), "{context}");
    for (i, (&a, b)) in a.iter().zip(b).enumerate() {
        let b = parse_numeric_input(b.as_str().unwrap()).unwrap();
        if a == b
            || (a.is_nan()
                && b.is_nan()
                && (a.to_bits() as u32 == 1954) == (b.to_bits() as u32 == 1954))
        {
            continue;
        }
        if a.is_finite() && b.is_finite() && (a - b).abs() <= absolute + relative * b.abs() {
            continue;
        }
        failures.push(format!(
            "{context}[{i}]: {a:.17e} != {b:.17e}; error {:.3e}",
            (a - b).abs()
        ));
    }
}
