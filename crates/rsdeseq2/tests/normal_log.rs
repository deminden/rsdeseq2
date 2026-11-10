use rsdeseq2::math::{normal_log_cdf, normal_log_upper_tail};
use serde::Deserialize;
#[derive(Deserialize)]
struct Reference {
    arguments: Vec<String>,
    log_cdf: Vec<String>,
    log_upper: Vec<String>,
}
fn number(s: &str) -> f64 {
    match s {
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        "NaN" => f64::NAN,
        _ => s.parse().unwrap(),
    }
}
fn within_ulps(a: f64, b: f64, context: &str) {
    if b.is_nan() {
        assert!(a.is_nan(), "{context}");
        return;
    }
    if a == b {
        return;
    }
    assert!(
        a.is_finite() && b.is_finite() && a.to_bits().abs_diff(b.to_bits()) <= 16,
        "{context}: {a:.17e} != {b:.17e}, distance {}",
        a.to_bits().abs_diff(b.to_bits())
    );
}
#[test]
fn standalone_log_normal_tails_match_original_r_across_transition_and_underflow() {
    let r: Reference = serde_json::from_str(include_str!("data/normal_log_r462.json")).unwrap();
    assert!(r.arguments.len() > 650);
    for (i, x) in r.arguments.iter().enumerate() {
        let x = number(x);
        within_ulps(
            normal_log_cdf(x),
            number(&r.log_cdf[i]),
            &format!("cdf({x})"),
        );
        within_ulps(
            normal_log_upper_tail(x),
            number(&r.log_upper[i]),
            &format!("upper({x})"),
        );
    }
    assert!(normal_log_cdf(-40.).is_finite());
    assert!(normal_log_cdf(-1.4e154).is_finite());
}
