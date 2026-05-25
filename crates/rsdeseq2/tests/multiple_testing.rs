use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

#[test]
fn bh_adjust_known_example() {
    let adjusted = bh_adjust_f64(&[0.01, 0.04, 0.03, 0.002]).unwrap();
    let expected = [0.02, 0.04, 0.04, 0.008];
    for (actual, expected) in adjusted.iter().zip(expected) {
        assert_relative_eq!(*actual, expected, epsilon = 1e-12);
    }
}

#[test]
fn bh_adjust_preserves_missing_values() {
    let adjusted = bh_adjust(&[Some(0.01), None, Some(0.04), Some(0.03)]);
    assert_eq!(adjusted[1], None);
    assert_relative_eq!(adjusted[0].unwrap(), 0.03, epsilon = 1e-12);
    assert_relative_eq!(adjusted[2].unwrap(), 0.04, epsilon = 1e-12);
    assert_relative_eq!(adjusted[3].unwrap(), 0.04, epsilon = 1e-12);
}

#[test]
fn original_bh_adjust_handles_ties_endpoints_and_omitted_values() {
    let adjusted = bh_adjust(&[
        Some(0.0),
        None,
        Some(0.02),
        Some(0.02),
        Some(1.0),
        Some(0.5),
        Some(f64::NAN),
        Some(1.2),
    ]);

    assert_eq!(adjusted[0], Some(0.0));
    assert_eq!(adjusted[1], None);
    assert_relative_eq!(adjusted[2].unwrap(), 1.0 / 30.0, epsilon = 1e-15);
    assert_relative_eq!(adjusted[3].unwrap(), 1.0 / 30.0, epsilon = 1e-15);
    assert_eq!(adjusted[4], Some(1.0));
    assert_relative_eq!(adjusted[5].unwrap(), 0.625, epsilon = 1e-15);
    assert_eq!(adjusted[6], None);
    assert_eq!(adjusted[7], None);
}

#[test]
fn bh_adjust_rejects_invalid_f64_input() {
    assert!(bh_adjust_f64(&[0.1, f64::NAN]).is_err());
    assert!(bh_adjust_f64(&[0.1, 1.2]).is_err());
}
