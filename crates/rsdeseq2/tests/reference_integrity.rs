//! Regression checks for the fail-closed numerical oracle helpers.
mod common;

#[test]
#[should_panic(expected = "required DESeq2 reference missing")]
fn missing_oracle_is_a_failure() {
    common::reference_file("this_required_reference_does_not_exist.tsv");
}

#[test]
fn numerical_assertions_reject_nonfinite_values_and_both_missingness_directions() {
    for (a, b) in [
        (f64::NAN, 1.0),
        (f64::INFINITY, 1.0),
        (1.0, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
        (1.0, f64::NAN),
    ] {
        assert!(
            std::panic::catch_unwind(|| common::assert_float_close(
                a,
                b,
                1e-6,
                1e-6,
                "invalid pair"
            ))
            .is_err()
        );
    }
    for (a, b) in [
        (None, Some(1.0)),
        (Some(1.0), None),
        (Some(f64::NAN), Some(1.0)),
    ] {
        assert!(
            std::panic::catch_unwind(|| common::assert_option_close(
                a,
                b,
                1e-6,
                1e-6,
                "invalid missingness"
            ))
            .is_err()
        );
    }
    common::assert_option_close(None, None, 1e-6, 1e-6, "matching missingness");
    common::assert_float_close(1.0, 1.0, 1e-6, 1e-6, "matching values");
}
