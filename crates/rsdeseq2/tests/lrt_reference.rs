mod common;

use common::*;
use rsdeseq2::prelude::*;

fn assert_f64_or_missing(actual: f64, expected: Option<f64>, label: &str) {
    assert_reference_value(actual, expected, 1e-5, 1e-5, label);
}

#[path = "lrt_reference/fixed.rs"]
mod fixed;

#[path = "lrt_reference/native.rs"]
mod native;

#[path = "lrt_reference/weighted.rs"]
mod weighted;
