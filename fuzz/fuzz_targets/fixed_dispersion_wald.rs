#![no_main]

// The body is shared with the stable smoke test (tests/fuzz_smoke.rs).
#[allow(dead_code)]
#[path = "../../crates/rsdeseq2/tests/fuzz/bodies.rs"]
mod bodies;

libfuzzer_sys::fuzz_target!(|data: &[u8]| bodies::fixed_dispersion_wald(data));
