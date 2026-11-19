#![no_main]

libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    rsdeseq2::cli::fuzz_output_transaction(data).expect("temporary filesystem fixture");
});
