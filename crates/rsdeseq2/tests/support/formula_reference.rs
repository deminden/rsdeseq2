// Exact independent references for the active numeric profile's R arithmetic:
// the `_linux_x86_64` files come from R on Linux x86_64, the others from R on macOS.
macro_rules! formula_reference {
    ($stem:literal) => {
        if (rsdeseq2::numeric_profile() == rsdeseq2::NumericProfile::RLinuxX86) {
            include_str!(concat!("data/", $stem, "_linux_x86_64.json"))
        } else {
            include_str!(concat!("data/", $stem, ".json"))
        }
    };
}
