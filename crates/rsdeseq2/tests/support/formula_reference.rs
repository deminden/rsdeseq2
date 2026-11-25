// Exact independent R references. They were recorded with R 4.6.2 on macOS,
// so they apply to the macOS R profile and to the portable profile; other
// hosts have no recorded reference and the tests that use one return early.
macro_rules! formula_reference {
    ($stem:literal) => {
        match rsdeseq2::numeric_profile() {
            rsdeseq2::NumericProfile::RMacos | rsdeseq2::NumericProfile::Portable => {
                Some(include_str!(concat!("data/", $stem, ".json")))
            }
            _ => {
                eprintln!("no {} reference for this numeric profile", $stem);
                None
            }
        }
    };
}

pub(crate) use formula_reference;
