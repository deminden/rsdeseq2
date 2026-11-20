// Exact independent R references selected by the active arithmetic profile.
macro_rules! formula_reference {
    ($stem:literal) => {
        match rsdeseq2::numeric_profile() {
            rsdeseq2::NumericProfile::RLinuxX86 => {
                include_str!(concat!("data/", $stem, "_linux_x86_64.json"))
            }
            rsdeseq2::NumericProfile::RLinuxArm64 => {
                include_str!(concat!("data/", $stem, "_linux_aarch64.json"))
            }
            _ => include_str!(concat!("data/", $stem, ".json")),
        }
    };
}

pub(crate) use formula_reference;
