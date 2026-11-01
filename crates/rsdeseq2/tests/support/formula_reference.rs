// Exact independent references for the platform's R arithmetic.
macro_rules! formula_reference {
    ($stem:literal) => {
        if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
            include_str!(concat!("data/", $stem, "_linux_x86_64.json"))
        } else {
            include_str!(concat!("data/", $stem, ".json"))
        }
    };
}
