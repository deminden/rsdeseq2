use std::ffi::c_char;
use std::ptr;

const DIAGNOSTIC_SCHEMA_NAMES: [&[u8]; 7] = [
    b"betaConv\0",
    b"fullBetaConv\0",
    b"reducedBetaConv\0",
    b"betaIter\0",
    b"reducedBetaIter\0",
    b"deviance\0",
    b"maxCooks\0",
];

#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_placeholder() -> i32 {
    0
}

/// Number of DESeq2-style diagnostic schema names exported by the R bridge.
#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_diagnostic_schema_len() -> usize {
    DIAGNOSTIC_SCHEMA_NAMES.len()
}

/// Return one null-terminated diagnostic schema name by index.
///
/// The returned pointer has static lifetime. A null pointer indicates an
/// out-of-range index.
#[unsafe(no_mangle)]
pub extern "C" fn rsdeseq2_r_diagnostic_schema_name(index: usize) -> *const c_char {
    DIAGNOSTIC_SCHEMA_NAMES
        .get(index)
        .map_or(ptr::null(), |name| name.as_ptr().cast())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn diagnostic_schema_contract_is_stable() {
        let names = (0..rsdeseq2_r_diagnostic_schema_len())
            .map(|index| {
                let pointer = rsdeseq2_r_diagnostic_schema_name(index);
                assert!(!pointer.is_null());
                unsafe { CStr::from_ptr(pointer) }
                    .to_str()
                    .unwrap()
                    .to_string()
            })
            .collect::<Vec<_>>();

        assert_eq!(
            names,
            vec![
                "betaConv",
                "fullBetaConv",
                "reducedBetaConv",
                "betaIter",
                "reducedBetaIter",
                "deviance",
                "maxCooks"
            ]
        );
        assert!(rsdeseq2_r_diagnostic_schema_name(names.len()).is_null());
    }
}
