//! Runtime selection of platform-specific R arithmetic.
//!
//! A few formula and shrinkage kernels can imitate one particular R build:
//! R on Linux x86_64, which accumulates scalar reductions in the x87 80-bit
//! extended format and uses OpenBLAS's Haswell dot-product schedule. Other R
//! builds (for example on macOS) round every step to binary64. Results under
//! the two profiles can differ in the last bits of formula transforms such as
//! `scale()` and `poly()`, in NaN payloads of R summaries, and in far-tail
//! posterior moments.
//!
//! The profile is a process-wide setting. Unless `RSDESEQ2_NUMERIC_PROFILE`
//! (`portable` or `r-linux-x86_64`) is set when it is first read, it follows
//! the host, as the compile-time switch it replaces did: `RLinuxX86` on Linux
//! x86_64 and `Portable` elsewhere. Selecting `Portable` on Linux x86_64 gives
//! binary64 reductions there too; transcendental functions still come from the
//! host's libm, so `Portable` output is identical across hosts only for the
//! reductions, not for `log10`, `expm1`, `pow` and similar.

use std::sync::atomic::{AtomicU8, Ordering};

use crate::errors::DeseqError;

/// Which R arithmetic the platform-sensitive kernels reproduce.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum NumericProfile {
    /// Binary64 reductions. Matches R builds without extended-precision
    /// reductions, such as R on macOS (with that host's libm).
    Portable,
    /// R on Linux x86_64 with OpenBLAS: x87 extended reductions and the
    /// Haswell dot schedule. Available only on x86_64 hosts.
    RLinuxX86,
}

impl NumericProfile {
    /// Name used by `RSDESEQ2_NUMERIC_PROFILE` and in manifests.
    pub fn name(self) -> &'static str {
        match self {
            Self::Portable => "portable",
            Self::RLinuxX86 => "r-linux-x86_64",
        }
    }

    /// Parse a profile name as accepted by `RSDESEQ2_NUMERIC_PROFILE`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "portable" => Some(Self::Portable),
            "r-linux-x86_64" => Some(Self::RLinuxX86),
            _ => None,
        }
    }

    /// The profile matching R on this host: `RLinuxX86` on Linux x86_64,
    /// `Portable` elsewhere.
    pub fn host_default() -> Self {
        if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
            Self::RLinuxX86
        } else {
            Self::Portable
        }
    }

    /// Whether this host can run the profile.
    pub fn is_available(self) -> bool {
        match self {
            Self::Portable => true,
            Self::RLinuxX86 => cfg!(target_arch = "x86_64"),
        }
    }
}

const UNSET: u8 = u8::MAX;
static PROFILE: AtomicU8 = AtomicU8::new(UNSET);

fn encode(profile: NumericProfile) -> u8 {
    match profile {
        NumericProfile::Portable => 0,
        NumericProfile::RLinuxX86 => 1,
    }
}

/// The active profile. The first call reads `RSDESEQ2_NUMERIC_PROFILE`; when it
/// is unset, unknown or unavailable, [`NumericProfile::host_default`] is used.
pub fn numeric_profile() -> NumericProfile {
    match PROFILE.load(Ordering::Relaxed) {
        0 => NumericProfile::Portable,
        1 => NumericProfile::RLinuxX86,
        _ => {
            let profile = std::env::var("RSDESEQ2_NUMERIC_PROFILE")
                .ok()
                .and_then(|name| NumericProfile::from_name(&name))
                .filter(|profile| profile.is_available())
                .unwrap_or_else(NumericProfile::host_default);
            // Another thread may have set the profile meanwhile; keep theirs.
            let _ = PROFILE.compare_exchange(
                UNSET,
                encode(profile),
                Ordering::Relaxed,
                Ordering::Relaxed,
            );
            numeric_profile()
        }
    }
}

/// Select the profile for the rest of the process. Call it before any
/// analysis; changing it between analyses is allowed but makes their outputs
/// incomparable.
pub fn set_numeric_profile(profile: NumericProfile) -> Result<(), DeseqError> {
    if !profile.is_available() {
        return Err(DeseqError::UnsupportedFeature {
            feature: format!("numeric profile {} on this host", profile.name()),
        });
    }
    PROFILE.store(encode(profile), Ordering::Relaxed);
    Ok(())
}

/// Whether the x87/OpenBLAS kernels are active.
#[inline]
pub(crate) fn r_linux_x86() -> bool {
    cfg!(target_arch = "x86_64") && numeric_profile() == NumericProfile::RLinuxX86
}

#[cfg(test)]
mod tests {
    use super::NumericProfile;

    // The process-wide setting is shared by every test, so these tests only
    // exercise the pure helpers and never call `set_numeric_profile`.
    #[test]
    fn names_round_trip_and_host_default_is_available() {
        for profile in [NumericProfile::Portable, NumericProfile::RLinuxX86] {
            assert_eq!(NumericProfile::from_name(profile.name()), Some(profile));
        }
        assert_eq!(NumericProfile::from_name("linux"), None);
        assert!(NumericProfile::host_default().is_available());
        assert!(NumericProfile::Portable.is_available());
        assert_eq!(
            NumericProfile::RLinuxX86.is_available(),
            cfg!(target_arch = "x86_64")
        );
    }
}
