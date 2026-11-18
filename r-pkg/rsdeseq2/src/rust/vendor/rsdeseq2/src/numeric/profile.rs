//! Runtime selection of the numeric profile.
//!
//! Results that depend on floating-point details can follow one of three
//! profiles:
//!
//! - `Portable` (the default): binary64 reductions and pure-Rust maths (the
//!   `libm` crate and ports of the special functions), so the same input gives
//!   the same bits on every operating system and architecture.
//! - `RLinuxX86`: R on Linux x86_64 with OpenBLAS. Scalar reductions use the
//!   x87 80-bit extended format and the Haswell dot-product schedule, and
//!   `log`, `exp`, `pow` and similar come from the host C library (glibc), as
//!   in R. Available on x86_64 hosts.
//! - `RMacos`: R on macOS. Binary64 reductions with the host C library
//!   (Apple's libm), as in R there. Available on macOS.
//!
//! The profile is a process-wide setting. When `RSDESEQ2_NUMERIC_PROFILE` is
//! set when it is first read (`portable`, `r-linux-x86_64`, `r-macos`, or
//! `r-host` for the R profile of this host), that profile is used; otherwise
//! `Portable`. The R package selects the R profile of its host.

use std::sync::atomic::{AtomicU8, Ordering};

use crate::errors::DeseqError;

/// Which floating-point arithmetic the platform-sensitive kernels follow.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum NumericProfile {
    /// Binary64 reductions and pure-Rust maths: identical bits on every host.
    Portable,
    /// R on Linux x86_64 with OpenBLAS: x87 extended reductions, the Haswell
    /// dot schedule and glibc maths. Available only on x86_64 hosts.
    RLinuxX86,
    /// R on macOS: binary64 reductions and the host C maths library.
    /// Available only on macOS.
    RMacos,
}

impl NumericProfile {
    /// Name used by `RSDESEQ2_NUMERIC_PROFILE` and in manifests.
    pub fn name(self) -> &'static str {
        match self {
            Self::Portable => "portable",
            Self::RLinuxX86 => "r-linux-x86_64",
            Self::RMacos => "r-macos",
        }
    }

    /// Parse a profile name as accepted by `RSDESEQ2_NUMERIC_PROFILE`;
    /// `r-host` names [`NumericProfile::r_host`].
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "portable" => Some(Self::Portable),
            "r-linux-x86_64" => Some(Self::RLinuxX86),
            "r-macos" => Some(Self::RMacos),
            "r-host" => Some(Self::r_host()),
            _ => None,
        }
    }

    /// The profile matching R on this host: `RLinuxX86` on Linux x86_64,
    /// `RMacos` on macOS and `Portable` elsewhere.
    pub fn r_host() -> Self {
        if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
            Self::RLinuxX86
        } else if cfg!(target_os = "macos") {
            Self::RMacos
        } else {
            Self::Portable
        }
    }

    /// The profile used when none is selected: `Portable`.
    pub fn host_default() -> Self {
        Self::Portable
    }

    /// Whether this host can run the profile.
    pub fn is_available(self) -> bool {
        match self {
            Self::Portable => true,
            Self::RLinuxX86 => cfg!(target_arch = "x86_64"),
            Self::RMacos => cfg!(target_os = "macos"),
        }
    }

    /// Whether `log`, `exp`, `pow` and the other elementary functions come
    /// from the host C library (the R profiles) rather than pure Rust.
    pub fn uses_host_libm(self) -> bool {
        !matches!(self, Self::Portable)
    }
}

const UNSET: u8 = u8::MAX;
static PROFILE: AtomicU8 = AtomicU8::new(UNSET);

fn encode(profile: NumericProfile) -> u8 {
    match profile {
        NumericProfile::Portable => 0,
        NumericProfile::RLinuxX86 => 1,
        NumericProfile::RMacos => 2,
    }
}

/// The active profile. The first call reads `RSDESEQ2_NUMERIC_PROFILE`; when it
/// is unset, unknown or unavailable, [`NumericProfile::host_default`]
/// (`Portable`) is used.
pub fn numeric_profile() -> NumericProfile {
    match PROFILE.load(Ordering::Relaxed) {
        0 => NumericProfile::Portable,
        1 => NumericProfile::RLinuxX86,
        2 => NumericProfile::RMacos,
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

/// Whether elementary functions come from the host C library.
#[inline]
pub(crate) fn host_libm() -> bool {
    numeric_profile().uses_host_libm()
}

#[cfg(test)]
mod tests {
    use super::NumericProfile;

    // The process-wide setting is shared by every test, so these tests only
    // exercise the pure helpers and never call `set_numeric_profile`.
    #[test]
    fn names_round_trip_and_host_default_is_available() {
        for profile in [
            NumericProfile::Portable,
            NumericProfile::RLinuxX86,
            NumericProfile::RMacos,
        ] {
            assert_eq!(NumericProfile::from_name(profile.name()), Some(profile));
        }
        assert_eq!(NumericProfile::from_name("linux"), None);
        assert!(NumericProfile::host_default().is_available());
        assert!(NumericProfile::r_host().is_available());
        assert_eq!(
            NumericProfile::from_name("r-host"),
            Some(NumericProfile::r_host())
        );
        assert!(NumericProfile::Portable.is_available());
        assert_eq!(
            NumericProfile::RLinuxX86.is_available(),
            cfg!(target_arch = "x86_64")
        );
    }
}
