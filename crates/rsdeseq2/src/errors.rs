//! The error type shared by every rsdeseq2 stage: `rcompat_numeric`'s
//! error, under its rsdeseq2 name.

pub use rcompat_numeric::errors::Error as DeseqError;
pub use rcompat_numeric::errors::{ErrorKind, NegativeBinomialStage, Requirement};
#[doc(hidden)]
pub use rcompat_numeric::errors::{
    invalid_dimensions, invalid_dispersion, invalid_formula, invalid_options, invalid_parameter,
    invariant, missing_prerequisite, scoped_invalid_options,
};
