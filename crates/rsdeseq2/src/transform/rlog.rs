use crate::errors::DeseqError;

/// Explicit placeholder for DESeq2's regularized-log transformation.
///
/// The Rust core does not yet implement rlog's dispersion-prior and
/// sample-effect shrinkage workflow. Keep this as a public unsupported marker
/// so callers can distinguish missing rlog support from accidental absence in
/// the transform namespace.
pub fn rlog() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "regularized-log transformation".to_string(),
    })
}
