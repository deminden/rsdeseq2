use super::*;
// Result rows share one construction path; only private pipeline stages
// that finalize Cook/filter policy may defer the initial BH calculation.

#[derive(Clone, Copy)]
pub(super) struct RowAssembly<'a> {
    gene_names: Option<&'a [String]>,
    dispersions: Option<&'a [f64]>,
    adjust: bool,
}
impl<'a> RowAssembly<'a> {
    fn adjusted(gene_names: Option<&'a [String]>, dispersions: Option<&'a [f64]>) -> Self {
        Self {
            gene_names,
            dispersions,
            adjust: true,
        }
    }
    fn unadjusted(gene_names: Option<&'a [String]>, dispersions: Option<&'a [f64]>) -> Self {
        Self {
            gene_names,
            dispersions,
            adjust: false,
        }
    }
}

// Wald and LRT rows retain stable public constructor contracts.
#[path = "row_builders/wald.rs"]
pub(super) mod wald;
pub use wald::*;
#[path = "row_builders/lrt.rs"]
pub(super) mod lrt;
pub use lrt::*;
