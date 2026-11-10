// Result rows share one construction path; only private pipeline stages
// that finalize Cook/filter policy may defer the initial BH calculation.

#[derive(Clone, Copy)]
struct RowAssembly<'a> {
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
include!("row_builders/wald.rs");
include!("row_builders/lrt.rs");
