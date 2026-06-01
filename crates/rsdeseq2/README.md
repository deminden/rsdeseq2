# rsdeseq2 Core Crate

This crate contains the Rust numerical core for `rsdeseq2`.

Project source and status docs live on GitHub:
[repository][repo], [gap analysis][gap-analysis],
[compatibility notes][compatibility], and [benchmarks][benchmarks].

The current API implements count matrix storage, design matrix validation,
size-factor estimation and caller-supplied size factors, normalized counts,
base means, BH adjusted p-values, fixed-dispersion NB GLM fitting for the
initial unweighted paths, selected-coefficient Wald results with optional
t-distribution p-values and LFC-threshold alternatives, supplied-dispersion
Wald/LRT pipelines, limited native-dispersion Wald/LRT branches, and an
inspectable fit-state skeleton. The pipeline code
skips all-zero rows for GLM fitting, expands them back as missing outputs, and records Cook's distances
plus `maxCooks`. Result rows can apply DESeq2-style Cook's cutoff p-value
masking, including the explicit two-group low-count Cook's heuristic helper,
and base-mean independent filtering before BH adjustment, including an R
`stats::lowess`-shaped smoother for the default DESeq2 threshold grid.
`DeseqResults::column_names()` exposes the current core and diagnostic result
column contract. A
primitive Cook's outlier replacement-count transform and replacement-refit
planning helper are available, along with a limited replacement-refit path for
the GLM-mu native Wald and LRT branches. The current
dispersion foundation includes linear-mu starts,
Cox-Reid objective scoring, the parametric trend form
`asymptDisp + extraPois / mean`, deterministic dispersion prior variance, and
the first no-weight MAP dispersion stage. A limited native Wald pipeline wires
that linear-mu/parametric/MAP subset into GLM fitting and result assembly.
Full DESeq2 dispersion estimation and GLM parity are intentionally left as
explicit future work.

```rust
use rsdeseq2::prelude::*;

let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16])?;
let fit = DeseqBuilder::new().fit_size_factors_and_base_means(&counts)?;
# Ok::<(), DeseqError>(())
```

[repo]: https://github.com/deminden/rsdeseq2
[gap-analysis]: https://github.com/deminden/rsdeseq2/blob/main/docs/deseq2-gap-analysis.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/main/docs/compatibility.md
[benchmarks]: https://github.com/deminden/rsdeseq2/blob/main/docs/benchmarks.md
