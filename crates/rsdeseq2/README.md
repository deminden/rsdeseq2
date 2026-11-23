# rsdeseq2

DESeq2-compatible RNA-seq differential expression in Rust, available as a
library and standalone command-line tool. Normalization, dispersion estimation,
negative-binomial fitting, statistical tests, transformations and supported
shrinkage all run in Rust. **R and DESeq2 are not runtime dependencies.**

Use the CLI with count/design tables, or embed the library to inspect and reuse
individual analysis stages. An optional R package wraps the same numerical core.

## Install

Requires **Rust 1.96.0 or newer**. From the root of the 0.6.0 checkout,
install the CLI with:

```sh
cargo install --path crates/rsdeseq2 --locked
```

To use the checkout as a library, run
`cargo add rsdeseq2 --path /absolute/path/to/rsdeseq2/crates/rsdeseq2`
from your application. [crates.io](https://crates.io/crates/rsdeseq2)
lists published versions. See the [release procedure][release] for how
archives are checked and published.

## Workflows

- **Normalize counts:** ratio/`poscounts` size factors, normalization factors,
  normalized counts and weighted base metadata.
- **Fit and test:** gene-wise/trend/MAP dispersions; negative-binomial Wald and
  likelihood-ratio tests; supported observation weights, offsets and contrasts.
- **Assemble results:** Cook's diagnostics and replacement/refitting,
  independent filtering, adjusted p-values and fitted-state reuse.
- **Transform and shrink:** VST, rlog, and supported normal, `ashr` and `apeglm`
  methods, with shared Rust execution through the R interface.

The implementation is independently validated against original R/DESeq2 outputs.
It does not implement the complete Bioconductor object API. `glmGamPoi`, plotting,
custom formula functions and method-specific unsupported options remain outside
scope. Read [compatibility and limits][compatibility] and the
[supported contract][contract] before selecting a workflow. The public Rust API
is pre-1.0 and may change between minor releases.

## Validation and Performance

The [validation summary][comparison] reports 0.6.0 precision and
performance for eleven full GTEx v8 Kidney Cortex contrasts. Independent
reference tests cover numerical stages and complete workflows. Separate
tests check that standalone and R-wrapped results agree.
[Testing and parity][testing] covers ARM release checks, R, sanitizers
and the full GTEx replay.

## Rust Usage

A minimal normalization example:

```rust
use rsdeseq2::prelude::*;

fn main() -> Result<(), DeseqError> {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 12, 20, 24,
            0,  0,  5,  7,
            100, 80, 90, 120,
        ],
    )?;

    let output = DeseqBuilder::new()
        .size_factor_method(SizeFactorMethod::Ratio)
        .run(&Analysis::without_design(&counts))?;
    let fit = output.fit();

    println!("{:?}", fit.size_factors);
    println!("{:?}", fit.base_mean);
    Ok(())
}
```

Cargo features, all enabled by default:

| Feature | Provides |
| --- | --- |
| `cli` | The `rsdeseq2` binary and `rsdeseq2::cli` (clap, sha2, same-file) |
| `formula` | ICU locale collation for formula factor levels; without it only `"C"` collation is accepted |
| `shrinkage` | `rsdeseq2::shrinkage` (ashr- and apeglm-style posteriors) |

A library-only checkout dependency can use `default-features = false`
with `features = ["formula", "shrinkage"]`.

## CLI

Run a Wald analysis with your count and numeric design tables:

```sh
rsdeseq2 wald \
  --counts counts.tsv \
  --design design.tsv \
  --fit-type mean \
  --coefficient 1 \
  --output results.tsv
```

`--fit-type mean` explicitly selects a constant mean dispersion trend.
Coefficient indices are zero-based. See the [CLI guide][cli] for input schemas,
trend choices, named coefficients, contrasts, LRT and transformation commands.

## Support and License

Report reproducible problems through [GitHub Issues](https://github.com/deminden/rsdeseq2/issues).
Include the version, inputs or a minimal reproducer, and reference-tool versions
for numerical differences. [Development instructions][development] describe
local checks and independent reference generation.

Project code is MIT-licensed; adapted numerical components additionally carry
BSD-3-Clause and BSD-2-Clause-Views terms. The crate includes `LICENSE` and
`THIRD_PARTY_NOTICES.md`; see [license][license] and [third-party notices][notices].

[contract]: https://github.com/deminden/rsdeseq2/blob/main/docs/native-supported-contract.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/main/docs/compatibility.md
[comparison]: https://github.com/deminden/rsdeseq2/blob/main/README.md#validation-and-performance
[release]: https://github.com/deminden/rsdeseq2/blob/main/docs/release.md
[testing]: https://github.com/deminden/rsdeseq2/blob/main/docs/testing-parity.md
[cli]: https://github.com/deminden/rsdeseq2/blob/main/docs/cli.md
[development]: https://github.com/deminden/rsdeseq2/blob/main/docs/development.md
[license]: https://github.com/deminden/rsdeseq2/blob/main/LICENSE
[notices]: https://github.com/deminden/rsdeseq2/blob/main/THIRD_PARTY_NOTICES.md
