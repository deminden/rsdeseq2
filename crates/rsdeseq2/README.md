# rsdeseq2

DESeq2-compatible RNA-seq differential expression in Rust, available as a
library and standalone command-line tool. Normalization, dispersion estimation,
negative-binomial fitting, statistical tests, transformations and supported
shrinkage all run in Rust. **R and DESeq2 are not runtime dependencies.**

Use the CLI with count/design tables, or embed the library to inspect and reuse
individual analysis stages. An optional R package wraps the same numerical core.

## Install

Requires **Rust 1.97.1 or newer**.

```sh
# Embed the library.
cargo add rsdeseq2@0.5.0

# Install the CLI.
cargo install rsdeseq2 --version 0.5.0 --locked
```

For source installation from a repository checkout:

```sh
cargo install --path crates/rsdeseq2 --locked
```

[Release procedure][release] describes checked artifacts, version tags and
publication verification.

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

The [main comparison][comparison] reports precision, elapsed-time ratios and
peak memory against DESeq2 for eleven full GTEx v8 Kidney Cortex contrasts.
The [measurement report][gtex-benchmark] provides per-contrast values, executable
and environment identities, and reproduction instructions. These are measured
results for the documented workload, not a guarantee for every design.

Independent reference tests cover numerical stages and complete workflows;
separate tests require identical standalone and R-wrapped results. See
[testing and parity][testing] for platform, R, sanitizer and full-GTEx checks.

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

    let fit = DeseqBuilder::new()
        .size_factor_method(SizeFactorMethod::Ratio)
        .execution_mode(ExecutionMode::Strict)
        .fit_size_factors_and_base_means(&counts)?;

    println!("{:?}", fit.size_factors);
    println!("{:?}", fit.base_mean);
    Ok(())
}
```

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
[gtex-benchmark]: https://github.com/deminden/rsdeseq2/blob/main/docs/reports/gtex-native-wald-benchmark.md
[comparison]: https://github.com/deminden/rsdeseq2/blob/main/README.md#validation-and-performance
[release]: https://github.com/deminden/rsdeseq2/blob/main/docs/release.md
[testing]: https://github.com/deminden/rsdeseq2/blob/main/docs/testing-parity.md
[cli]: https://github.com/deminden/rsdeseq2/blob/main/docs/cli.md
[development]: https://github.com/deminden/rsdeseq2/blob/main/docs/development.md
[license]: https://github.com/deminden/rsdeseq2/blob/main/LICENSE
[notices]: https://github.com/deminden/rsdeseq2/blob/main/THIRD_PARTY_NOTICES.md
