# rsdeseq2

[![License: MIT + BSD notices](https://img.shields.io/badge/license-MIT%20%2B%20BSD-blue.svg)](THIRD_PARTY_NOTICES.md)

`rsdeseq2` runs DESeq2-compatible RNA-seq differential expression analyses in
Rust. Use it from the command line, as a Rust library, or as an R package.
The CLI and Rust library run without R or DESeq2.

## Installation

The steps below build **0.6.0** from this repository. You’ll need
**Rust 1.96.0 or newer** with Cargo. For published versions, see
[crates.io](https://crates.io/crates/rsdeseq2) and
[GitHub Releases](https://github.com/deminden/rsdeseq2/releases).

### Command line

From the repository root, install the `rsdeseq2` command.

```sh
cargo install --path crates/rsdeseq2 --locked
```

### Rust library

From your application, add the library using the path to your local copy.

```sh
cargo add rsdeseq2 --path /path/to/rsdeseq2/crates/rsdeseq2
```

Use `DeseqBuilder::run(&Analysis)` to run a full analysis or stop after a
specific stage. [Rust example](crates/rsdeseq2/README.md#rust-usage) ·
[Analysis API](docs/design/run-entry-point.md)

### R package

With R and its native build tools installed, run this from the repository root.

```sh
R CMD INSTALL r-pkg/rsdeseq2
```

The R package accepts matrices, lists, formulas, and supported `DESeqDataSet`
inputs. [R usage guide](docs/r-wrapper.md)

## Run a Wald test

Use `counts.tsv` for raw, nonnegative integer counts. Put gene IDs in the
first column (`gene`) and sample IDs in the remaining column headers. In
`design.tsv`, put sample IDs in the first column (`sample`), followed by the
numeric design columns.

```sh
rsdeseq2 wald \
  --counts counts.tsv \
  --design design.tsv \
  --fit-type mean \
  --coefficient 1 \
  --output results.tsv
```

This tests the second design coefficient and writes the results to
`results.tsv`. The example uses a constant dispersion trend (`--fit-type mean`).
[CLI commands, file formats, and contrasts](docs/cli.md)

## Supported Workflows

- Estimate size factors with the ratio or `poscounts` method, or supply
  normalization factors.
- Estimate dispersions and run negative-binomial Wald or likelihood-ratio
  tests. Test coefficients or contrasts and compute adjusted p-values.
- Apply Cook’s distance diagnostics, outlier replacement and refitting, and
  independent filtering.
- Run variance-stabilizing (VST) and regularized log (rlog) transformations,
  or shrink log fold changes with the normal method.

R-style formula evaluation and `ashr`- and `apeglm`-style shrinkage are
experimental. The R package supports selected Bioconductor result objects.
`glmGamPoi`, plotting, and custom formula functions are outside the supported
workflows.
[Compatibility guide](docs/compatibility.md) ·
[Supported methods by interface](docs/native-supported-contract.md)

## Validation and Performance

The 0.6.0 acceptance checks at `a0808b2` passed on macOS arm64 and Linux
arm64. Both numeric profiles passed the full 56,200-gene GTEx replay.

In a separate comparison against R 4.6.2 and DESeq2 1.52.0, both profiles
passed all 130 precision checks across eleven GTEx v8 Kidney Cortex Wald
contrasts and agreed with DESeq2 at adjusted p-value cutoffs of 0.05 and
0.1. The runs estimated normalization and dispersions with Cook’s masking
and independent filtering disabled. The table shows the largest absolute
differences from DESeq2.

| Numeric profile | Log fold change | Standard error | p-value | Adjusted p-value |
| --- | ---: | ---: | ---: | ---: |
| Portable | `3.65e-6` | `6.18e-5` | `7.25e-5` | `3.57e-5` |
| Host R | `7.15e-6` | `2.75e-5` | `3.74e-5` | `2.87e-5` |

With four workers on macOS arm64, the portable CLI ran all eleven contrasts
in 15.85 seconds and used 259.41 MiB of peak resident memory in one run.
[Benchmark setup and results](docs/reports/hardening-progress.md#batched-acceptance-after-the-portable-arithmetic-changes)

Both numeric profiles also passed the macOS arm64 integration tests and
matched the reference outputs after the builder API cleanup.
[API cleanup validation](docs/reports/hardening-progress.md#analysis-api-cleanup)

## Documentation

- [0.6.0 release notes](docs/releases/0.6.0.md)
- [Algorithms](docs/algorithms.md)
- [Benchmarks](docs/benchmarks.md) and [reproducibility](docs/reproducibility.md)
- [Development guide](docs/development.md) and [release checklist](docs/release.md)

## Support and Citation

When reporting a bug or a difference from DESeq2, include the rsdeseq2 version,
the command or API call, the input dimensions, and the comparison versions of
R and DESeq2. [GitHub Issues](https://github.com/deminden/rsdeseq2/issues)

Cite the exact rsdeseq2 release you used. For the statistical method, cite
Love, Huber, and Anders (2014).
[Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8)

## License

Project code is MIT-licensed. Adapted numerical components also carry
BSD-3-Clause and BSD-2-Clause-Views terms.
[License](LICENSE) · [Third-party notices](THIRD_PARTY_NOTICES.md)
