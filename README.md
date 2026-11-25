# rsdeseq2

[![Project license MIT](https://img.shields.io/badge/project%20license-MIT-blue.svg)](LICENSE)

`rsdeseq2` runs DESeq2-compatible RNA-seq differential expression analyses in
Rust. Use it from the command line, as a Rust library, or as an R package.
The CLI and Rust library run without R or DESeq2.

In a single-worker benchmark on an Apple M5, the rsdeseq2 0.6.0 CLI ran **4× as
fast as DESeq2** and used **86% less peak RAM** for an end-to-end
differential expression analysis of 56,200 genes across 61 GTEx samples.
[Benchmark details](#single-worker-benchmark)

## Installation

The steps below build **0.7.0** from this repository. You’ll need
**Rust 1.96.0 or newer** with Cargo. For published versions, see
[crates.io](https://crates.io/crates/rsdeseq2) and
[GitHub Releases](https://github.com/deminden/rsdeseq2/releases).

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

### Command line

From the repository root, install the `rsdeseq2` command.

```sh
cargo install --path crates/rsdeseq2 --locked
```

## Compare treated and control samples

In `counts.tsv`, put gene IDs in the first column (`gene`) and raw counts
in one column per sample. In `samples.tsv`, use two columns, `sample` and
`condition`. Give each sample a matching ID and label it `control` or
`treated`.

Both examples write `differential_expression.tsv` with per-gene log2 fold
changes and adjusted p-values. Positive log2 fold changes mean higher
expression in the treated group.

### Rust

Save this as `src/main.rs` in your Rust project. Run `cargo run` from the
project directory with both TSV files there.

```rust
use rsdeseq2::io::{
    align_sample_levels_to_samples, read_count_matrix_tsv, read_sample_levels_tsv,
    write_deseq_results_tsv,
};
use rsdeseq2::prelude::{Analysis, DeseqBuilder, DesignMatrix, ResultColumn};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let counts = read_count_matrix_tsv("counts.tsv")?;
    let samples = read_sample_levels_tsv("samples.tsv")?;
    let sample_names = counts.sample_names().ok_or("missing sample names")?;
    let conditions = align_sample_levels_to_samples(&samples, sample_names)?;

    // Each row is an intercept followed by the treated-group indicator.
    let mut values = Vec::with_capacity(2 * conditions.len());
    for condition in conditions {
        let treated = match condition.as_str() {
            "control" => 0.0,
            "treated" => 1.0,
            _ => return Err(format!("unknown condition: {condition}").into()),
        };
        values.extend([1.0, treated]);
    }
    let design = DesignMatrix::from_row_major(
        counts.n_samples(),
        2,
        values,
        Some(vec!["Intercept".into(), "treated_vs_control".into()]),
    )?;

    let analysis = Analysis::new(&counts, &design)
        .result(ResultColumn::Name("treated_vs_control"));
    let output = DeseqBuilder::new().run(&analysis)?;
    let results = output.results().ok_or("analysis returned no results")?;
    write_deseq_results_tsv("differential_expression.tsv", results)?;
    Ok(())
}
```

### R

The R package runs the same comparison.

```r
counts <- as.matrix(read.delim("counts.tsv", row.names = 1, check.names = FALSE))
samples <- read.delim("samples.tsv", row.names = 1)
stopifnot(all(colnames(counts) %in% rownames(samples)))
samples <- samples[colnames(counts), , drop = FALSE]
samples$condition <- factor(samples$condition, levels = c("control", "treated"))

design <- model.matrix(~ condition, data = samples)
fit <- rsdeseq2::DESeq(counts, design = design)
res <- as.data.frame(rsdeseq2::results(fit, name = "conditiontreated"))

write.table(
  data.frame(gene = rownames(res), res),
  "differential_expression.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
```

For the same workflow from the command line, see the [CLI guide](docs/cli.md).

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

### Single-worker benchmark

On Apple M5 with macOS 27.2, the portable 0.6.0 CLI (Cargo `local` profile)
took 5.99 seconds and used 239 MiB of peak resident memory. DESeq2 1.52.0
on R 4.6.2 took 23.86 seconds and used 1,709 MiB. Each tool ran once with
one worker. The timings include startup, loading the inputs, fitting, and
writing the results.

Both tools fit the same blocked Wald contrast with parametric dispersions
and with Cook’s masking, outlier replacement, and independent filtering
disabled. They agreed on significant genes at adjusted p-value cutoffs of
0.05 and 0.1.
[GTEx fixture](crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md) ·
[Benchmark procedure](docs/benchmarks.md#full-quantification)

## Documentation

- [0.7.0 release notes](docs/releases/0.7.0.md)
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

The project’s source code is [MIT-licensed](LICENSE). It is an independent
implementation: the [provenance ledger](docs/provenance/ledger.md) records,
file by file, the papers, documentation and black-box R runs each part was
written from. The [third-party notices](THIRD_PARTY_NOTICES.md) name the
published coefficients and public-domain tables it uses. With default
features every Rust dependency is available under MIT, the Unlicense or Zlib
([dependency report](docs/provenance/dependencies.md)); the optional
`formula` feature and the R package add the Unicode-3.0 ICU crates.
External validation datasets retain their own terms.
