# rsdeseq2

Run DESeq2-compatible RNA-seq differential expression in Rust. Use the crate
in your application or run the standalone command-line tool. Both run without
R or DESeq2.

In a single-worker benchmark on an Apple M5, the 0.6.0 CLI ran **4× as fast as
DESeq2** and used **86% less peak RAM** for an end-to-end analysis of 56,200
genes across 61 GTEx samples. [Benchmark details][benchmark]

## Install

The 0.7.0 release requires **Rust 1.96.0 or newer**.

```sh
# Add the library to your Rust project.
cargo add rsdeseq2@0.7.0

# Install the standalone command.
cargo install rsdeseq2 --version 0.7.0 --locked
```

Building from a checkout? Follow the [source installation instructions][source-install].

## Rust Usage

Compare treated samples with controls using raw gene counts and a numeric
design matrix.

- In `counts.tsv`, put gene IDs in the first column (`gene`) and raw counts
  in one column per sample.
- In `design.tsv`, put one row per sample, with columns `sample`,
  `Intercept` (always `1`), and `treated` (`0` for controls, `1` for treated
  samples).

Save this as `src/main.rs` and run `cargo run` with both files in the project
directory.

```rust
use rsdeseq2::io::{
    align_design_matrix_to_samples, read_count_matrix_tsv, read_labeled_design_matrix_tsv,
    write_deseq_results_tsv,
};
use rsdeseq2::prelude::{Analysis, DeseqBuilder, ResultColumn};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let counts = read_count_matrix_tsv("counts.tsv")?;
    let design = read_labeled_design_matrix_tsv("design.tsv")?;
    let sample_names = counts.sample_names().ok_or("missing sample names")?;
    let design = align_design_matrix_to_samples(design, sample_names)?;

    let analysis = Analysis::new(&counts, &design).result(ResultColumn::Name("treated"));
    let output = DeseqBuilder::new().run(&analysis)?;
    let results = output.results().ok_or("analysis returned no results")?;
    write_deseq_results_tsv("differential_expression.tsv", results)?;
    Ok(())
}
```

The output contains per-gene log2 fold changes and adjusted p-values. A positive
log2 fold change means higher expression in the treated group.
[Analysis API][analysis] · [API reference][api]

## What You Can Do

- Estimate size factors, supply normalization factors, and export normalized counts.
- Estimate dispersions and run Wald or likelihood-ratio tests with coefficients,
  contrasts, offsets, and observation weights.
- Apply Cook’s outlier diagnostics, replacement and refitting, and independent
  filtering. Reuse fitted models for additional result queries.
- Run VST and rlog transformations or shrink log fold changes. R-style formulas
  and `ashr`- and `apeglm`-style shrinkage are experimental.

The [CLI guide][cli] covers file formats and commands. The [R package guide][r]
covers running the same Rust engine from R. Check supported options and
Bioconductor object coverage in the [compatibility guide][compatibility].

## Cargo Features

`cli` and `shrinkage` are enabled by default. Formula factor levels always
sort in R's collation (`C` or the Unicode root collation) through
`rcompat-collate`.

| Feature | Adds |
| --- | --- |
| `cli` | The standalone binary and `rsdeseq2::cli` |
| `shrinkage` | `rsdeseq2::shrinkage` with `ashr`- and `apeglm`-style methods |

To embed the library with no CLI dependencies, disable default features and
enable only the features you need.

## Support and License

When reporting a bug or a numerical difference, include the rsdeseq2 version,
a minimal reproducer, and the reference tool version. [GitHub Issues][issues]

The rsdeseq2 Rust code is [MIT-licensed][license]. The CLI and R distributions
include [third-party notices][notices] for their dependencies.

[analysis]: https://github.com/deminden/rsdeseq2/blob/main/docs/design/run-entry-point.md
[api]: https://docs.rs/rsdeseq2
[benchmark]: https://github.com/deminden/rsdeseq2/blob/main/README.md#single-worker-benchmark
[cli]: https://github.com/deminden/rsdeseq2/blob/main/docs/cli.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/main/docs/compatibility.md
[issues]: https://github.com/deminden/rsdeseq2/issues
[license]: https://github.com/deminden/rsdeseq2/blob/main/LICENSE
[notices]: https://github.com/deminden/rsdeseq2/blob/main/THIRD_PARTY_NOTICES.md
[r]: https://github.com/deminden/rsdeseq2/blob/main/docs/r-wrapper.md
[source-install]: https://github.com/deminden/rsdeseq2/blob/main/README.md#installation
