use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::errors::DeseqError;
use crate::io::{read_count_matrix_tsv, write_base_mean_tsv, write_size_factors_tsv};
use crate::normalization::{base_mean, estimate_size_factors, normalized_counts};
use crate::options::SizeFactorMethod;

/// Command-line arguments for the minimal `rsdeseq2` CLI.
#[derive(Debug, Parser)]
#[command(name = "rsdeseq2")]
#[command(about = "Early DESeq2-compatible Rust workflow stages")]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Estimate sample size factors.
    SizeFactors {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Estimate size factors, normalized counts, and base means.
    BaseMean {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SizeFactorMethodArg {
    Ratio,
    Poscounts,
}

impl From<SizeFactorMethodArg> for SizeFactorMethod {
    fn from(value: SizeFactorMethodArg) -> Self {
        match value {
            SizeFactorMethodArg::Ratio => Self::Ratio,
            SizeFactorMethodArg::Poscounts => Self::PosCounts,
        }
    }
}

/// Parse process arguments and run the CLI.
pub fn run_cli() -> Result<(), DeseqError> {
    run(Cli::parse())
}

fn run(cli: Cli) -> Result<(), DeseqError> {
    match cli.command {
        Commands::SizeFactors {
            counts,
            method,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let size_factors = estimate_size_factors(&counts, method.into())?;
            write_size_factors_tsv(output, counts.sample_names(), &size_factors)
        }
        Commands::BaseMean {
            counts,
            method,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let size_factors = estimate_size_factors(&counts, method.into())?;
            let normalized = normalized_counts(&counts, &size_factors)?;
            let base_mean = base_mean(&normalized)?;
            write_base_mean_tsv(output, counts.gene_names(), &base_mean)
        }
    }
}
