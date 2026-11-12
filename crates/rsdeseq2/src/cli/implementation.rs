//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// CLI imports, argument definitions, output path types, and option conversions.
include!("implementation/args.rs");
include!("implementation/paths.rs");

// Top-level CLI dispatch and command execution.
#[path = "implementation/run.rs"]
mod run;
pub use run::run_cli;
mod run_lrt;
mod run_wald;

// CLI result adapters and beta-prior analysis builders.
#[path = "implementation/analysis.rs"]
mod analysis;
use analysis::*;

// Sidecar and diagnostic output writers.
#[path = "implementation/writers.rs"]
mod writers;
use writers::*;

// Input readers, normalization controls, contrast parsing, and small helpers.
#[path = "implementation/helpers.rs"]
mod helpers;
use helpers::*;

// Unit tests for CLI export helpers.
include!("implementation/tests.rs");
