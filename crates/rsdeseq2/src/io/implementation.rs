//! Private implementation behind the stable public module API.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Input row-label types and TSV readers.
#[path = "implementation/read.rs"]
mod read;
pub use read::*;

// Label alignment helpers for sample and gene keyed inputs.
#[path = "implementation/align.rs"]
mod align;
pub use align::*;

// TSV field parsers and primitive matrix writers.
#[path = "implementation/write_matrices.rs"]
mod write_matrices;
pub use write_matrices::*;

// DESeq result, diagnostics, and independent-filtering TSV writers.
#[path = "implementation/write_results.rs"]
mod write_results;
pub use write_results::*;

// Formatting and export validation helpers.
#[path = "implementation/format_validation.rs"]
mod format_validation;
use format_validation::*;

// Unit tests for I/O helpers.
#[cfg(test)]
#[path = "implementation/tests.rs"]
mod tests;
