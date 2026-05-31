//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// Input row-label types and TSV readers.
include!("implementation/read.rs");

// Label alignment helpers for sample and gene keyed inputs.
include!("implementation/align.rs");

// TSV field parsers and primitive matrix writers.
include!("implementation/write_matrices.rs");

// DESeq result, diagnostics, and independent-filtering TSV writers.
include!("implementation/write_results.rs");

// Formatting and export validation helpers.
include!("implementation/format_validation.rs");

// Unit tests for I/O helpers.
include!("implementation/tests.rs");
