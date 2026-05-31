//! Private implementation for the stable public module surface.
//!
//! The source is split by workflow so the public module can stay small
//! while preserving existing paths through re-exports.

// CLI imports, argument definitions, output path types, and option conversions.
include!("implementation/args.rs");

// Top-level CLI dispatch and command execution.
include!("implementation/run.rs");

// CLI result adapters and beta-prior analysis builders.
include!("implementation/analysis.rs");

// Sidecar and diagnostic output writers.
include!("implementation/writers.rs");

// Input readers, normalization controls, contrast parsing, and small helpers.
include!("implementation/helpers.rs");

// Unit tests for CLI export helpers.
include!("implementation/tests.rs");
