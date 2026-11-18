//! Tests of the `rsdeseq2` binary. They are written against the core crate's
//! fixtures and helpers (`crates/rsdeseq2/tests`), which they include by path.
#![allow(clippy::duplicate_mod, clippy::module_inception)]
#![allow(deprecated)]

#[path = "../../rsdeseq2/tests/cli.rs"]
mod cli;
#[path = "../../rsdeseq2/tests/cli_help.rs"]
mod cli_help;
