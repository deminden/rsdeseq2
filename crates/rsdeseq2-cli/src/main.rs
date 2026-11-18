//! The `rsdeseq2` command-line tool.
//!
//! In 0.6.0 the implementation is still `rsdeseq2::cli` (deprecated there); it
//! moves into this crate in 0.7.0.
#![allow(deprecated)]

use std::process::ExitCode;

fn main() -> ExitCode {
    match rsdeseq2::cli::run_cli() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}
