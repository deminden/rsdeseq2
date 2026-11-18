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
