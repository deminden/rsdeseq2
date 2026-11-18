//! `--help` snapshots for every CLI command, taken before the 4.4 argument
//! regrouping so any change to the user-visible interface shows up in review.
//! Regenerate with `RSDESEQ2_UPDATE_HELP=1` after an intended change.

use std::path::PathBuf;
use std::process::Command;

const COMMANDS: [&str; 8] = [
    "",
    "size-factors",
    "base-mean",
    "normalized-counts",
    "vst",
    "rlog",
    "wald",
    "lrt",
];

fn help(command: &str) -> String {
    let mut process = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"));
    if !command.is_empty() {
        process.arg(command);
    }
    let output = process.arg("--help").output().expect("run rsdeseq2 --help");
    assert!(output.status.success(), "{command} --help failed");
    String::from_utf8(output.stdout).expect("UTF-8 help")
}

#[test]
fn help_text_matches_snapshots() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/cli_help");
    let update = std::env::var_os("RSDESEQ2_UPDATE_HELP").is_some();
    let mut changed = Vec::new();
    for command in COMMANDS {
        let name = if command.is_empty() {
            "rsdeseq2"
        } else {
            command
        };
        let path = dir.join(format!("{name}.txt"));
        let actual = help(command);
        if update {
            std::fs::write(&path, &actual).expect("write help snapshot");
            continue;
        }
        let expected = std::fs::read_to_string(&path).unwrap_or_default();
        if actual != expected {
            changed.push(name);
        }
    }
    assert!(
        changed.is_empty(),
        "help text changed for {changed:?}; review and regenerate"
    );
}
