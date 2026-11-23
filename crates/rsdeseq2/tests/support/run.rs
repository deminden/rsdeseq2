//! Extract the expected variant of an analysis that has already been run.
use rsdeseq2::prelude::*;

#[track_caller]
pub(crate) fn run_fit(output: Result<RunOutput, DeseqError>) -> Result<DeseqFit, DeseqError> {
    match output? {
        RunOutput::Fit(fit) | RunOutput::Results { fit, .. } => Ok(fit),
        _ => panic!("expected a fit without replacement"),
    }
}

#[track_caller]
pub(crate) fn run_results(
    output: Result<RunOutput, DeseqError>,
) -> Result<(DeseqFit, DeseqResults), DeseqError> {
    match output? {
        RunOutput::Results { fit, results } => Ok((fit, results)),
        _ => panic!("expected result rows without replacement"),
    }
}

#[track_caller]
pub(crate) fn run_replacement(
    output: Result<RunOutput, DeseqError>,
) -> Result<CooksReplacementTestOutput, DeseqError> {
    match output? {
        RunOutput::Replacement(output) => Ok(output),
        _ => panic!("expected Cook's replacement"),
    }
}

#[track_caller]
pub(crate) fn run_wald_replacement(
    output: Result<RunOutput, DeseqError>,
) -> Result<CooksReplacementWaldOutput, DeseqError> {
    match run_replacement(output)? {
        CooksReplacementTestOutput::Wald(output) => Ok(output),
        _ => panic!("expected Wald replacement"),
    }
}

#[track_caller]
pub(crate) fn run_lrt_replacement(
    output: Result<RunOutput, DeseqError>,
) -> Result<CooksReplacementLrtOutput, DeseqError> {
    match run_replacement(output)? {
        CooksReplacementTestOutput::Lrt(output) => Ok(output),
        _ => panic!("expected LRT replacement"),
    }
}
