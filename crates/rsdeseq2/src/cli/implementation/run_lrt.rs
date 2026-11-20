use super::*;

// Execute the LRT command after output staging.
pub(super) fn run_lrt(args: LrtArgs) -> Result<(), DeseqError> {
    let LrtArgs {
        counts,
        design,
        reduced_design,
        normalization,
        target,
        contrast,
        filtering,
        outputs,
    } = args;
    let counts = read_count_matrix_tsv(counts)?;
    let design = read_cli_design_matrix(design, &counts)?;
    let reduced_design = read_cli_design_matrix(reduced_design, &counts)?;
    let builder = DeseqBuilder::new()
        .automatic_dispersion_means()
        .size_factor_method(normalization.method.into())
        .fit_type(normalization.fit_type.into());
    let builder = configure_cli_results(builder, &counts, normalization, filtering)?;
    validate_cli_result_target("LRT", &target, &contrast)?;
    let cutoff = resolve_cooks_cutoff(
        builder.current_cooks_cutoff(),
        design.n_samples(),
        design.n_coefficients(),
    )?;
    let analysis = run_cli_model(
        &builder,
        &counts,
        &design,
        Some(&reduced_design),
        target,
        contrast,
        cutoff,
    )?;
    write_cli_analysis_outputs(outputs, &counts, &analysis)
}
