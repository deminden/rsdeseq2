use super::*;

#[test]
fn cli_wald_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("wald-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_ambiguous_coefficient_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--coefficient-name",
            "(Intercept)",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient '(Intercept)' resolves ambiguously after R-style cleanup",
    );
}

#[test]
fn cli_wald_rejects_factor_level_contrast_without_sample_levels() {
    let dir = temp_dir("wald-contrast-factor-missing-levels");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "condition",
            "--contrast-numerator",
            "B",
            "--contrast-denominator",
            "A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: factor-level results contrast requires --contrast-sample-levels",
    );
}

#[test]
fn cli_wald_rejects_ambiguous_cleaned_factor_level_alias() {
    let dir = temp_dir("wald-contrast-factor-ambiguous-level");
    let design = dir.join("design.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("wald.tsv");
    write_cell_type_contrast_design_fixture(&design);
    write_ambiguous_cell_type_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "cell.type",
            "--contrast-numerator",
            "B.cell",
            "--contrast-denominator",
            "A.cell",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "denominator level 'A.cell' resolves ambiguously",
    );
}

#[test]
fn cli_lrt_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("lrt-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let output = dir.join("lrt.tsv");
    write_ambiguous_coefficient_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reduced_design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--coefficient-name",
            "(Intercept)",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient '(Intercept)' resolves ambiguously after R-style cleanup",
    );
}

#[test]
fn cli_lrt_rejects_factor_level_contrast_without_sample_levels() {
    let dir = temp_dir("lrt-contrast-factor-missing-levels");
    let design = dir.join("design.tsv");
    let output = dir.join("lrt.tsv");
    write_standard_contrast_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reference_data_path("design_reduced.tsv").to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "condition",
            "--contrast-numerator",
            "B",
            "--contrast-denominator",
            "A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: factor-level results contrast requires --contrast-sample-levels",
    );
}

#[test]
fn cli_lrt_rejects_ambiguous_cleaned_factor_level_alias() {
    let dir = temp_dir("lrt-contrast-factor-ambiguous-level");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_cell_type_contrast_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);
    write_ambiguous_cell_type_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reduced_design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "cell.type",
            "--contrast-numerator",
            "B.cell",
            "--contrast-denominator",
            "A.cell",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "denominator level 'A.cell' resolves ambiguously",
    );
}

#[test]
fn cli_lrt_rejects_sample_levels_without_factor_level_contrast() {
    let dir = temp_dir("lrt-sample-levels-without-factor");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            reference_data_path("design_full.tsv").to_str().unwrap(),
            "--reduced-design",
            reference_data_path("design_reduced.tsv").to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for factor-level contrast inputs",
    );
}
