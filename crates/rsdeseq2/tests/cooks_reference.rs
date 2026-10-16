mod common;

use common::*;
use rsdeseq2::prelude::*;

fn read_count_matrix(name: &str) -> CountMatrix {
    rsdeseq2::io::read_count_matrix_tsv(reference_file(name)).unwrap()
}

fn read_design_matrix(name: &str) -> DesignMatrix {
    rsdeseq2::io::read_design_matrix_tsv(reference_file(name)).unwrap()
}

fn read_replacement_options() -> CooksReplacementOptions {
    let rows = read_reference_tsv("cooks_replacement_options.tsv");
    let row = &rows[0];
    CooksReplacementOptions {
        trim: parse_required_f64(row, "trim"),
        cooks_cutoff: parse_required_f64(row, "cooksCutoff"),
        min_replicates: parse_required_f64(row, "minReplicates") as usize,
        which_samples: Some(read_replacement_size_factor_rows().1),
    }
}

fn read_replacement_size_factor_rows() -> (Vec<f64>, Vec<bool>) {
    let rows = read_reference_tsv("cooks_replacement_size_factors.tsv");
    (
        rows.iter()
            .map(|row| parse_required_f64(row, "size_factor"))
            .collect(),
        rows.iter()
            .map(|row| parse_required_bool(row, "replaceable"))
            .collect(),
    )
}

#[test]
fn cooks_replacement_counts_match_deseq2_reference() {
    let counts = read_count_matrix("cooks_replacement_counts.tsv");
    let expected_candidates = read_count_matrix("cooks_replacement_candidate_counts.tsv");
    let expected_replaced = read_count_matrix("cooks_replacement_replaced_counts.tsv");
    let cooks = read_reference_matrix("cooks_replacement_cooks.tsv");
    let design = read_design_matrix("cooks_replacement_design.tsv");
    let (size_factors, replaceable) = read_replacement_size_factor_rows();
    let options = read_replacement_options();

    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &size_factors,
        None,
        &cooks,
        &design,
        &options,
    )
    .unwrap();

    assert_eq!(
        output.candidate_replacement_counts.as_slice(),
        expected_candidates.as_slice()
    );
    assert_eq!(
        output.replaced_counts.as_slice(),
        expected_replaced.as_slice()
    );
    assert_eq!(output.replaceable_samples, replaceable);
}

#[test]
fn cooks_replacement_refit_plan_matches_deseq2_reference() {
    let counts = read_count_matrix("cooks_replacement_counts.tsv");
    let rows = read_reference_tsv("cooks_replacement_rows.tsv");
    let cooks = read_reference_matrix("cooks_replacement_cooks.tsv");
    let design = read_design_matrix("cooks_replacement_design.tsv");
    let (size_factors, _replaceable) = read_replacement_size_factor_rows();
    let options = read_replacement_options();

    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let plan = prepare_cooks_replacement_refit(
        &counts,
        &normalized,
        &size_factors,
        None,
        &cooks,
        &design,
        &options,
    )
    .unwrap();

    assert_eq!(plan.replacement.replace.len(), rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(
            plan.replacement.replace[gene],
            Some(parse_required_bool(row, "replace"))
        );
        assert_eq!(
            plan.replaced_all_zero[gene],
            parse_required_bool(row, "allZero")
        );
        assert_float_close(
            plan.replaced_base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("replacement baseMean gene {gene}"),
        );
        assert_float_close(
            plan.replaced_base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("replacement baseVar gene {gene}"),
        );
        assert_eq!(
            plan.new_all_zero_rows.contains(&gene),
            parse_required_bool(row, "newAllZero")
        );
        assert_eq!(
            plan.refit_rows.contains(&gene),
            parse_required_bool(row, "refitReplace")
        );
        assert_option_close(
            plan.post_refit_max_cooks[gene],
            parse_optional_f64(row, "postRefitMaxCooks"),
            1e-12,
            1e-12,
            &format!("post-refit maxCooks gene {gene}"),
        );
    }
}
