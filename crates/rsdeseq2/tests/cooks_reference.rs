mod common;

use common::*;
use rsdeseq2::prelude::*;

fn read_count_matrix(name: &str) -> Option<CountMatrix> {
    let rows = read_optional_tsv(name)?;
    let samples = reference_sample_names();
    let genes = rows
        .iter()
        .map(|row| {
            row.get("gene")
                .unwrap_or_else(|| panic!("missing gene column in {name}"))
                .to_string()
        })
        .collect::<Vec<_>>();
    let values = rows
        .iter()
        .flat_map(|row| {
            samples.iter().map(|sample| {
                let value = parse_required_f64(row, sample);
                assert!(value >= 0.0 && value.fract() == 0.0);
                value as u32
            })
        })
        .collect::<Vec<_>>();
    Some(
        CountMatrix::from_row_major_u32_with_names(
            rows.len(),
            samples.len(),
            values,
            Some(genes),
            Some(samples),
        )
        .unwrap(),
    )
}

fn read_design_matrix(name: &str) -> Option<DesignMatrix> {
    let rows = read_optional_tsv(name)?;
    let coefficient_names = rows[0]
        .keys()
        .filter(|name| name.as_str() != "sample")
        .cloned()
        .collect::<Vec<_>>();
    let values = rows
        .iter()
        .flat_map(|row| {
            coefficient_names
                .iter()
                .map(|coefficient| parse_required_f64(row, coefficient))
        })
        .collect::<Vec<_>>();
    Some(
        DesignMatrix::from_row_major(
            rows.len(),
            coefficient_names.len(),
            values,
            Some(coefficient_names),
        )
        .unwrap(),
    )
}

fn read_replacement_options() -> Option<CooksReplacementOptions> {
    let rows = read_optional_tsv("cooks_replacement_options.tsv")?;
    let row = &rows[0];
    Some(CooksReplacementOptions {
        trim: parse_required_f64(row, "trim"),
        cooks_cutoff: parse_required_f64(row, "cooksCutoff"),
        min_replicates: parse_required_f64(row, "minReplicates") as usize,
        which_samples: Some(read_replacement_size_factor_rows()?.1),
    })
}

fn read_replacement_size_factor_rows() -> Option<(Vec<f64>, Vec<bool>)> {
    let rows = read_optional_tsv("cooks_replacement_size_factors.tsv")?;
    Some((
        rows.iter()
            .map(|row| parse_required_f64(row, "size_factor"))
            .collect(),
        rows.iter()
            .map(|row| parse_required_bool(row, "replaceable"))
            .collect(),
    ))
}

#[test]
fn cooks_replacement_counts_match_optional_deseq2_reference() {
    let Some(counts) = read_count_matrix("cooks_replacement_counts.tsv") else {
        return;
    };
    let Some(expected_candidates) = read_count_matrix("cooks_replacement_candidate_counts.tsv")
    else {
        return;
    };
    let Some(expected_replaced) = read_count_matrix("cooks_replacement_replaced_counts.tsv") else {
        return;
    };
    let Some(cooks) = read_reference_matrix("cooks_replacement_cooks.tsv") else {
        return;
    };
    let Some(design) = read_design_matrix("cooks_replacement_design.tsv") else {
        return;
    };
    let Some((size_factors, replaceable)) = read_replacement_size_factor_rows() else {
        return;
    };
    let Some(options) = read_replacement_options() else {
        return;
    };

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
fn cooks_replacement_refit_plan_matches_optional_deseq2_reference() {
    let Some(counts) = read_count_matrix("cooks_replacement_counts.tsv") else {
        return;
    };
    let Some(rows) = read_optional_tsv("cooks_replacement_rows.tsv") else {
        return;
    };
    let Some(cooks) = read_reference_matrix("cooks_replacement_cooks.tsv") else {
        return;
    };
    let Some(design) = read_design_matrix("cooks_replacement_design.tsv") else {
        return;
    };
    let Some((size_factors, _replaceable)) = read_replacement_size_factor_rows() else {
        return;
    };
    let Some(options) = read_replacement_options() else {
        return;
    };

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
