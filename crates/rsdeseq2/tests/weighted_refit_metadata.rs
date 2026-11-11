//! Original replacement counts and raw-weight base metadata, independent of optimizers.
use rsdeseq2::{cooks::prepare_cooks_replacement_refit_with_weights, prelude::*};
use serde_json::Value;
fn number(v: &Value) -> f64 {
    let s = v.as_str().unwrap();
    if let Some(h) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(h, 16).unwrap())
    } else {
        s.parse().unwrap()
    }
}
fn matrix(v: &Value) -> RowMajorMatrix<f64> {
    let r = v["rows"].as_u64().unwrap() as usize;
    let c = v["columns"].as_u64().unwrap() as usize;
    let values = v["values"]
        .as_array()
        .unwrap()
        .iter()
        .map(number)
        .collect::<Vec<_>>();
    RowMajorMatrix::from_row_major(
        r,
        c,
        (0..r)
            .flat_map(|i| {
                (0..c).map({
                    let values = &values;
                    move |j| values[i + j * r]
                })
            })
            .collect(),
    )
    .unwrap()
}
#[test]
fn replacement_metadata_matches_original_deseq2_with_raw_weights_and_both_normalizations() {
    let reference: Value =
        serde_json::from_str(include_str!("data/weighted_refit_metadata_r462.json")).unwrap();
    for c in reference["cases"].as_array().unwrap() {
        let counts = matrix(&c["counts"]);
        let counts = CountMatrix::from_row_major_u32(
            counts.n_rows(),
            counts.n_cols(),
            counts.as_slice().iter().map(|v| *v as u32).collect(),
        )
        .unwrap();
        let m = matrix(&c["design"]);
        let design =
            DesignMatrix::from_row_major(m.n_rows(), m.n_cols(), m.into_values(), None).unwrap();
        let sf = c["size_factors"]
            .as_array()
            .unwrap()
            .iter()
            .map(number)
            .collect::<Vec<_>>();
        let nf =
            (!c["normalization_factors"].is_null()).then(|| matrix(&c["normalization_factors"]));
        let weights = (!c["weights"].is_null()).then(|| matrix(&c["weights"]));
        let normalized = if let Some(nf) = &nf {
            normalized_counts_with_factors(&counts, nf).unwrap()
        } else {
            normalized_counts(&counts, &sf).unwrap()
        };
        let mut options = CooksReplacementOptions::new(number(&c["cutoff"]));
        options.min_replicates = 3;
        options.which_samples = Some(
            c["which_samples"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_bool().unwrap())
                .collect(),
        );
        let plan = prepare_cooks_replacement_refit_with_weights(
            &counts,
            &normalized,
            &sf,
            nf.as_ref(),
            &matrix(&c["cooks"]),
            &design,
            &options,
            weights.as_ref(),
        )
        .unwrap();
        assert_eq!(
            plan.replacement
                .replaced_counts
                .as_slice()
                .iter()
                .map(|v| f64::from(*v))
                .collect::<Vec<_>>(),
            matrix(&c["replaced_counts"]).as_slice(),
            "{} counts",
            c["name"]
        );
        assert_eq!(plan.n_refit, c["n_refit"].as_u64().unwrap() as usize);
        assert_eq!(
            plan.replaced_all_zero,
            c["all_zero"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_bool().unwrap())
                .collect::<Vec<_>>()
        );
        for (field, actual) in [
            ("base_mean", &plan.replaced_base_mean),
            ("base_var", &plan.replaced_base_var),
        ] {
            for (i, (a, b)) in actual
                .iter()
                .zip(c[field].as_array().unwrap().iter().map(number))
                .enumerate()
            {
                assert!(
                    a.is_finite() && (a - b).abs() <= 1e-10 + 1e-12 * b.abs(),
                    "{} {field} row {i}: {a} != {b}",
                    c["name"]
                );
            }
        }
    }
}
