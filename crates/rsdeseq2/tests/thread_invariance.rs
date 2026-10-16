mod common;
use rsdeseq2::prelude::*;

fn assert_numbers(a: &[f64], b: &[f64], label: &str) {
    assert_eq!(a.len(), b.len(), "{label}");
    for (i, (a, b)) in a.iter().zip(b).enumerate() {
        assert!(
            a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
            "{label}[{i}]: {a} != {b}"
        );
    }
}

fn assert_fit_equal(a: &DeseqFit, b: &DeseqFit) {
    assert_numbers(&a.size_factors, &b.size_factors, "size factors");
    assert_numbers(&a.base_mean, &b.base_mean, "base means");
    assert_numbers(&a.base_var, &b.base_var, "base variance");
    macro_rules! vectors {
        ($($field:ident),*) => {$(
            match (&a.$field, &b.$field) {
                (Some(a),Some(b)) => assert_numbers(a,b,stringify!($field)),
                (None,None) => {},
                _ => panic!("{} missingness mismatch", stringify!($field)),
            }
        )*};
    }
    vectors!(
        disp_gene_est,
        disp_fit,
        disp_map,
        dispersion,
        beta_optim_iter,
        beta_optim_start_objective,
        beta_optim_objective,
        beta_optim_gradient_norm,
        log_like,
        full_deviance,
        reduced_log_like
    );
    macro_rules! matrices {
        ($($field:ident),*) => {$(
            match (&a.$field, &b.$field) {
                (Some(a),Some(b)) => {
                    assert_eq!((a.n_rows(),a.n_cols()),(b.n_rows(),b.n_cols()));
                    assert_numbers(a.as_slice(),b.as_slice(),stringify!($field));
                },
                (None,None) => {},
                _ => panic!("{} missingness mismatch", stringify!($field)),
            }
        )*};
    }
    matrices!(
        beta,
        beta_se,
        beta_covariance,
        beta_optim_start,
        mu,
        hat_diagonal,
        cooks,
        reduced_mu,
        reduced_hat_diagonal,
        observation_weights
    );
    assert_eq!(a.all_zero, b.all_zero);
    assert_eq!(a.weights_fail, b.weights_fail);
    assert_eq!(a.disp_gene_iter, b.disp_gene_iter);
    assert_eq!(a.disp_iter, b.disp_iter);
    assert_eq!(a.disp_outlier, b.disp_outlier);
    assert_eq!(a.dispersion_converged, b.dispersion_converged);
    assert_eq!(a.beta_iter, b.beta_iter);
    assert_eq!(a.beta_converged, b.beta_converged);
    assert_eq!(a.reduced_beta_iter, b.reduced_beta_iter);
    assert_eq!(a.reduced_beta_converged, b.reduced_beta_converged);
}

#[test]
fn native_wald_and_lrt_are_identical_across_thread_counts_with_and_without_weights() {
    let dir =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gtex_auto_dispersion");
    let counts = rsdeseq2::io::read_count_matrix_tsv(dir.join("counts.tsv")).unwrap();
    let design = rsdeseq2::io::read_design_matrix_tsv(dir.join("design_matrix.tsv")).unwrap();
    let reduced = DesignMatrix::from_row_major(
        counts.n_samples(),
        1,
        vec![1.0; counts.n_samples()],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();
    let factors = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("reference_size_factors.tsv"))
        .unwrap()
        .records()
        .map(|r| r.unwrap()[1].parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    for weighted in [false, true] {
        let mut builder = DeseqBuilder::new()
            .size_factors(factors.clone())
            .fit_type(FitType::Mean)
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        if weighted {
            let values = (0..counts.n_genes() * counts.n_samples())
                .map(|i| {
                    if i % 17 == 0 {
                        0.0
                    } else {
                        0.25 + (i % 7) as f64 / 8.0
                    }
                })
                .collect();
            builder = builder.observation_weights(
                RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
                    .unwrap(),
            );
        }
        for lrt in [false, true] {
            let run = |threads| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                pool.install(|| {
                    assert_eq!(rayon::current_num_threads(), threads);
                    let b = builder.clone().threads(threads);
                    let normalization = DeseqBuilder::new()
                        .threads(threads)
                        .fit_size_factors_and_base_means(&common::reference_counts())
                        .unwrap();
                    let fitted = if lrt {
                        b.fit_lrt_with_results(&counts, &design, &reduced)
                    } else {
                        b.fit_with_results(&counts, &design)
                    }
                    .unwrap();
                    (normalization, fitted)
                })
            };
            let (initial_a, (a, ar)) = run(1);
            let (initial_b, (b, br)) = run(4);
            assert_fit_equal(&initial_a, &initial_b);
            assert_fit_equal(&a, &b);
            // Results use None for missing fields, so structural equality also
            // checks all six numerical fields, filtering and Cook's metadata.
            assert_eq!(ar, br, "weighted={weighted}, lrt={lrt}");
            assert!(ar.rows.iter().any(|r| r.pvalue.is_some()));
        }
    }
}
