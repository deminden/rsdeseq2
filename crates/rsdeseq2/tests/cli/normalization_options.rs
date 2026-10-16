use super::numerical_assertions::*;
use super::*;

fn shuffled_assay(source: &Path, target: &Path) {
    let (headers, rows) = read_export(source);
    let order = [0, 3, 1, 4, 2];
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(target)
        .unwrap();
    writer.write_record(order.map(|i| &headers[i])).unwrap();
    for i in [3, 1, 0, 2] {
        writer.write_record(order.map(|j| &rows[i][j])).unwrap();
    }
    writer.flush().unwrap();
}

struct NormalizationCase {
    name: &'static str,
    args: Vec<String>,
    builder: DeseqBuilder,
    commands: &'static [&'static str],
}

fn normalization_cases(dir: &Path) -> Vec<NormalizationCase> {
    let sf = dir.join("size_factors.tsv");
    write_size_factor_fixture(&sf);
    let geo = dir.join("geometric_means.tsv");
    fs::write(
        &geo,
        "gene\tgeo_mean\ngene4\t8\ngene2\t4\ngene1\t3\ngene3\t6\n",
    )
    .unwrap();
    // Gene/sample factors retain overdispersion for native Wald/LRT fitting.
    let nf = dir.join("normalization_factors.tsv");
    let canonical_nf = dir.join("canonical_factors.tsv");
    let factors = RowMajorMatrix::from_row_major(
        4,
        4,
        [1.0, 2.0, 0.5, 1.5]
            .into_iter()
            .flat_map(|gene| [1.0, 0.5, 1.5, 2.0].map(|sample| gene * sample))
            .collect(),
    )
    .unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&canonical_nf)
        .unwrap();
    writer
        .write_record(["gene", "sample1", "sample2", "sample3", "sample4"])
        .unwrap();
    for (i, gene) in common::reference_gene_names().into_iter().enumerate() {
        writer
            .write_record(
                std::iter::once(gene).chain(factors.row(i).unwrap().iter().map(f64::to_string)),
            )
            .unwrap();
    }
    writer.flush().unwrap();
    shuffled_assay(&canonical_nf, &nf);
    let weights = dir.join("weights.tsv");
    shuffled_assay(&reference_data_path("observation_weights.tsv"), &weights);
    let file_args =
        |flag: &str, path: &Path| vec![flag.to_owned(), path.to_str().unwrap().to_owned()];
    vec![
        NormalizationCase {
            name: "control subset",
            args: vec!["--control-genes".into(), "0,3".into()],
            builder: test_builder().control_genes(vec![0, 3]),
            commands: &[
                "size-factors",
                "base-mean",
                "normalized-counts",
                "vst",
                "rlog",
                "wald",
                "lrt",
            ],
        },
        NormalizationCase {
            name: "shuffled geometric means",
            args: file_args("--geometric-means", &geo),
            builder: test_builder().geometric_means(vec![3.0, 4.0, 6.0, 8.0]),
            commands: &[
                "size-factors",
                "base-mean",
                "normalized-counts",
                "vst",
                "rlog",
                "wald",
                "lrt",
            ],
        },
        NormalizationCase {
            name: "shuffled supplied size factors",
            args: file_args("--size-factors", &sf),
            builder: test_builder().size_factors(vec![1.0, 0.5, 1.5, 2.0]),
            commands: &[
                "base-mean",
                "normalized-counts",
                "vst",
                "rlog",
                "wald",
                "lrt",
            ],
        },
        NormalizationCase {
            name: "shuffled normalization factors",
            args: file_args("--normalization-factors", &nf),
            builder: test_builder().normalization_factors(factors),
            commands: &[
                "base-mean",
                "normalized-counts",
                "vst",
                "rlog",
                "wald",
                "lrt",
            ],
        },
        NormalizationCase {
            name: "shuffled nonuniform and zero weights",
            args: file_args("--observation-weights", &weights),
            builder: test_builder().observation_weights(
                rsdeseq2::io::read_observation_weights_tsv(reference_data_path(
                    "observation_weights.tsv",
                ))
                .unwrap(),
            ),
            commands: &["base-mean"],
        },
    ]
}

// Compare cells directly with in-memory values, independently of the TSV writer.
enum ExpectedExport {
    Numbers {
        headers: Vec<String>,
        names: Vec<String>,
        values: Vec<Option<f64>>,
    },
    Results(Box<DeseqResults>),
}

impl ExpectedExport {
    fn values(&self) -> Vec<Option<f64>> {
        match self {
            Self::Numbers { values, .. } => values.clone(),
            Self::Results(results) => result_values(results),
        }
    }

    fn check(&self, path: &Path, label: &str) {
        match self {
            Self::Results(results) => assert_result_export(path, results, label),
            Self::Numbers {
                headers: expected_header,
                names,
                values,
            } => {
                let (headers, rows) = read_export(path);
                assert_eq!(&headers, expected_header, "{label}: header");
                assert_eq!(rows.len(), names.len(), "{label}: row count");
                for (row, name) in rows.iter().zip(names) {
                    assert_eq!(&row[0], name, "{label}: label/order");
                }
                let actual = rows
                    .iter()
                    .flat_map(|row| row.iter().skip(1).map(|s| number(s)))
                    .collect::<Vec<_>>();
                assert_eq!(actual.len(), values.len(), "{label}: number of cells");
                for (a, b) in actual.into_iter().zip(values) {
                    common::assert_option_close(a, *b, 1e-12, 1e-12, label);
                }
            }
        }
    }
}

fn expected_export(
    command: &str,
    builder: &DeseqBuilder,
    design: Option<&DesignMatrix>,
) -> ExpectedExport {
    let counts = common::reference_counts();
    let fit = || builder.fit_size_factors_and_base_means(&counts).unwrap();
    let values = match command {
        "size-factors" => fit().size_factors,
        "base-mean" => fit().base_mean,
        "normalized-counts" => fit().normalized_counts(&counts).unwrap().into_values(),
        "vst" => match design {
            Some(design) => builder.default_vst_glm_mu_auto(&counts, design),
            None => builder.default_blind_vst_glm_mu_auto(&counts),
        }
        .unwrap()
        .transformed
        .into_values(),
        "rlog" => match design {
            Some(design) => builder.rlog_glm_mu(&counts, design),
            None => builder.blind_rlog_glm_mu(&counts),
        }
        .unwrap()
        .transformed
        .into_values(),
        "wald" => {
            return ExpectedExport::Results(Box::new(
                builder
                    .fit_with_results(&counts, &common::reference_full_design())
                    .unwrap()
                    .1,
            ));
        }
        "lrt" => {
            return ExpectedExport::Results(Box::new(
                builder
                    .fit_lrt_with_results(
                        &counts,
                        &common::reference_full_design(),
                        &common::reference_reduced_design(),
                    )
                    .unwrap()
                    .1,
            ));
        }
        _ => unreachable!(),
    };
    let headers = match command {
        "size-factors" => vec!["sample".to_owned(), "size_factor".to_owned()],
        "base-mean" => vec!["gene".to_owned(), "base_mean".to_owned()],
        _ => std::iter::once("gene".to_owned())
            .chain(common::reference_sample_names())
            .collect(),
    };
    let names = if command == "size-factors" {
        common::reference_sample_names()
    } else {
        common::reference_gene_names()
    };
    ExpectedExport::Numbers {
        headers,
        names,
        values: values.into_iter().map(Some).collect(),
    }
}

fn check_normalization_options(command: &str) {
    let dir = temp_dir("normalization-routing-table");
    let cases = normalization_cases(&dir);
    let counts = reference_data_path("counts.tsv");
    let design = reference_data_path("design_full.tsv");
    let reduced = reference_data_path("design_reduced.tsv");
    let baseline = expected_export(command, &test_builder(), None);
    let output = dir.join("results.tsv");
    let mut args = basic_args(command, &counts, &output);
    if !matches!(command, "size-factors" | "base-mean" | "normalized-counts") {
        model_args(&mut args, command, &design, &reduced);
    }
    invoke(&format!("{command} default"), &args);
    baseline.check(&output, command);
    for case in &cases {
        if !case.commands.contains(&command) {
            continue;
        }
        let label = format!("{command}: {}", case.name);
        eprintln!("Preparing {label}");
        let expected = expected_export(command, &case.builder, None);
        assert!(
            numerically_distinct(&expected.values(), &baseline.values()),
            "{label}: fixture cannot detect an ignored flag"
        );
        let mut selected = args.clone();
        selected.extend(case.args.clone());
        invoke(&label, &selected);
        expected.check(&output, &label);
    }
    if matches!(command, "vst" | "rlog") {
        let label = format!("{command}: design-aware transform");
        let expected = expected_export(
            command,
            &test_builder(),
            Some(&common::reference_full_design()),
        );
        assert!(
            numerically_distinct(&expected.values(), &baseline.values()),
            "{label}: fixture cannot detect ignored --blind=false"
        );
        let mut aware = args.clone();
        aware.extend(["--blind=false", "--design", design.to_str().unwrap()].map(str::to_owned));
        invoke(&label, &aware);
        expected.check(&output, &label);
    }
}

#[test]
fn normalization_routes_size_factors() {
    check_normalization_options("size-factors");
}
#[test]
fn normalization_routes_base_means() {
    check_normalization_options("base-mean");
}
#[test]
fn normalization_routes_normalized_counts() {
    check_normalization_options("normalized-counts");
}
#[test]
fn normalization_routes_vst() {
    check_normalization_options("vst");
}
#[test]
fn normalization_routes_rlog() {
    check_normalization_options("rlog");
}
#[test]
fn normalization_routes_wald() {
    check_normalization_options("wald");
}
#[test]
fn normalization_routes_lrt() {
    check_normalization_options("lrt");
}
