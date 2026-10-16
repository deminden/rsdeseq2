use super::numerical_assertions::*;
use super::*;

fn filter_fixture(dir: &Path) -> (CountMatrix, DesignMatrix, DesignMatrix) {
    // Four replicates per group make Cook's masking eligible. The outlier and
    // separated gene means distinguish cutoff and independent-filter options.
    let counts = CountMatrix::from_row_major_u32_with_names(
        6,
        8,
        vec![
            100_000, 12, 18, 10, 30, 20, 40, 24, 4, 8, 3, 6, 40, 18, 60, 24, 80, 120, 60, 100, 8,
            16, 5, 12, 40, 80, 20, 60, 70, 30, 100, 40, 400, 800, 200, 600, 700, 300, 1000, 400, 1,
            3, 0, 2, 4, 1, 5, 2,
        ],
        Some((1..=6).map(|i| format!("gene{i}")).collect()),
        Some((1..=8).map(|i| format!("sample{i}")).collect()),
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        (0..8)
            .flat_map(|i| [1.0, if i >= 4 { 1.0 } else { 0.0 }])
            .collect(),
        Some(vec!["Intercept".into(), "conditionB".into()]),
    )
    .unwrap();
    let reduced = DesignMatrix::intercept_only(8).unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("counts.tsv"))
        .unwrap();
    writer
        .write_record(
            std::iter::once("gene")
                .chain(counts.sample_names().unwrap().iter().map(String::as_str)),
        )
        .unwrap();
    for (i, gene) in counts.gene_names().unwrap().iter().enumerate() {
        writer
            .write_record(
                std::iter::once(gene.clone())
                    .chain(counts.row(i).unwrap().iter().map(u32::to_string)),
            )
            .unwrap();
    }
    writer.flush().unwrap();
    for (filename, full) in [("design.tsv", true), ("reduced.tsv", false)] {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(dir.join(filename))
            .unwrap();
        writer
            .write_record(if full {
                vec!["sample", "Intercept", "conditionB"]
            } else {
                vec!["sample", "Intercept"]
            })
            .unwrap();
        for i in [4, 0, 6, 2, 7, 1, 5, 3] {
            let mut row = vec![format!("sample{}", i + 1), "1".into()];
            if full {
                row.push(if i >= 4 { "1" } else { "0" }.into());
            }
            writer.write_record(row).unwrap();
        }
        writer.flush().unwrap();
    }
    (counts, design, reduced)
}

fn assert_filter_exports(dir: &Path, expected: &DeseqResults, label: &str) {
    let filtering = expected.independent_filtering.as_ref().unwrap();
    let (headers, rows) = read_export(&dir.join("metadata.tsv"));
    assert_eq!(headers, ["name", "value"]);
    let scalars = filtering.scalar_metadata();
    assert_eq!(rows.len(), scalars.len(), "{label}: scalar metadata");
    for (row, expected) in rows.iter().zip(scalars) {
        assert_eq!(row[0], expected.name, "{label}: scalar name");
        common::assert_option_close(number(&row[1]), Some(expected.value), 1e-12, 1e-12, label);
    }
    let (headers, rows) = read_export(&dir.join("rejections.tsv"));
    assert_eq!(headers, ["theta", "numRej"]);
    assert_eq!(rows.len(), filtering.theta.len());
    for (i, row) in rows.iter().enumerate() {
        common::assert_option_close(
            number(&row[0]),
            Some(filtering.theta[i]),
            1e-12,
            1e-12,
            label,
        );
        assert_eq!(
            row[1].parse::<usize>().unwrap(),
            filtering.num_rejections[i],
            "{label}: rejections"
        );
    }
    let (headers, rows) = read_export(&dir.join("lowess.tsv"));
    assert_eq!(headers, ["x", "y"]);
    let lowess = filtering.lowess_fit_table();
    assert_eq!(rows.len(), lowess.len());
    for (row, expected) in rows.iter().zip(lowess) {
        common::assert_option_close(number(&row[0]), Some(expected.theta), 1e-12, 1e-12, label);
        common::assert_option_close(
            number(&row[1]),
            Some(expected.fitted_rejections),
            1e-12,
            1e-12,
            label,
        );
    }
}

#[test]
fn wald_and_lrt_filter_controls_match_library_results_and_metadata() {
    let dir = temp_dir("filter-options-numerical");
    let (counts, design, reduced) = filter_fixture(&dir);
    for command in ["wald", "lrt"] {
        let result = |builder: &DeseqBuilder| {
            if command == "wald" {
                builder.fit_with_results(&counts, &design)
            } else {
                builder.fit_lrt_with_results(&counts, &design, &reduced)
            }
            .unwrap()
            .1
        };
        let base = test_builder();
        let unfiltered = result(&base);
        let forced = IndependentFilteringOptions {
            enabled: true,
            alpha: 0.05,
            theta: Some(vec![1.0, 1.0]),
        };
        let cases = [
            (
                "explicit cutoff and disabled filtering",
                vec![
                    "--cooks-cutoff",
                    "10",
                    "--disable-independent-filtering",
                    "--independent-filtering-alpha",
                    "0.05",
                    "--independent-filtering-theta",
                    "0,0.5,1",
                ],
                base.clone()
                    .cooks_cutoff_threshold(10.0)
                    .independent_filtering_alpha(0.05)
                    .independent_filtering_theta(vec![0.0, 0.5, 1.0]),
            ),
            (
                "disabled Cook's and filtering",
                vec![
                    "--disable-cooks-cutoff",
                    "--disable-independent-filtering",
                    "--independent-filtering-theta",
                    "0,0.5,1",
                ],
                base.clone()
                    .independent_filtering_theta(vec![0.0, 0.5, 1.0]),
            ),
            (
                "forced independent filtering",
                vec![
                    "--disable-cooks-cutoff",
                    "--independent-filtering-alpha",
                    "0.05",
                    "--independent-filtering-theta",
                    "1,1",
                ],
                base.clone().independent_filtering_options(forced),
            ),
            (
                "zero Cook's cutoff",
                vec!["--cooks-cutoff", "0", "--disable-independent-filtering"],
                base.clone().cooks_cutoff_threshold(0.0),
            ),
            (
                "active filtering grid and alpha",
                vec![
                    "--disable-cooks-cutoff",
                    "--independent-filtering-alpha",
                    "0.5",
                    "--independent-filtering-theta",
                    "0,0.5,1",
                ],
                base.clone()
                    .independent_filtering_options(IndependentFilteringOptions {
                        enabled: true,
                        alpha: 0.5,
                        theta: Some(vec![0.0, 0.5, 1.0]),
                    }),
            ),
        ];
        for (name, flags, builder) in cases {
            let label = format!("{command}: {name}");
            let expected = result(&builder);
            if name == "forced independent filtering" || name == "zero Cook's cutoff" {
                assert!(
                    numerically_distinct(&result_values(&expected), &result_values(&unfiltered)),
                    "{label}: ignored option undetectable"
                );
            }
            let output = dir.join("results.tsv");
            let mut args = basic_args(command, &dir.join("counts.tsv"), &output);
            args.extend([
                "--design".into(),
                dir.join("design.tsv").to_str().unwrap().into(),
                "--fit-type".into(),
                "mean".into(),
                "--coefficient".into(),
                "1".into(),
            ]);
            if command == "lrt" {
                args.extend([
                    "--reduced-design".into(),
                    dir.join("reduced.tsv").to_str().unwrap().into(),
                ]);
            }
            for (flag, filename) in [
                ("--independent-filter-metadata-output", "metadata.tsv"),
                ("--independent-filter-num-rej-output", "rejections.tsv"),
                ("--independent-filter-lowess-output", "lowess.tsv"),
            ] {
                args.extend([flag.into(), dir.join(filename).to_str().unwrap().into()]);
            }
            args.extend(flags.into_iter().map(str::to_owned));
            invoke(&label, &args);
            assert_result_export(&output, &expected, &label);
            assert_filter_exports(&dir, &expected, &label);
        }
    }
}
