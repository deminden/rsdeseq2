use super::numerical_assertions::*;
use super::*;

struct ContrastFixture {
    counts: CountMatrix,
    design: DesignMatrix,
    reduced: DesignMatrix,
    counts_path: PathBuf,
    design_path: PathBuf,
    reduced_path: PathBuf,
}

fn contrast_fixture(dir: &Path, effect_name: &str) -> ContrastFixture {
    let counts = CountMatrix::from_row_major_u32_with_names(
        4,
        8,
        vec![
            10, 20, 40, 80, 30, 10, 18, 160, 80, 16, 20, 4, 22, 46, 8, 18, 50, 40, 60, 44, 120, 9,
            170, 15, 6, 12, 4, 8, 30, 4, 18, 3,
        ],
        Some(common::reference_gene_names()),
        Some((1..=8).map(|i| format!("sample{i}")).collect()),
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        3,
        (0..8)
            .flat_map(|i| [1.0, if i % 4 >= 2 { 1.0 } else { 0.0 }, (i % 2) as f64])
            .collect(),
        Some(vec![
            "Intercept".into(),
            effect_name.into(),
            "batchY".into(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::intercept_only(8).unwrap();
    let counts_path = dir.join("counts.tsv");
    let design_path = dir.join("design.tsv");
    let reduced_path = dir.join("reduced.tsv");
    let samples = counts.sample_names().unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&counts_path)
        .unwrap();
    writer
        .write_record(std::iter::once("gene").chain(samples.iter().map(String::as_str)))
        .unwrap();
    for (i, name) in counts.gene_names().unwrap().iter().enumerate() {
        writer
            .write_record(
                std::iter::once(name.clone())
                    .chain(counts.row(i).unwrap().iter().map(u32::to_string)),
            )
            .unwrap();
    }
    writer.flush().unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&design_path)
        .unwrap();
    writer
        .write_record(["sample", "Intercept", effect_name, "batchY"])
        .unwrap();
    let mut reduced_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&reduced_path)
        .unwrap();
    reduced_writer
        .write_record(["sample", "Intercept"])
        .unwrap();
    // Independent order in both design files catches ignored sample labels.
    for i in [3, 0, 7, 2, 6, 1, 5, 4] {
        writer
            .write_record([
                samples[i].clone(),
                "1".into(),
                (if i % 4 >= 2 { "1" } else { "0" }).into(),
                (i % 2).to_string(),
            ])
            .unwrap();
        reduced_writer
            .write_record([samples[i].as_str(), "1"])
            .unwrap();
    }
    writer.flush().unwrap();
    reduced_writer.flush().unwrap();
    ContrastFixture {
        counts,
        design,
        reduced,
        counts_path,
        design_path,
        reduced_path,
    }
}

fn contrast_result(command: &str, fixture: &ContrastFixture, contrast: &[f64]) -> DeseqResults {
    if command == "wald" {
        test_builder().fit_with_results_contrast(&fixture.counts, &fixture.design, contrast)
    } else {
        test_builder().fit_lrt_with_results_contrast(
            &fixture.counts,
            &fixture.design,
            &fixture.reduced,
            contrast,
        )
    }
    .unwrap()
    .1
}

fn check_coefficient_options(command: &str) {
    for effect_name in ["conditionB", "if."] {
        let dir = temp_dir("contrast-routing-table");
        let fixture = contrast_fixture(&dir, effect_name);
        let first = contrast_result(command, &fixture, &[0.0, 1.0, 0.0]);
        let last = contrast_result(command, &fixture, &[0.0, 0.0, 1.0]);
        assert!(
            numerically_distinct(&result_values(&first), &result_values(&last)),
            "{command}: selectors must produce distinct numbers"
        );
        let mut cases = Vec::new();
        for (index, name, contrast) in [
            ("1", effect_name, [0.0, 1.0, 0.0]),
            ("2", "batchY", [0.0, 0.0, 1.0]),
        ] {
            for (flag, value) in [
                ("--coefficient", index),
                ("--coefficient-name", name),
                ("--contrast-name", name),
            ] {
                // Numeric indices and the second ordinary name already have
                // coverage above; reserved names need their own name lookups.
                if effect_name == "if." && (flag == "--coefficient" || index == "2") {
                    continue;
                }
                cases.push((
                    format!("{flag} {value}"),
                    vec![flag.to_owned(), value.to_owned()],
                    contrast,
                    None,
                ));
            }
        }
        if effect_name == "if." {
            for flag in ["--coefficient-name", "--contrast-name"] {
                cases.push((
                    format!("{flag} cleaned reserved alias"),
                    vec![flag.into(), "if".into()],
                    [0.0, 1.0, 0.0],
                    None,
                ));
            }
        }
        if effect_name != "if." {
            cases.push((
                "numeric contrast".into(),
                vec!["--contrast=0,1,-0.5".into()],
                [0.0, 1.0, -0.5],
                Some("primitive numeric contrast".to_owned()),
            ));
        }
        let list_name = if effect_name == "if." {
            "if"
        } else {
            effect_name
        };
        cases.push((
            "weighted list contrast".into(),
            vec![
                "--contrast-positive".into(),
                list_name.into(),
                "--contrast-negative".into(),
                "(Intercept)".into(),
                "--contrast-positive-weight".into(),
                "1".into(),
                "--contrast-negative-weight=-0.5".into(),
            ],
            [-0.5, 1.0, 0.0],
            Some(format!(
                "coefficient list contrast: {list_name} vs 0.5 (Intercept)"
            )),
        ));
        for (name, flags, contrast, comparison) in cases {
            let label = format!("{command}/{effect_name}/{name}");
            let output = dir.join("results.tsv");
            let metadata = dir.join("metadata.tsv");
            let mut args = basic_args(command, &fixture.counts_path, &output);
            model_args(
                &mut args,
                command,
                &fixture.design_path,
                &fixture.reduced_path,
            );
            args.extend([
                "--result-table-metadata-output".to_owned(),
                metadata.to_str().unwrap().to_owned(),
            ]);
            args.extend(flags);
            invoke(&label, &args);
            let expected = if contrast == [0.0, 1.0, 0.0] {
                first.clone()
            } else if contrast == [0.0, 0.0, 1.0] {
                last.clone()
            } else {
                contrast_result(command, &fixture, &contrast)
            };
            assert_result_export(&output, &expected, &label);
            assert_tsv_contains(
                &metadata,
                &format!(
                    "testType\t{}\n",
                    if command == "wald" { "Wald" } else { "LRT" }
                ),
            );
            if let Some(comparison) = comparison {
                assert_tsv_contains(&metadata, "resultName\tcontrast\n");
                assert_tsv_contains(&metadata, &format!("comparison\t{comparison}\n"));
                let vector = contrast
                    .iter()
                    .map(f64::to_string)
                    .collect::<Vec<_>>()
                    .join(",");
                assert_tsv_contains(&metadata, &format!("contrast\t{vector}\n"));
                assert!(
                    numerically_distinct(&result_values(&expected), &result_values(&last)),
                    "{label}: contrast must differ from default"
                );
            }
        }
    }
}

fn check_factor_options(command: &str) {
    for (effect, factor, numerator, denominator, reference, raw_a, raw_b, comparison) in [
        (
            "condition_B_vs_A",
            "condition",
            "B",
            "A",
            None,
            "A",
            "B",
            "factor-level contrast: condition B vs A",
        ),
        (
            "condition_if._vs_TRUE.",
            "condition",
            "if",
            "TRUE",
            None,
            "TRUE",
            "if",
            "factor-level contrast: condition if vs TRUE",
        ),
        (
            "cell.type_B.cell_vs_A.cell",
            "cell.type",
            "B.cell",
            "A.cell",
            Some("A.cell"),
            "A cell",
            "B cell",
            "factor-level contrast: cell.type B cell vs A cell",
        ),
    ] {
        let dir = temp_dir("factor-routing-table");
        let fixture = contrast_fixture(&dir, effect);
        let levels = dir.join("levels.tsv");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&levels)
            .unwrap();
        writer.write_record(["sample", factor]).unwrap();
        for i in [6, 1, 4, 0, 7, 2, 5, 3] {
            writer
                .write_record([
                    format!("sample{}", i + 1),
                    if i % 4 >= 2 {
                        raw_b.to_owned()
                    } else {
                        raw_a.to_owned()
                    },
                ])
                .unwrap();
        }
        writer.flush().unwrap();
        let expected = contrast_result(command, &fixture, &[0.0, 1.0, 0.0]);
        let default = contrast_result(command, &fixture, &[0.0, 0.0, 1.0]);
        assert!(numerically_distinct(
            &result_values(&expected),
            &result_values(&default)
        ));
        let output = dir.join("results.tsv");
        let metadata = dir.join("metadata.tsv");
        let mut args = basic_args(command, &fixture.counts_path, &output);
        model_args(
            &mut args,
            command,
            &fixture.design_path,
            &fixture.reduced_path,
        );
        args.extend(
            [
                "--contrast-factor",
                factor,
                "--contrast-numerator",
                numerator,
                "--contrast-denominator",
                denominator,
                "--contrast-sample-levels",
                levels.to_str().unwrap(),
                "--result-table-metadata-output",
                metadata.to_str().unwrap(),
            ]
            .map(str::to_owned),
        );
        if let Some(reference) = reference {
            args.extend(["--contrast-reference", reference].map(str::to_owned));
        }
        let label = format!("{command}/{effect}");
        invoke(&label, &args);
        assert_result_export(&output, &expected, &label);
        assert_tsv_contains(&metadata, &format!("comparison\t{comparison}\n"));
        assert_tsv_contains(&metadata, &format!("resultName\t{effect}\n"));
        assert_tsv_contains(&metadata, "contrast\t0,1,0\n");
    }
}

#[test]
fn wald_coefficient_names_aliases_and_contrasts() {
    check_coefficient_options("wald");
}
#[test]
fn lrt_coefficient_names_aliases_and_contrasts() {
    check_coefficient_options("lrt");
}
#[test]
fn wald_factor_levels_and_cleaned_references() {
    check_factor_options("wald");
}
#[test]
fn lrt_factor_levels_and_cleaned_references() {
    check_factor_options("lrt");
}
