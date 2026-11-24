//! Every TSV reader with a `_from_reader` form returns the same value, or the
//! same kind of error, for the same bytes read from a file or from memory.

use std::fmt::Debug;

use rsdeseq2::{errors::DeseqError, io::*};

fn temp_file(name: &str, bytes: &[u8]) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!("rsdeseq2-{}-{name}.tsv", std::process::id()));
    std::fs::write(&path, bytes).unwrap();
    path
}

fn same<T: Debug + PartialEq>(
    name: &str,
    bytes: &[u8],
    from_path: impl Fn(&std::path::Path) -> Result<T, rsdeseq2::errors::DeseqError>,
    from_reader: impl Fn(&[u8]) -> Result<T, rsdeseq2::errors::DeseqError>,
) {
    let path = temp_file(name, bytes);
    let file = from_path(&path);
    let memory = from_reader(bytes);
    std::fs::remove_file(&path).unwrap();
    match (file, memory) {
        (Ok(a), Ok(b)) => assert_eq!(a, b, "{name}"),
        (Err(a), Err(b)) => assert_eq!(a.kind(), b.kind(), "{name}: {a} vs {b}"),
        (a, b) => panic!("{name}: {a:?} vs {b:?}"),
    }
}

#[test]
fn reader_forms_match_path_forms() {
    let counts = b"gene\ts1\ts2\ng1\t1\t2\ng2\t0\t9\n";
    let bad_counts = b"gene\ts1\ts2\ng1\t1\t-2\n";
    for (name, bytes) in [("counts", &counts[..]), ("bad_counts", &bad_counts[..])] {
        same(
            name,
            bytes,
            |p| read_count_matrix_tsv(p),
            |b| read_count_matrix_tsv_from_reader(b),
        );
    }
    let design = b"sample\tintercept\tconditionB\ns1\t1\t0\ns2\t1\t1\n";
    same(
        "design",
        design,
        |p| read_labeled_design_matrix_tsv(p),
        |b| read_labeled_design_matrix_tsv_from_reader(b),
    );
    let values = b"sample\tvalue\ns1\t1.5\ns2\t0.75\n";
    same(
        "size_factors",
        values,
        |p| read_labeled_size_factors_tsv(p),
        |b| read_labeled_size_factors_tsv_from_reader(b),
    );
    same(
        "sample_numeric",
        values,
        |p| read_labeled_sample_numeric_tsv(p, "x"),
        |b| read_labeled_sample_numeric_tsv_from_reader(b, "x"),
    );
    same(
        "wald_t",
        values,
        |p| read_labeled_wald_t_degrees_of_freedom_tsv(p),
        |b| read_labeled_wald_t_degrees_of_freedom_tsv_from_reader(b),
    );
    let genes = b"gene\tvalue\ng1\t2.5\ng2\t3\n";
    same(
        "geometric_means",
        genes,
        |p| read_labeled_geometric_means_tsv(p),
        |b| read_labeled_geometric_means_tsv_from_reader(b),
    );
    same(
        "gene_numeric",
        genes,
        |p| read_labeled_gene_numeric_tsv(p, "x"),
        |b| read_labeled_gene_numeric_tsv_from_reader(b, "x"),
    );
    let levels = b"sample\tcondition\ns1\tA\ns2\tB\n";
    same(
        "levels",
        levels,
        |p| read_sample_levels_tsv(p),
        |b| read_sample_levels_tsv_from_reader(b),
    );
}

#[test]
fn reused_matrix_records_keep_error_coordinates() {
    let counts = b"gene\ts1\ts2\ng1\t1\t2\ng2\t3\n";
    let error = read_count_matrix_tsv_from_reader(&counts[..]).unwrap_err();
    std::assert_matches!(
        error,
        DeseqError::InvalidTableWidth {
            line: 3,
            expected: 3,
            actual: 2
        }
    );

    let design = b"sample\tintercept\tcondition\n\"s\t1\"\t1\t0\ns2\t1\tNaN\n";
    let error = read_labeled_design_matrix_tsv_from_reader(&design[..]).unwrap_err();
    std::assert_matches!(
        error,
        DeseqError::NonFiniteValue { context, index: Some(3), value }
            if context == "design matrix" && value.is_nan()
    );
}

#[test]
fn quoted_records_cross_buffer_boundaries_and_validate_field_encoding() {
    let label = format!("{}\tline\n\"quoted\"", "é".repeat(5000));
    let escaped = label.replace('"', "\"\"");
    let input = format!("sample\tintercept\n\"{escaped}\"\t1\n");
    let result = read_labeled_design_matrix_tsv_from_reader(input.as_bytes()).unwrap();
    assert_eq!(result.sample_names, [label]);

    // These bytes form a valid character only if the field delimiter is ignored.
    let bad = b"sample\tintercept\n\xc3\t\xa9\n";
    std::assert_matches!(
        read_labeled_design_matrix_tsv_from_reader(&bad[..]),
        Err(DeseqError::InvalidTableEncoding { line: 2 })
    );
}
