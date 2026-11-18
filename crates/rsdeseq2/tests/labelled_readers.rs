//! Error messages of the two-column labelled-table readers, which share one
//! parser. Each reader must keep its own wording for every failure.

use rsdeseq2::errors::DeseqError;
use rsdeseq2::io::{
    read_labeled_gene_numeric_tsv_from_reader, read_labeled_geometric_means_tsv_from_reader,
    read_labeled_sample_numeric_tsv_from_reader, read_labeled_size_factors_tsv_from_reader,
    read_labeled_wald_t_degrees_of_freedom_tsv_from_reader, read_sample_levels_tsv_from_reader,
};

type Reader = fn(&[u8]) -> Result<usize, DeseqError>;

fn readers() -> Vec<(&'static str, &'static str, Reader)> {
    vec![
        (
            "size-factor",
            "size-factor sample names must not be empty",
            |b| read_labeled_size_factors_tsv_from_reader(b).map(|v| v.len()),
        ),
        (
            "geometric-mean",
            "geometric-mean gene names must not be empty",
            |b| read_labeled_geometric_means_tsv_from_reader(b).map(|v| v.len()),
        ),
        ("offsets", "offsets sample names must not be empty", |b| {
            read_labeled_sample_numeric_tsv_from_reader(b, "offsets").map(|v| v.len())
        }),
        (
            "dispersions",
            "dispersions gene names must not be empty",
            |b| read_labeled_gene_numeric_tsv_from_reader(b, "dispersions").map(|v| v.len()),
        ),
        (
            "Wald t degrees-of-freedom",
            "Wald t degrees-of-freedom gene names must not be empty",
            |b| read_labeled_wald_t_degrees_of_freedom_tsv_from_reader(b).map(|v| v.len()),
        ),
        ("sample-level", "sample names must not be empty", |b| {
            read_sample_levels_tsv_from_reader(b).map(|v| v.len())
        }),
    ]
}

fn dimensions(error: DeseqError) -> (String, usize, usize) {
    match error {
        DeseqError::InvalidDimensions {
            context,
            expected,
            actual,
        } => (context, expected, actual),
        other => panic!("expected InvalidDimensions, got {other:?}"),
    }
}

#[test]
fn labelled_readers_keep_their_error_wording() {
    for (what, empty_label, read) in readers() {
        assert_eq!(
            dimensions(read(b"a\tb\tc\nx\t1\t2\n").unwrap_err()),
            (format!("{what} columns"), 2, 3)
        );
        // A ragged row is rejected by the TSV parser itself.
        assert!(matches!(read(b"a\tb\nx\t1\t2\n"), Err(DeseqError::Csv(_))));
        assert_eq!(
            dimensions(read(b"a\tb\n").unwrap_err()),
            (format!("{what} rows"), 1, 0)
        );
        match read(b"a\tb\n\t1\n").unwrap_err() {
            DeseqError::InvalidOptions { reason } => assert_eq!(reason, empty_label),
            other => panic!("{what}: {other:?}"),
        }
        assert_eq!(read(b"a\tb\nx\t1\ny\t2\n").unwrap(), 2, "{what}");
    }
}

#[test]
fn labelled_readers_keep_their_value_rules() {
    let size_factors = read_labeled_size_factors_tsv_from_reader(&b"s\tv\nx\t0\n"[..]);
    assert!(size_factors.is_err());
    let geometric = read_labeled_geometric_means_tsv_from_reader(&b"g\tv\nx\t-1\n"[..]);
    assert!(matches!(
        geometric,
        Err(DeseqError::InvalidSizeFactors { .. }
            | DeseqError::InvalidParameter {
                domain: rsdeseq2::errors::ErrorKind::Normalization,
                ..
            })
    ));
    let levels = read_sample_levels_tsv_from_reader(&b"s\tl\nx\t\n"[..]).unwrap_err();
    assert!(
        levels
            .to_string()
            .contains("sample levels must not be empty")
    );
    let numeric = read_labeled_gene_numeric_tsv_from_reader(&b"g\tv\nx\tnope\n"[..], "dispersions");
    assert!(numeric.is_err());
}
