use std::path::Path;

use csv::{ReaderBuilder, WriterBuilder};

use crate::core::CountMatrix;
use crate::errors::DeseqError;

/// Read a tab-delimited count matrix with gene IDs in the first column.
pub fn read_count_matrix_tsv(path: impl AsRef<Path>) -> Result<CountMatrix, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidCounts {
            reason: "count table must have a gene column and at least one sample column"
                .to_string(),
        });
    }
    let sample_names = headers
        .iter()
        .skip(1)
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    let mut gene_names = Vec::new();
    let mut values = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidCounts {
                reason: format!(
                    "record has {} columns but header has {} columns",
                    record.len(),
                    headers.len()
                ),
            });
        }
        gene_names.push(record.get(0).unwrap_or_default().to_string());
        for field in record.iter().skip(1) {
            values.push(parse_count_field(field, path_ref)?);
        }
    }
    CountMatrix::from_row_major_u32_with_names(
        gene_names.len(),
        sample_names.len(),
        values,
        Some(gene_names),
        Some(sample_names),
    )
}

fn parse_count_field(field: &str, path: &Path) -> Result<u32, DeseqError> {
    if let Ok(value) = field.parse::<u32>() {
        return Ok(value);
    }

    let value = field.parse::<f64>().map_err(|_| DeseqError::ParseInt {
        context: path.display().to_string(),
        value: field.to_string(),
    })?;
    if !value.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: path.display().to_string(),
            index: None,
            value,
        });
    }
    if value < 0.0 || value.fract() != 0.0 || value > f64::from(u32::MAX) {
        return Err(DeseqError::InvalidCounts {
            reason: format!("count value '{field}' must be a non-negative integer"),
        });
    }
    Ok(value as u32)
}

/// Write sample-level size factors to a tab-delimited file.
pub fn write_size_factors_tsv(
    path: impl AsRef<Path>,
    sample_names: Option<&[String]>,
    size_factors: &[f64],
) -> Result<(), DeseqError> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["sample", "size_factor"])?;
    for (idx, size_factor) in size_factors.iter().copied().enumerate() {
        let sample = sample_names
            .and_then(|names| names.get(idx))
            .cloned()
            .unwrap_or_else(|| format!("sample{}", idx + 1));
        writer.write_record([sample, size_factor.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write gene-level base means to a tab-delimited file.
pub fn write_base_mean_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    base_mean: &[f64],
) -> Result<(), DeseqError> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["gene", "base_mean"])?;
    for (idx, value) in base_mean.iter().copied().enumerate() {
        let gene = gene_names
            .and_then(|names| names.get(idx))
            .cloned()
            .unwrap_or_else(|| format!("gene{}", idx + 1));
        writer.write_record([gene, value.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_count_field_accepts_integer_scientific_notation() {
        let path = Path::new("counts.tsv");
        assert_eq!(parse_count_field("1e+05", path).unwrap(), 100_000);
        assert_eq!(parse_count_field("42.0", path).unwrap(), 42);
    }

    #[test]
    fn parse_count_field_rejects_non_integer_numeric_values() {
        let path = Path::new("counts.tsv");
        assert!(matches!(
            parse_count_field("1.5", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
        assert!(matches!(
            parse_count_field("-1", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
    }
}
