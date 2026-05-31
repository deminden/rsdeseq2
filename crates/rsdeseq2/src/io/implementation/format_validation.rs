fn diagnostic_frame_row_count(
    columns: &[Deseq2McolsDiagnosticColumn],
) -> Result<usize, DeseqError> {
    let Some(first) = columns.first() else {
        return Ok(0);
    };
    let expected = diagnostic_column_len(first);
    for column in columns.iter().skip(1) {
        let actual = diagnostic_column_len(column);
        if actual != expected {
            return Err(invalid_dimensions(
                format!("diagnostic column {}", column.name),
                expected,
                actual,
            ));
        }
    }
    Ok(expected)
}

fn diagnostic_column_len(column: &Deseq2McolsDiagnosticColumn) -> usize {
    match &column.values {
        Deseq2McolsDiagnosticValues::Numeric(values) => values.len(),
        Deseq2McolsDiagnosticValues::OptionalNumeric(values) => values.len(),
        Deseq2McolsDiagnosticValues::Integer(values) => values.len(),
        Deseq2McolsDiagnosticValues::Logical(values) => values.len(),
    }
}

fn format_result_column_value(values: &DeseqResultColumnValues, row_idx: usize) -> String {
    match values {
        DeseqResultColumnValues::Numeric(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        DeseqResultColumnValues::Logical(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .map(|value| {
                if value {
                    "TRUE".to_string()
                } else {
                    "FALSE".to_string()
                }
            })
            .unwrap_or_else(|| "NA".to_string()),
    }
}

fn format_diagnostic_column_value(values: &Deseq2McolsDiagnosticValues, row_idx: usize) -> String {
    match values {
        Deseq2McolsDiagnosticValues::Numeric(values) => values
            .get(row_idx)
            .copied()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::OptionalNumeric(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::Integer(values) => values
            .get(row_idx)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::Logical(values) => values
            .get(row_idx)
            .map(|value| {
                if *value {
                    "TRUE".to_string()
                } else {
                    "FALSE".to_string()
                }
            })
            .unwrap_or_else(|| "NA".to_string()),
    }
}
