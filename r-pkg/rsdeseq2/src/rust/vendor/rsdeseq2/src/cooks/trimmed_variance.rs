use super::*;

pub(super) fn trimmed_cell_variance(
    normalized_counts: &RowMajorMatrix<f64>,
    model_matrix: &DesignMatrix,
) -> Result<Vec<f64>, DeseqError> {
    let groups = model_matrix_groups(model_matrix)?;
    let mut group_ids = unique_group_ids(&groups);
    group_ids.retain(|group| {
        groups
            .iter()
            .filter(|candidate| *candidate == group)
            .count()
            >= 3
    });

    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let mut group_variances = Vec::with_capacity(group_ids.len());
        for group in &group_ids {
            let values = groups
                .iter()
                .enumerate()
                .filter_map(|(sample, candidate)| (candidate == group).then_some(row[sample]))
                .collect::<Vec<_>>();
            let bin = trim_bin(values.len());
            let mean = trimmed_mean(values.clone(), trim_ratio(bin), "Cook's trimmed cell mean")?;
            let sq_errors = values
                .into_iter()
                .enumerate()
                .map(|(sample, value)| {
                    let residual = value - mean;
                    checked_mul(
                        residual,
                        residual,
                        sample,
                        "Cook's trimmed cell residual square",
                    )
                })
                .collect::<Result<Vec<_>, DeseqError>>()?;
            group_variances.push(checked_mul(
                trim_scale(bin),
                trimmed_mean(
                    sq_errors,
                    trim_ratio(bin),
                    "Cook's trimmed cell variance mean",
                )?,
                gene,
                "Cook's trimmed cell variance",
            )?);
        }
        variances.push(row_max(&group_variances));
    }
    Ok(variances)
}

pub(super) fn trimmed_variance(
    normalized_counts: &RowMajorMatrix<f64>,
) -> Result<Vec<f64>, DeseqError> {
    let mut variances = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let row = normalized_counts.row(gene)?;
        let mean = trimmed_mean(row.to_vec(), 1.0 / 8.0, "Cook's trimmed mean")?;
        let sq_errors = row
            .iter()
            .copied()
            .enumerate()
            .map(|(sample, value)| {
                let residual = value - mean;
                checked_mul(residual, residual, sample, "Cook's trimmed residual square")
            })
            .collect::<Result<Vec<_>, DeseqError>>()?;
        variances.push(checked_mul(
            1.51,
            trimmed_mean(sq_errors, 1.0 / 8.0, "Cook's trimmed variance mean")?,
            gene,
            "Cook's trimmed variance",
        )?);
    }
    Ok(variances)
}
