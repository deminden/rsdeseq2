//! Optional deterministic spectral approximation of the augmented likelihood.
use super::{DMatrix, Problem, invalid};
use crate::DeseqError;
impl Problem {
    pub(super) fn truncated(&self, tolerance: f64) -> Result<Option<Self>, DeseqError> {
        let k = self.components;
        if tolerance == 0.0 || k <= 4 {
            return Ok(None);
        }
        let rows = self
            .prior_counts
            .iter()
            .filter(|&&x| x > 0.0)
            .count()
            .checked_add(
                self.observation_weights
                    .iter()
                    .filter(|&&x| x > 0.0)
                    .count(),
            )
            .ok_or_else(|| invalid("augmented likelihood row count overflows"))?;
        let cells = rows
            .checked_mul(k)
            .ok_or_else(|| invalid("augmented likelihood size overflows"))?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(cells)
            .map_err(|_| invalid("unable to reserve augmented likelihood storage"))?;
        let mut weights = Vec::with_capacity(rows);
        for (j, &count) in self.prior_counts.iter().enumerate() {
            if count > 0.0 {
                let start = values.len();
                values.resize(start + k, 0.0);
                values[start + j] = 1.0;
                weights.push(count);
            }
        }
        for (row, &weight) in self
            .likelihoods
            .chunks_exact(k)
            .zip(&self.observation_weights)
        {
            if weight > 0.0 {
                values.extend_from_slice(row);
                weights.push(weight);
            }
        }
        let rows = weights.len();
        let dimension = rows.min(k);
        if dimension <= 2 {
            return Ok(None);
        }
        let matrix = DMatrix::from_row_slice(rows, k, &values);
        let svd = nalgebra::linalg::SVD::try_new(matrix, true, true, f64::EPSILON, 1000)
            .ok_or_else(|| invalid("likelihood SVD did not converge"))?;
        // Original rank discovery doubles its requested singular subspace.
        // Its partial solver warns above half the smaller matrix dimension,
        // and mixsqp then falls back to the full matrix. Preserve that policy
        // using deterministic singular values instead of a randomized solver.
        let mut probe = 2;
        loop {
            if probe > dimension / 2 {
                return Ok(None);
            }
            if svd.singular_values[probe - 1] < tolerance {
                break;
            }
            probe = probe.saturating_mul(2).min(dimension);
        }
        let rank = svd
            .singular_values
            .iter()
            .filter(|&&d| d > tolerance)
            .count()
            .max(2);
        if rank >= dimension {
            return Ok(None);
        }
        let u = svd
            .u
            .ok_or_else(|| invalid("missing left likelihood singular vectors"))?;
        let vt = svd
            .v_t
            .ok_or_else(|| invalid("missing right likelihood singular vectors"))?;
        let mut minimum = 0.0_f64;
        for i in 0..rows {
            for j in 0..k {
                let value = (0..rank)
                    .map(|r| u[(i, r)] * svd.singular_values[r] * vt[(r, j)])
                    .sum::<f64>();
                if !value.is_finite() {
                    return Err(invalid("nonfinite truncated likelihood"));
                }
                values[i * k + j] = value;
                minimum = minimum.min(value);
            }
        }
        Ok(Some(Self {
            likelihoods: values,
            row_maxima: vec![0.0; rows],
            observation_weights: weights,
            prior_counts: vec![0.0; k],
            mass: self.mass,
            components: k,
            offset: self.offset - minimum,
            low_rank: true,
        }))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rank_two_likelihood_is_reconstructed_without_changing_its_values() {
        // A_ij = (i+1) + (j+1)/10 has rank two by construction.
        let values = (0..12)
            .flat_map(|i| (0..10).map(move |j| (i + 1) as f64 + (j + 1) as f64 / 10.0))
            .collect::<Vec<_>>();
        let p = Problem {
            likelihoods: values.clone(),
            row_maxima: vec![0.0; 12],
            observation_weights: vec![1.0; 12],
            prior_counts: vec![0.0; 10],
            mass: 12.0,
            components: 10,
            offset: 0.0,
            low_rank: false,
        };
        let q = p.truncated(1e-8).unwrap().unwrap();
        assert!(
            q.likelihoods
                .iter()
                .zip(&values)
                .all(|(a, b)| (a - b).abs() < 1e-12)
        );
        assert!(p.truncated(0.0).unwrap().is_none());
    }
}
