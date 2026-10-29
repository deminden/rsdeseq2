//! Borrowed inference inputs; unrelated fitting scratch is never reconstructed.
use super::NbinomGlmFit;
use crate::design::DesignMatrix;
use crate::matrix::RowMajorMatrix;

#[derive(Clone, Copy)]
pub(crate) struct GlmInferenceView<'a> {
    pub beta: &'a RowMajorMatrix<f64>,
    pub beta_se: &'a RowMajorMatrix<f64>,
    pub beta_covariance: Option<&'a RowMajorMatrix<f64>>,
    pub beta_converged: &'a [bool],
    pub model_matrix: &'a DesignMatrix,
}

impl<'a> From<&'a NbinomGlmFit> for GlmInferenceView<'a> {
    fn from(fit: &'a NbinomGlmFit) -> Self {
        Self {
            beta: &fit.beta,
            beta_se: &fit.beta_se,
            beta_covariance: fit.beta_covariance.as_ref(),
            beta_converged: &fit.beta_converged,
            model_matrix: &fit.model_matrix,
        }
    }
}

impl<'a> GlmInferenceView<'a> {
    pub(crate) fn estimates(
        &self,
    ) -> Result<super::InferenceEstimates<'a>, crate::errors::DeseqError> {
        let covariance = self
            .beta_covariance
            .map(|matrix| super::CovarianceView::new(self.beta.n_cols(), matrix.into()))
            .transpose()?;
        super::InferenceEstimates::new(self.beta.into(), Some(self.beta_se.into()), covariance)
    }
}
