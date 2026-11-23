//! Builder normalization stages.
use super::*;

impl DeseqBuilder {
    pub(super) fn leaf_size_factors_and_base_means(
        &self,
        counts: &CountMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.in_pool(|| {
            let stages = self.normalization_stages(counts)?;
            Ok(self.base_fit(counts, None, stages.into_base_fit_input()))
        })
    }

    pub(super) fn leaf_size_factors_and_base_means_with_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        self.in_pool(|| {
            let stages = self.normalization_stages_for_design(counts, design)?;
            Ok(self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input()))
        })
    }
}
