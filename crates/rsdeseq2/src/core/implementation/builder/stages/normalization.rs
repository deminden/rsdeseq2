//! Builder normalization stages.
use super::*;

impl DeseqBuilder {
    /// Run only the implemented initial normalization stages.
    pub fn fit_size_factors_and_base_means(
        &self,
        counts: &CountMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_size_factors_and_base_means(counts));
        }

        let stages = self.normalization_stages(counts)?;
        Ok(self.base_fit(counts, None, stages.into_base_fit_input()))
    }

    /// Run initial normalization stages with design-aware observation-weight checks.
    ///
    /// This is useful for parity checks against DESeq2's early metadata when a
    /// `weights` assay is present: raw weights are used for `baseMean` and
    /// `baseVar`, weights are row-normalized for fitting checks, design/rank
    /// failures are recorded in `weights_fail`, and failed rows are marked in
    /// `all_zero` for downstream skipping.
    pub fn fit_size_factors_and_base_means_with_design(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_size_factors_and_base_means_with_design(counts, design));
        }

        let stages = self.normalization_stages_for_design(counts, design)?;
        Ok(self.base_fit(counts, Some(design.clone()), stages.into_base_fit_input()))
    }
}
